"""
run_compare.py — 다중 케이스 비교 해석

여러 JSON 설정 파일을 받아 Gap에 따른 A/B/C 모듈 결과를 오버레이 비교

사용법:
  python run_compare.py case1.json case2.json case3.json
  python run_compare.py *.json

출력:
  compare_a.png  — Module A 비교 (Cap, Q_net, SHR, 손실, dT, 제습)
  compare_b.png  — Module B 비교 (Flux, P_reach, Q_carry)
  compare_c.png  — Module C 비교 (dP_total, dP 분해, Gap 비율)
  compare_results.csv — 전 케이스 통합 비교 데이터
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os, sys, json, csv, warnings
warnings.filterwarnings("ignore")

from common import (
    FinTubeSpec, MCHXSpec, RefrigerantState, DARK, C,
    compute_ft_geometry, compute_mchx_geometry, compute_UA,
    FLUX_CAUTION, FLUX_DANGER, FLUX_SEVERE,
)
from module_a import GapParams, simulate_gap, sweep
from module_b import (analyze_combined, compute_carry_penalty, apply_carry_penalty,
                      CarryoverSpec, we_ch, we_crit, v_onset, monte_carlo, eta_co)
from module_c import InletCondition, sweep_dp
from run_single import load_config, build_spec

_NOTO = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
if os.path.exists(_NOTO):
    from matplotlib import font_manager as fm
    fm.fontManager.addfont(_NOTO)
    matplotlib.rcParams['font.family'] = fm.FontProperties(fname=_NOTO).get_name()
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['mathtext.fontset'] = 'cm'

# 케이스별 색상 (최대 8개)
CASE_COLORS = ['#0077cc', '#d63031', '#00875a', '#8e44ad',
               '#e67e22', '#e17055', '#2d3436', '#6c5ce7']
CASE_LS     = ['-', '--', '-.', ':', '-', '--', '-.', ':']


# ═══════════════════════════════════════════════════════════════
#  단일 케이스 해석 (run_single 핵심 로직 재사용)
# ═══════════════════════════════════════════════════════════════

def analyze_case(cfg, label=None):
    """단일 케이스 해석 → 결과 dict 반환"""

    DEF_EVAP = dict(W=0.30, H=0.25, D=0.045, fin_pitch=1.8, fin_thickness=0.10,
        fin_height=25.4, tube_do=9.52, tube_di=8.52, tube_rows=2, tube_cols=8,
        tube_pitch_t=25.4, tube_pitch_l=22.0, fin_type="plain", corr_mode="etype",
        tube_layout="staggered", wavy_height=1.5, wavy_angle=17.5,
        ft_louver_pitch=1.7, ft_louver_angle=28.0, slit_num=6, slit_height=1.0,
        fin_material="Al", tube_material="Cu", h_drain=25.0)
    DEF_COND = dict(W=0.35, H=0.28, D=0.050, fin_pitch=1.5, fin_thickness=0.10,
        fin_height=25.4, tube_do=9.52, tube_di=8.52, tube_rows=2, tube_cols=9,
        tube_pitch_t=25.4, tube_pitch_l=22.0, fin_type="plain", corr_mode="etype",
        tube_layout="staggered", wavy_height=1.5, wavy_angle=17.5,
        ft_louver_pitch=1.7, ft_louver_angle=28.0, slit_num=6, slit_height=1.0,
        fin_material="Al", tube_material="Cu", h_drain=25.0)

    if cfg is None:
        cfg = {}
    evap_d = {**DEF_EVAP, **cfg.get('evap', {})}
    cond_d = {**DEF_COND, **cfg.get('cond', {})}
    ref_d  = {**dict(refrigerant="R410A", T_sat_evap=5.0, T_sat_cond=45.0), **cfg.get('ref', {})}
    gap_d  = {**dict(T_amb=25.0, RH_in=0.50, CMM=6.75, gap_mode="semi",
                     seal_fraction=0.7, mode="forced",
                     frame_material="Al", A_frame=4.0), **cfg.get('gap', {})}
    sim_d  = {**dict(gap_min=1, gap_max=100, gap_points=25), **cfg.get('sim', {})}

    evap = build_spec(evap_d)
    cond = build_spec(cond_d)
    ref  = RefrigerantState(refrigerant=ref_d['refrigerant'],
                            T_sat_evap=ref_d['T_sat_evap'],
                            T_sat_cond=ref_d['T_sat_cond'])
    gp = GapParams(T_amb=gap_d['T_amb'], RH_in=gap_d['RH_in'], CMM=gap_d['CMM'],
                   gap_mode=gap_d['gap_mode'], seal_fraction=gap_d['seal_fraction'],
                   mode=gap_d['mode'], frame_material=gap_d['frame_material'],
                   A_frame=gap_d['A_frame']*1e-4)
    gaps = np.linspace(sim_d['gap_min'], sim_d['gap_max'], int(sim_d['gap_points']))

    geo_e = compute_ft_geometry(evap)
    geo_c = compute_ft_geometry(cond)
    V_face = gp.CMM / (60.0 * evap.W * evap.H)
    ua_e = compute_UA(evap, geo_e, ref, V_face, 'evap')
    ua_c = compute_UA(cond, geo_c, ref, V_face, 'cond')

    if label is None:
        label = (f"{evap.fin_type.capitalize()} {evap.tube_layout[0].upper()} | "
                 f"CMM={gp.CMM} | {gp.gap_mode} sf={gap_d['seal_fraction']}")

    # Module A
    results_a = sweep(gaps, evap, cond, geo_e, geo_c, ua_e, ua_c, ref, gp)

    # Module B
    case_b = dict(name=label, T_in=gap_d['T_amb'], RH_in=gap_d['RH_in'],
                  CMM=gap_d['CMM'], T_wall=ref.T_sat_evap,
                  gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
                  gap_mode=gap_d['gap_mode'], seal_fraction=gap_d['seal_fraction'],
                  hx_type="FT", tube_layout=evap.tube_layout)
    result_b = analyze_combined(case_b, gaps, evap, geo_e, cond, geo_c, ua_e, ua_c, ref)

    # A+B 통합
    carry_penalty = compute_carry_penalty(result_b, evap, geo_e, gaps, ref.T_sat_cond)
    apply_carry_penalty(results_a, carry_penalty)

    # Module C
    inlet = InletCondition(T_in=gap_d['T_amb'], RH_in=gap_d['RH_in'],
                           CMM=gap_d['CMM'], T_wall_evap=ref.T_sat_evap)
    result_c = sweep_dp(gaps, evap, cond, geo_e, geo_c, inlet)

    # Module B 추가 데이터
    co_spec = CarryoverSpec(evap, geo_e)
    V_onset_val = v_onset(co_spec)
    We_ch_val = we_ch(V_face, co_spec)
    We_crit_val = we_crit(co_spec)
    eta_val = eta_co(V_face, co_spec, 5e-4) if V_face > V_onset_val else 0
    P_reach_arr = np.array([
        monte_carlo(g, V_face, co_spec, 0, 0.75, N=40, seed=42)['P_reach']
        for g in gaps])

    return dict(
        label=label, gaps=gaps, cfg=cfg, gap_d=gap_d,
        evap=evap, cond=cond, ref=ref, gp=gp,
        ua_e=ua_e, ua_c=ua_c, V_face=V_face,
        results_a=results_a,
        result_b=result_b,
        result_c=result_c,
        carry_penalty=carry_penalty,
        P_reach_arr=P_reach_arr,
        V_onset=V_onset_val,
        q_cond=result_b['cr']['q_cond'],
        We_ratio=We_ch_val / (We_crit_val + 1e-9),
        eta_co=eta_val,
        co_spec=co_spec,
    )


# ═══════════════════════════════════════════════════════════════
#  비교 가시화
# ═══════════════════════════════════════════════════════════════

def _styl(ax):
    ax.set_facecolor(DARK["panel"])
    ax.tick_params(colors=DARK["text"], labelsize=9)
    for s in ax.spines.values(): s.set_color(DARK["border"])
    ax.grid(True, color=DARK["grid"], alpha=0.4, lw=0.5)


def _legend(ax, ncol=1):
    ax.legend(fontsize=8, loc='best', framealpha=0.5, ncol=ncol,
              labelcolor=DARK["text"], facecolor=DARK["panel"],
              edgecolor=DARK["border"])


def compare_fig_a(cases):
    """
    Module A 비교 — 2×3
    (0,0) Cap Retention + Cap*
    (0,1) Q 분해 (Q_sensible / Q_latent / Q_net)
    (0,2) 손실 분해 (Q_rad + Q_cond + Q_mix + Q_carry)
    (1,0) ΔT_recir
    (1,1) SHR
    (1,2) 제습량
    """
    fig, axes = plt.subplots(2, 3, figsize=(22, 11), facecolor=DARK["bg"])
    fig.subplots_adjust(hspace=0.34, wspace=0.32, left=0.06, right=0.97, top=0.89, bottom=0.07)

    for ci, case in enumerate(cases):
        cc = CASE_COLORS[ci % len(CASE_COLORS)]
        ls = CASE_LS[ci % len(CASE_LS)]
        gaps = case['gaps']
        ra = case['results_a']
        lbl = case['label']

        cap   = np.array([r['cap_ratio'] for r in ra])
        cap_c = np.array([r.get('cap_corrected', r['cap_ratio']) for r in ra])
        Qnet  = np.array([r['Q_net'] for r in ra])
        Qsen  = np.array([r['Q_sensible'] for r in ra])
        Qlat  = np.array([r['Q_latent'] for r in ra])
        Qrad  = np.array([r['Q_rad'] for r in ra])
        Qcond = np.array([r['Q_cond'] for r in ra])
        Qmix  = np.array([r['Q_mix'] for r in ra])
        Qcarry = np.array([r.get('Q_carry_penalty', 0) for r in ra])
        Qtot_loss = Qrad + Qcond + Qmix + Qcarry
        dTr   = np.array([r['dT_recir'] for r in ra])
        SHR   = np.array([r['SHR'] for r in ra])
        deh   = np.array([r['dehumid_rate'] for r in ra])

        # (0,0) Cap + Cap*
        ax = axes[0,0]; _styl(ax)
        ax.plot(gaps, cap, color=cc, lw=2, ls=ls, label=lbl)
        if np.max(np.abs(cap - cap_c)) > 0.1:
            ax.plot(gaps, cap_c, color=cc, lw=1.2, ls=':', alpha=0.6)

        # (0,1) Q_sensible + Q_latent + Q_net
        ax = axes[0,1]; _styl(ax)
        ax.plot(gaps, Qnet, color=cc, lw=2, ls=ls, label=f'{lbl} (net)')
        ax.plot(gaps, Qsen+Qlat, color=cc, lw=1, ls=':', alpha=0.4)

        # (0,2) Total loss
        ax = axes[0,2]; _styl(ax)
        ax.plot(gaps, Qtot_loss, color=cc, lw=2, ls=ls, label=lbl)

        # (1,0) dT_recir
        ax = axes[1,0]; _styl(ax)
        ax.plot(gaps, dTr, color=cc, lw=2, ls=ls, label=lbl)

        # (1,1) SHR
        ax = axes[1,1]; _styl(ax)
        ax.plot(gaps, SHR, color=cc, lw=2, ls=ls, label=lbl)

        # (1,2) 제습량
        ax = axes[1,2]; _styl(ax)
        ax.plot(gaps, deh, color=cc, lw=2, ls=ls, label=lbl)

    # 제목/라벨
    axes[0,0].axhline(90, color=DARK["dim"], ls=':', alpha=0.5)
    axes[0,0].text(cases[0]['gaps'][-1]*0.95, 91, '90%', fontsize=8,
                   color=DARK["dim"], ha='right')
    titles = ['Capacity Retention', 'Heat Transfer (net / total)',
              'Total Loss (rad+cond+mix+carry)',
              'Recirculation Temperature Rise', 'Sensible Heat Ratio',
              'Dehumidification Rate']
    ylabels = ['Cap Retention [%]', 'Q [W]', '$Q_{loss}$ [W]',
               r'$\Delta T_{recir}$ [°C]', 'SHR [-]', 'Dehumidification [g/h]']
    axes[1,1].set_ylim(0, 1.05)
    axes[1,1].axhline(0.75, color=DARK["dim"], ls=':', alpha=0.3)

    ncol = 1 if len(cases) <= 3 else 2
    for i, ax in enumerate(axes.flat):
        ax.set_xlabel('Gap [mm]', color=DARK["text"])
        ax.set_title(titles[i], color=DARK["text"], fontsize=11, fontweight='bold')
        ax.set_ylabel(ylabels[i], color=DARK["text"], fontsize=10)
        _legend(ax, ncol)

    fig.suptitle(f'Module A — Cooling Performance Comparison  ({len(cases)} cases)',
                 fontsize=14, color=DARK["text"], fontweight='bold', y=0.97)
    fig.savefig("compare_a.png", dpi=150, bbox_inches='tight', facecolor=DARK["bg"])
    plt.close(fig)
    print("  → compare_a.png")


def compare_fig_b(cases):
    """
    Module B 비교 — 2×3
    (0,0) Cap vs Flux (dual axis)
    (0,1) Carryover Flux vs Gap
    (0,2) P_reach vs Gap
    (1,0) Weber ratio vs V_face (각 케이스 운전점 마커)
    (1,1) Risk Zone (grouped bar)
    (1,2) Q_carry penalty
    """
    fig, axes = plt.subplots(2, 3, figsize=(22, 11), facecolor=DARK["bg"])
    fig.subplots_adjust(hspace=0.34, wspace=0.35, left=0.06, right=0.97, top=0.89, bottom=0.07)

    # 4단계: FLUX_CAUTION=1.0 / FLUX_DANGER=4.0 / FLUX_SEVERE=12.0 (from common)

    # (1,0) Weber sweep 준비 — 공유 x축
    V_max = max(c['V_face'] for c in cases) * 2.5
    V_sweep = np.linspace(0.3, max(V_max, 6.0), 60)

    for ci, case in enumerate(cases):
        cc = CASE_COLORS[ci % len(CASE_COLORS)]
        ls = CASE_LS[ci % len(CASE_LS)]
        gaps = case['gaps']
        lbl = case['label']
        flux = case['result_b']['flux_arr']
        cap_arr = case['result_b']['cap_arr']
        P_r = case['P_reach_arr']
        cp = case['carry_penalty']

        # (0,0) Cap vs Flux dual axis
        ax = axes[0,0]; _styl(ax)
        ax.plot(gaps, cap_arr, color=cc, lw=2, ls=ls, label=f'{lbl}')

        # (0,1) Flux
        ax = axes[0,1]; _styl(ax)
        ax.plot(gaps, flux, color=cc, lw=2, ls=ls, label=lbl)

        # (0,2) P_reach
        ax = axes[0,2]; _styl(ax)
        ax.plot(gaps, P_r, color=cc, lw=2, ls=ls, label=lbl)

        # (1,0) Weber sweep
        ax = axes[1,0]; _styl(ax)
        co_sp = case['co_spec']
        We_crit_val = we_crit(co_sp)
        We_sw = np.array([we_ch(v, co_sp) / We_crit_val for v in V_sweep])
        ax.plot(V_sweep, We_sw, color=cc, lw=1.5, ls=ls, alpha=0.7)
        # 운전점 마커
        ax.scatter([case['V_face']], [case['We_ratio']], s=120, c=cc,
                   marker='*', zorder=5, edgecolors=DARK["text"], lw=1)
        ax.annotate(f'{case["V_face"]:.1f}',
                    (case['V_face'], case['We_ratio']),
                    fontsize=7, color=cc, xytext=(4, 4), textcoords='offset points')

        # (1,1) Risk Zone — 루프 내에서는 아무것도 안 함 (루프 밖에서 한번에 그림)
        pass

        # (1,2) Q_carry
        ax = axes[1,2]; _styl(ax)
        ax.plot(gaps, cp, color=cc, lw=2, ls=ls, label=lbl)

    # ── (1,1) Risk Zone Map — 케이스별 행 분리 ──
    ax = axes[1,1]; _styl(ax)
    n_cases = len(cases)
    gaps0 = cases[0]['gaps']
    bw = (gaps0[-1] - gaps0[0]) / len(gaps0) * 0.92

    for ci, case in enumerate(cases):
        gaps = case['gaps']
        flux = case['result_b']['flux_arr']
        cc = CASE_COLORS[ci % len(CASE_COLORS)]
        y_base = n_cases - 1 - ci  # 위에서부터 아래로

        for gi, (g, fl) in enumerate(zip(gaps, flux)):
            if fl < FLUX_CAUTION:   fc = '#00875a'
            elif fl < FLUX_DANGER:  fc = '#e67e22'
            elif fl < FLUX_SEVERE:  fc = '#d63031'
            else:                   fc = '#7b0000'
            ax.barh(y_base, bw, left=g - bw/2, height=0.7, color=fc, alpha=0.65,
                    edgecolor=DARK["panel"], lw=0.3)

        # 케이스 라벨 (좌측)
        short_lbl = case['label'][:20]
        ax.text(gaps[0] - bw, y_base, short_lbl, ha='right', va='center',
                fontsize=7, color=cc, fontweight='bold')

    ax.set_ylim(-0.5, n_cases - 0.5)
    ax.set_yticks([])
    ax.set_xlim(gaps0[0] - bw*2, gaps0[-1] + bw)

    # 스타일 보정
    axes[0,0].axhline(90, color=DARK["dim"], ls=':', alpha=0.5)
    axes[0,0].set_ylabel('Cap Retention [%]', color=DARK["text"])
    axes[0,0].set_title('Capacity (from Module B)', color=DARK["text"], fontsize=11, fontweight='bold')

    axes[0,1].axhline(FLUX_CAUTION, color=C[3], ls=':', alpha=0.3)
    axes[0,1].axhline(FLUX_DANGER, color=C[2], ls=':', alpha=0.4)
    axes[0,1].axhline(FLUX_SEVERE, color=C[1], ls=':', alpha=0.3)
    axes[0,1].set_ylabel('Carryover Flux [-]', color=DARK["text"])
    axes[0,1].set_title('Carryover Flux', color=DARK["text"], fontsize=11, fontweight='bold')

    axes[0,2].set_ylabel('$P_{reach}$ [-]', color=DARK["text"])
    axes[0,2].set_ylim(0, 1.05)
    axes[0,2].set_title('Droplet Reach Probability', color=DARK["text"], fontsize=11, fontweight='bold')

    axes[1,0].axhline(1.0, color=C[2], ls='--', lw=1.5, alpha=0.5)
    axes[1,0].text(V_sweep[-1]*0.95, 1.08, '$We/We_{crit}=1$',
                   fontsize=8, color=C[2], ha='right')
    axes[1,0].fill_between(V_sweep, 0, 1, alpha=0.04, color=C[3])
    axes[1,0].set_ylabel('$We / We_{crit}$ [-]', color=DARK["text"])
    axes[1,0].set_xlabel('$V_{face}$ [m/s]', color=DARK["text"])
    axes[1,0].set_title('Weber Number Ratio', color=DARK["text"], fontsize=11, fontweight='bold')

    axes[1,1].set_yticks([])
    axes[1,1].set_title('Risk Zone Map', color=DARK["text"], fontsize=11, fontweight='bold')
    from matplotlib.patches import Patch
    axes[1,1].legend(handles=[
        Patch(color='#00875a', label='SAFE'),
        Patch(color='#e67e22', label='CAUTION'),
        Patch(color='#d63031', label='DANGER'),
        Patch(color='#7b0000', label='SEVERE')],
        fontsize=8, loc='upper right', framealpha=0.5,
        labelcolor=DARK["text"], facecolor=DARK["panel"])

    axes[1,2].set_ylabel('$Q_{carry}$ [W]', color=DARK["text"])
    axes[1,2].set_title('Condenser Carry Penalty', color=DARK["text"], fontsize=11, fontweight='bold')

    ncol = 1 if len(cases) <= 3 else 2
    for ax in [axes[0,0], axes[0,1], axes[0,2], axes[1,2]]:
        ax.set_xlabel('Gap [mm]', color=DARK["text"])
        _legend(ax, ncol)
    axes[1,1].set_xlabel('Gap [mm]', color=DARK["text"])

    fig.suptitle(f'Module B — Carryover Risk Comparison  ({len(cases)} cases)',
                 fontsize=14, color=DARK["text"], fontweight='bold', y=0.97)
    fig.savefig("compare_b.png", dpi=150, bbox_inches='tight', facecolor=DARK["bg"])
    plt.close(fig)
    print("  → compare_b.png")


def compare_fig_c(cases):
    """
    Module C 비교 — 2×3
    (0,0) ΔP_total vs Gap
    (0,1) ΔP_gap vs Gap
    (0,2) ΔP_gap / ΔP_total [%]
    (1,0) ΔP stacked bar @ 대표 Gap (20mm)
    (1,1) f_dry / f_wet 비교 (grouped bar)
    (1,2) 수치 요약 테이블
    """
    fig, axes = plt.subplots(2, 3, figsize=(22, 11), facecolor=DARK["bg"])
    fig.subplots_adjust(hspace=0.38, wspace=0.32, left=0.06, right=0.97, top=0.89, bottom=0.07)

    # 라인 플롯 (0,0~0,2)
    for ci, case in enumerate(cases):
        cc = CASE_COLORS[ci % len(CASE_COLORS)]
        ls = CASE_LS[ci % len(CASE_LS)]
        gaps = case['gaps']
        lbl = case['label']
        rc = case['result_c']
        frac = rc['dp_gap'] / (rc['dp_total'] + 1e-9) * 100

        ax = axes[0,0]; _styl(ax)
        ax.plot(gaps, rc['dp_total'], color=cc, lw=2, ls=ls, label=lbl)
        ax = axes[0,1]; _styl(ax)
        ax.plot(gaps, rc['dp_gap'], color=cc, lw=2, ls=ls, label=lbl)
        ax = axes[0,2]; _styl(ax)
        ax.plot(gaps, frac, color=cc, lw=2, ls=ls, label=lbl)

    axes[0,0].set_ylabel(r'$\Delta P_{total}$ [Pa]', color=DARK["text"])
    axes[0,0].set_title('System Pressure Drop', color=DARK["text"], fontsize=11, fontweight='bold')
    axes[0,1].set_ylabel(r'$\Delta P_{gap}$ [Pa]', color=DARK["text"])
    axes[0,1].set_title('Gap Pressure Drop', color=DARK["text"], fontsize=11, fontweight='bold')
    axes[0,2].set_ylabel(r'$\Delta P_{gap}$ / Total [%]', color=DARK["text"])
    axes[0,2].set_title('Gap ΔP Fraction', color=DARK["text"], fontsize=11, fontweight='bold')

    ncol = 1 if len(cases) <= 3 else 2
    for ax in axes[0,:]:
        ax.set_xlabel('Gap [mm]', color=DARK["text"])
        _legend(ax, ncol)

    # (1,0) ΔP stacked bar @ G=20mm
    ax = axes[1,0]; _styl(ax)
    n = len(cases)
    x_pos = np.arange(n)
    bar_w = 0.5
    dp_e_vals, dp_g_vals, dp_c_vals, labels = [], [], [], []
    for case in cases:
        rc = case['result_c']
        gaps = case['gaps']
        idx20 = np.argmin(np.abs(gaps - 20))
        dp_e_vals.append(rc['dp_evap'][idx20])
        dp_g_vals.append(rc['dp_gap'][idx20])
        dp_c_vals.append(rc['dp_cond'][idx20])
        labels.append(case['label'][:18])

    dp_e = np.array(dp_e_vals)
    dp_g = np.array(dp_g_vals)
    dp_c = np.array(dp_c_vals)
    ax.bar(x_pos, dp_e, bar_w, color=C[0], alpha=0.7, label=r'$\Delta P_{evap}$')
    ax.bar(x_pos, dp_g, bar_w, bottom=dp_e, color=C[2], alpha=0.7, label=r'$\Delta P_{gap}$')
    ax.bar(x_pos, dp_c, bar_w, bottom=dp_e+dp_g, color=C[1], alpha=0.7, label=r'$\Delta P_{cond}$')
    for i, total in enumerate(dp_e+dp_g+dp_c):
        ax.text(i, total+0.5, f'{total:.0f}', ha='center', fontsize=8,
                color=DARK["text"], fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=7, rotation=15, ha='right', color=DARK["text"])
    ax.set_ylabel(r'$\Delta P$ [Pa]', color=DARK["text"])
    ax.set_title(r'$\Delta P$ Breakdown @ G=20mm', color=DARK["text"], fontsize=11, fontweight='bold')
    _legend(ax)

    # (1,1) f_dry / f_wet grouped bar
    ax = axes[1,1]; _styl(ax)
    bar_w2 = 0.35
    f_dry_vals = [c['result_c']['f_evap_dry'] for c in cases]
    f_wet_vals = [c['result_c']['f_evap_wet'] for c in cases]
    ax.bar(x_pos - bar_w2/2, f_dry_vals, bar_w2, color=C[0], alpha=0.7, label='$f_{dry}$')
    ax.bar(x_pos + bar_w2/2, f_wet_vals, bar_w2, color=C[1], alpha=0.7, label='$f_{wet}$')
    for i in range(n):
        ratio = f_wet_vals[i] / (f_dry_vals[i] + 1e-9)
        ax.text(i, max(f_wet_vals[i], f_dry_vals[i]) + 0.003,
                f'{ratio:.2f}×', ha='center', fontsize=8, color=C[2], fontweight='bold')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(labels, fontsize=7, rotation=15, ha='right', color=DARK["text"])
    ax.set_ylabel('Fanning f [-]', color=DARK["text"])
    ax.set_title('Friction Factor: Dry vs Wet', color=DARK["text"], fontsize=11, fontweight='bold')
    _legend(ax)

    # (1,2) 수치 요약 테이블
    ax = axes[1,2]; _styl(ax)
    ax.axis('off')
    # 테이블 헤더
    col_headers = ['Case', 'V_face', 'f_dry', 'f_wet', 'f_w/f_d',
                   'ΔP@20', 'gap%']
    row_y = 0.92
    for j, h in enumerate(col_headers):
        ax.text(0.02 + j*0.14, row_y, h, transform=ax.transAxes,
                fontsize=8, color=DARK["dim"], fontweight='bold', va='top')
    row_y -= 0.08
    ax.plot([0.01, 0.99], [row_y+0.03, row_y+0.03], transform=ax.transAxes,
            color=DARK["border"], lw=0.5)

    for ci, case in enumerate(cases):
        cc = CASE_COLORS[ci % len(CASE_COLORS)]
        rc = case['result_c']
        gaps = case['gaps']
        idx20 = np.argmin(np.abs(gaps - 20))
        f_r = rc.get('f_wet_ratio', rc['f_evap_wet']/(rc['f_evap_dry']+1e-9))
        vals = [
            case['label'][:16],
            f"{case['V_face']:.2f}",
            f"{rc['f_evap_dry']:.4f}",
            f"{rc['f_evap_wet']:.4f}",
            f"{f_r:.3f}",
            f"{rc['dp_total'][idx20]:.1f}",
            f"{rc['dp_gap'][idx20]/(rc['dp_total'][idx20]+1e-9)*100:.0f}%",
        ]
        for j, v in enumerate(vals):
            ax.text(0.02 + j*0.14, row_y, v, transform=ax.transAxes,
                    fontsize=7.5, color=cc if j==0 else DARK["text"], va='top',
                    fontweight='bold' if j==0 else 'normal')
        row_y -= 0.1

    ax.set_title('Key Numbers', color=DARK["text"], fontsize=11, fontweight='bold')

    fig.suptitle(f'Module C — Pressure Drop Comparison  ({len(cases)} cases)',
                 fontsize=14, color=DARK["text"], fontweight='bold', y=0.97)
    fig.savefig("compare_c.png", dpi=150, bbox_inches='tight', facecolor=DARK["bg"])
    plt.close(fig)
    print("  → compare_c.png")


# ═══════════════════════════════════════════════════════════════
#  비교 CSV
# ═══════════════════════════════════════════════════════════════

def export_compare_csv(cases):
    """전 케이스 통합 비교 CSV"""
    with open('compare_results.csv', 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, lineterminator='\n')

        # 헤더
        w.writerow([
            'Case', 'Gap_mm',
            'Cap_%', 'Cap*_%', 'Q_net_W', 'Q_sensible_W', 'Q_latent_W',
            'SHR', 'Dehumid_g_per_h', 'dT_recir_C',
            'Q_rad_W', 'Q_cond_W', 'Q_mix_W', 'Q_carry_W',
            'Flux', 'P_reach', 'Risk',
            'dP_total_Pa', 'dP_evap_Pa', 'dP_gap_Pa', 'dP_cond_Pa', 'dP_gap_%',
        ])

        for case in cases:
            lbl = case['label']
            gaps = case['gaps']
            for i, g in enumerate(gaps):
                r = case['results_a'][i]
                flux = case['result_b']['flux_arr'][i]
                risk = 'SAFE' if flux < FLUX_CAUTION else ('CAUTION' if flux < FLUX_DANGER else ('DANGER' if flux < FLUX_SEVERE else 'SEVERE'))
                rc = case['result_c']
                dp_t = rc['dp_total'][i]
                dp_g = rc['dp_gap'][i]
                w.writerow([
                    lbl, f"{g:.1f}",
                    f"{r['cap_ratio']:.2f}", f"{r.get('cap_corrected', r['cap_ratio']):.2f}",
                    f"{r['Q_net']:.1f}", f"{r['Q_sensible']:.1f}", f"{r['Q_latent']:.1f}",
                    f"{r['SHR']:.4f}", f"{r['dehumid_rate']:.1f}", f"{r['dT_recir']:.3f}",
                    f"{r['Q_rad']:.2f}", f"{r['Q_cond']:.2f}", f"{r['Q_mix']:.2f}",
                    f"{r.get('Q_carry_penalty', 0):.2f}",
                    f"{flux:.4f}", f"{case['P_reach_arr'][i]:.4f}", risk,
                    f"{dp_t:.2f}", f"{rc['dp_evap'][i]:.2f}",
                    f"{dp_g:.2f}", f"{rc['dp_cond'][i]:.2f}",
                    f"{dp_g/(dp_t+1e-9)*100:.1f}",
                ])

        # 케이스별 설계 조건 요약
        w.writerow([])
        w.writerow(['# Case Conditions'])
        w.writerow(['Case', 'Refrigerant', 'T_evap', 'T_cond', 'T_amb', 'RH',
                     'CMM', 'gap_mode', 'seal_frac', 'fin_type', 'layout',
                     'UA_evap', 'UA_cond', 'V_face', 'V_onset',
                     'f_dry', 'f_wet', 'f_wet_ratio'])
        for case in cases:
            gd = case['gap_d']
            rc = case['result_c']
            w.writerow([
                case['label'], case['ref'].refrigerant,
                case['ref'].T_sat_evap, case['ref'].T_sat_cond,
                gd['T_amb'], gd['RH_in'], gd['CMM'],
                gd['gap_mode'], gd['seal_fraction'],
                case['evap'].fin_type, case['evap'].tube_layout,
                f"{case['ua_e']['UA']:.1f}", f"{case['ua_c']['UA']:.1f}",
                f"{case['V_face']:.3f}", f"{case['V_onset']:.3f}",
                f"{rc['f_evap_dry']:.5f}", f"{rc['f_evap_wet']:.5f}",
                f"{rc.get('f_wet_ratio', 0):.4f}",
            ])

    print("  → compare_results.csv")


# ═══════════════════════════════════════════════════════════════
#  콘솔 요약
# ═══════════════════════════════════════════════════════════════

def print_summary(cases):
    g_targets = [10, 20, 50, 100]

    print(f"\n  {'':>30} |", end='')
    for g in g_targets:
        print(f"  G={g:>3}mm  ", end='')
    print()
    print(f"  {'-'*78}")

    for case in cases:
        lbl = case['label'][:28]
        gaps = case['gaps']

        # Cap
        row = f"  {lbl:>30} |"
        for g in g_targets:
            if g > gaps[-1]: row += f"{'---':>10}"; continue
            idx = np.argmin(np.abs(gaps - g))
            row += f" {case['results_a'][idx]['cap_ratio']:>6.1f}%  "
        print(row + "  Cap")

        # dP
        row = f"  {'':>30} |"
        for g in g_targets:
            if g > gaps[-1]: row += f"{'---':>10}"; continue
            idx = np.argmin(np.abs(gaps - g))
            row += f" {case['result_c']['dp_total'][idx]:>6.1f}Pa "
        print(row + "  dP")
        print()


# ═══════════════════════════════════════════════════════════════
#  메인
# ═══════════════════════════════════════════════════════════════

def main():
    if len(sys.argv) < 3:
        print("사용법: python run_compare.py case1.json case2.json [case3.json ...]")
        print("  2개 이상의 JSON 설정 파일 필요")
        sys.exit(1)

    paths = sys.argv[1:]
    print(f"\n{'='*70}")
    print(f"  Multi-Case Gap Comparison  ({len(paths)} cases)")
    print(f"{'='*70}")

    cases = []
    for i, path in enumerate(paths):
        name = os.path.splitext(os.path.basename(path))[0]
        print(f"\n  [{i+1}/{len(paths)}] {path}...")
        cfg = load_config(path)
        if cfg is None:
            print(f"    SKIP: cannot load {path}")
            continue

        # JSON에 label 필드가 있으면 사용, 없으면 파일명
        label = cfg.get('label', name)
        case = analyze_case(cfg, label)
        cases.append(case)
        print(f"    {label}: UA_e={case['ua_e']['UA']:.1f} UA_c={case['ua_c']['UA']:.1f}"
              f"  V={case['V_face']:.2f}m/s")

    if len(cases) < 2:
        print("\n  비교를 위해 2개 이상의 유효 케이스가 필요합니다.")
        sys.exit(1)

    # 가시화
    print(f"\n  Generating comparison figures...")
    compare_fig_a(cases)
    compare_fig_b(cases)
    compare_fig_c(cases)

    # CSV
    print(f"\n  Exporting CSV...")
    export_compare_csv(cases)

    # 콘솔 요약
    print_summary(cases)

    print(f"\n  → compare_a.png, compare_b.png, compare_c.png, compare_results.csv")


if __name__ == '__main__':
    main()
