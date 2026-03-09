"""visualize.py — 시각화 (모듈 A/B Figure)"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from common import *
from module_a import GapParams, simulate_gap, sweep

def _combo_schematic(ax, evap_spec, evap_geo, cond_spec, cond_geo,
                     ua_ev, ua_cd, combo_label: str):
    """단일 조합 개략도 (증발기 → 간격 → 응축기)
    공기 주 흐름: 좌(EVAP 흡입) → Gap → 우(COND 배출)
    재순환:       우(COND 출구 고온 배출) → 역방향 곡선 → 좌(EVAP 흡입구)
    """
    ax.set_xlim(0,10); ax.set_ylim(0,5.5)
    ax.axis('off')
    ax.set_title(combo_label, fontsize=8.5, color="#ffc940", pad=4)

    def _hx_block(x0, spec, ua_d, side):
        clr = "#00c8ff" if side=="evap" else "#ff5e3a"
        tag = "EVAP" if side=="evap" else "COND"
        typ = spec.hx_type
        rect = mpatches.FancyBboxPatch((x0,0.8), 2.0, 3.4,
                   boxstyle="round,pad=0.06", fc=clr+"22", ec=clr, lw=1.2)
        ax.add_patch(rect)
        for y in np.linspace(1.0, 4.0, 7):
            ax.plot([x0+0.15, x0+1.85], [y,y], color=clr, lw=0.6, alpha=0.5)
        if typ == "MCHX":
            for cx in np.linspace(x0+0.3, x0+1.7, 4):
                ax.add_patch(plt.Circle((cx,2.5), 0.18, fc=clr+"55", ec=clr, lw=0.7))
        else:
            for cx in np.linspace(x0+0.35, x0+1.65, 3):
                ax.add_patch(plt.Circle((cx,2.5), 0.28, fc=clr+"44", ec=clr, lw=0.8))
        ax.text(x0+1.0, 4.45, f"{tag}\n{typ}", ha='center', fontsize=7,
                color=clr, fontweight='bold')
        ax.text(x0+1.0, 0.38,
                f"UA={ua_d['UA']:.0f}W/K\nh_o={ua_d['h_o']:.0f}  h_i={ua_d['h_i']:.0f}",
                ha='center', fontsize=5.8, color=DARK["dim"])

    _hx_block(0.3, evap_spec, ua_ev, "evap")
    _hx_block(7.7, cond_spec, ua_cd, "cond")

    # Gap 화살표
    ax.annotate('', xy=(7.65,2.5), xytext=(2.35,2.5),
                arrowprops=dict(arrowstyle='<->', color='#ffc940', lw=1.4))
    ax.text(5.0, 2.72, "Gap G", ha='center', fontsize=7.5, color='#ffc940', fontweight='bold')

    # 공기 주 흐름 화살표 (증발기 → 응축기)
    ax.annotate('', xy=(7.55,1.55), xytext=(2.45,1.55),
                arrowprops=dict(arrowstyle='->', color='#5599cc', lw=1.3))
    ax.text(5.0, 1.35, "공기 흐름 →", ha='center', fontsize=6.0, color='#5599cc')

    # 재순환 화살표 (응축기 출구 → 역방향 → 증발기 흡입구)
    ax.annotate('', xy=(2.35,3.8), xytext=(7.65,3.8),
                arrowprops=dict(arrowstyle='->', color='#cc7a30', lw=1.2,
                                connectionstyle='arc3,rad=0.4'))
    ax.text(5.0, 5.1, "재순환", ha='center', fontsize=6.5, color='#cc7a30')

    # 복사
    xr = np.linspace(2.4, 7.6, 50)
    ax.plot(xr, 2.1 + 0.07*np.sin(xr*9), color='#cc8830', lw=1.0, ls='--', alpha=0.7)
    ax.text(5.0, 1.8, "복사", ha='center', fontsize=6, color='#cc8830')

    # 전도
    ax.plot([2.35,7.65],[1.1,1.1], color='#ff4444', lw=2.2, solid_capstyle='round')
    ax.annotate('', xy=(7.65,1.1), xytext=(5.0,1.1),
                arrowprops=dict(arrowstyle='->', color='#ff4444', lw=1.2))
    ax.text(5.0, 0.7, "전도단락", ha='center', fontsize=6, color='#ff4444')


def _refr_panel(ax, ref):
    ax.set_facecolor(DARK["panel"]); ax.axis('off')
    ax.set_title(f"냉매  {ref.refrigerant}", fontsize=9, color="#d070ff")
    for i,(k,v) in enumerate([
        ("T_evap / T_cond", f"{ref.T_sat_evap:.0f} / {ref.T_sat_cond:.0f} °C"),
        ("P_evap / P_cond",  f"{ref.P_evap/1e5:.2f} / {ref.P_cond/1e5:.2f} bar"),
        ("h_fg (evap)",      f"{ref.h_fg_evap/1e3:.1f} kJ/kg"),
        ("ρ_l / ρ_v",        f"{ref.rho_l_evap:.1f} / {ref.rho_v_evap:.2f} kg/m³"),
        ("μ_l",              f"{ref.mu_l_evap*1e6:.1f} μPa·s"),
        ("k_l",              f"{ref.k_l_evap*1e3:.2f} mW/mK"),
        ("Pr_l",             f"{ref.Pr_l_evap:.2f}"),
    ]):
        ax.text(0.04, 0.90-i*0.12, k, transform=ax.transAxes, fontsize=7, color=DARK["dim"])
        ax.text(0.97, 0.90-i*0.12, v, transform=ax.transAxes,
                fontsize=7, color="#d070ff", ha='right', fontweight='bold')


def _ua_bar(ax, ua_dict):
    labels = list(ua_dict.keys())
    comp_keys  = ["R_o","R_wall","R_i"]
    comp_names = ["공기측 R_o","벽 R_wall","냉매측 R_i"]
    comp_clrs  = ["#3a7fc1","#ffc940","#e05a3a"]
    x = np.arange(len(labels)); bottoms = np.zeros(len(labels))
    for ck, clr, cnm in zip(comp_keys, comp_clrs, comp_names):
        fracs = [ua_dict[l][ck]/(ua_dict[l]['R_o']+ua_dict[l]['R_wall']+ua_dict[l]['R_i'])*100
                 for l in labels]
        ax.bar(x, fracs, 0.55, bottom=bottoms, color=clr, label=cnm, alpha=0.85)
        for xi,(f,b) in enumerate(zip(fracs,bottoms)):
            if f > 6: ax.text(xi, b+f/2, f"{f:.0f}%", ha='center', va='center',
                              fontsize=6.5, color=DARK["text"], fontweight='bold')
        bottoms += np.array(fracs)
    for xi,l in enumerate(labels):
        ax.text(xi, 103, f"UA={ua_dict[l]['UA']:.0f}\nh_o={ua_dict[l]['h_o']:.0f}",
                ha='center', fontsize=6, color=DARK["text"])
    ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=7.5)
    ax.set_ylabel("열저항 비율 [%]", fontsize=8); ax.set_ylim(0,118)
    ax.set_title("HX 별 열저항 분포 및 UA", fontsize=9, color=DARK["text"])
    ax.legend(fontsize=7, framealpha=0.3, loc='upper right')
    ax.grid(axis='y', alpha=0.3); ax.tick_params(labelsize=7)


# ─── 통합 핵심 패널: Cap Retention + Carryover Flux (Sweet spot) ──

def _sweet_spot_panel(ax, combined_results, case_colors):
    """
    [신규] Sweet-spot 분석 패널
    Cap Retention (좌축, %) 과 Carryover Flux (우축, Mid η) 를 동시 표시
    최적 Gap 구간: Cap≥95% AND Flux≤FLUX_DANGER
    """
    ax2 = ax.twinx()
    for res, cc in zip(combined_results, case_colors):
        gaps = res['gaps']
        ax.plot(gaps, res['cap_arr'], color=cc, lw=2.0,
                label=res['case']['name'][:12])
        ax2.plot(gaps, np.maximum(res['flux_arr'], 0.01), color=cc,
                 lw=1.5, ls='--', alpha=0.7)
        if res['sweet_gap'] is not None:
            ax.axvline(res['sweet_gap'], color=cc, lw=1.2, ls=':', alpha=0.6)
            ax.text(res['sweet_gap']+0.5, 12,
                    f"G*={res['sweet_gap']:.0f}mm", fontsize=6, color=cc, rotation=90)

    # 기준선
    ax.axhline(95, color='#ffc940', lw=0.9, ls=':', alpha=0.7)
    ax2.axhline(FLUX_DANGER, color='#ff7c43', lw=0.9, ls=':', alpha=0.7)
    ax2.text(98, FLUX_DANGER*1.1, "위험기준", fontsize=6.5, color='#ff7c43', ha='right')

    ax.set_xlabel("Gap G [mm]", fontsize=8)
    ax.set_ylabel("Cap Retention [%]", fontsize=8, color='#5599cc')
    ax2.set_ylabel("Carryover Flux [mg/m²s] (Mid η)", fontsize=8, color='#cc8830')
    ax.set_ylim(0, 110); ax2.set_ylim(0, max(
        max(res['flux_arr'].max() for res in combined_results)*1.3, 15))
    ax.set_title("통합 Sweet-spot 분석 (냉방성능 ↔ 비말동반 위험도)", fontsize=9, color=DARK["text"])
    ax.legend(fontsize=6.5, framealpha=0.3, loc='lower right')
    ax.grid(True); ax.tick_params(labelsize=7); ax2.tick_params(labelsize=7)


def _flux_gap_panel(ax, combined_results, case_colors):
    """Gap vs Carryover Flux (Mid η, 3케이스)"""
    for res, cc in zip(combined_results, case_colors):
        ax.semilogy(res['gaps'], np.maximum(res['flux_arr'], 1e-3),
                    color=cc, lw=1.8, label=res['case']['name'][:12])
        if res['sweet_gap'] is not None:
            ax.axvline(res['sweet_gap'], color=cc, lw=0.8, ls=':', alpha=0.5)
    for fv, lbl, clr in [(FLUX_CAUTION,"안전",C[3]),(FLUX_DANGER,"주의",C[2]),(FLUX_SEVERE,"위험",C[1])]:
        ax.axhline(fv, color=clr, lw=1.0, ls=':', alpha=0.7)
        ax.text(98, fv*1.12, lbl, ha='right', fontsize=6.5, color=clr)
    ax.set_title("Carryover Flux vs Gap  (Mid η)", fontsize=8.5, color=DARK["text"])
    ax.set_xlabel("Gap G [mm]", fontsize=8); ax.set_ylabel("Flux [mg/m²s]", fontsize=8)
    ax.set_ylim(0.001, 500); ax.grid(True, which='both', alpha=0.3)
    ax.legend(fontsize=6.5, framealpha=0.3); ax.tick_params(labelsize=7)


# ═══════════════════════════════════════════════════════════════════
#  모듈 A Figure — 냉방성능 (Cap Retention)
# ═══════════════════════════════════════════════════════════════════

def make_module_a_figure(combos_v3, gaps_v3, ref, CMM_ref):
    """
    모듈 A 독립 Figure:
      Row 0 : 6조합 개략도
      Row 1 : Cap Retention 6조합 | 냉매 물성
      Row 2 : Q_evap | Q_net | ΔT_recir | UA breakdown
      Row 3 : 민감도 히트맵 (Gap × CMM) × 2 | UA bar
    """
    apply_style()
    fig = plt.figure(figsize=(26, 18))
    fig.patch.set_facecolor(DARK["bg"])

    gs = gridspec.GridSpec(4, 6, figure=fig,
                           hspace=0.48, wspace=0.36,
                           top=0.93, bottom=0.04, left=0.05, right=0.97)

    ax_s   = [fig.add_subplot(gs[0, i]) for i in range(6)]
    ax_cap = fig.add_subplot(gs[1, :5])
    ax_ref = fig.add_subplot(gs[1, 5])
    ax_qev = fig.add_subplot(gs[2, 0:2])
    ax_qn  = fig.add_subplot(gs[2, 2:4])
    ax_dtr = fig.add_subplot(gs[2, 4:6])
    ax_hm0 = fig.add_subplot(gs[3, 0:2])
    ax_hm1 = fig.add_subplot(gs[3, 2:4])
    ax_ua  = fig.add_subplot(gs[3, 4:6])

    # ── Row 0: 개략도 ────────────────────────────────────────
    for i, cb in enumerate(combos_v3[:6]):
        _combo_schematic(ax_s[i], cb['evap_spec'], cb['evap_geo'],
                         cb['cond_spec'], cb['cond_geo'],
                         cb['ua_ev'], cb['ua_cd'], cb['label'])

    # ── Row 1: Cap Retention ─────────────────────────────────
    for cb in combos_v3:
        cap = [r['cap_ratio'] for r in cb['results']]
        ax_cap.plot(gaps_v3, cap, color=cb['color'], lw=2.0,
                    ls=cb['ls'], label=cb['label'])
        try:
            g95 = gaps_v3[next(i for i,c in enumerate(cap) if c>=95.0)]
            ax_cap.axvline(g95, color=cb['color'], lw=0.7, alpha=0.35, ls=':')
        except StopIteration:
            pass
    ax_cap.axhline(95, color="#ffc940", lw=0.9, ls=':', alpha=0.7)
    ax_cap.text(99, 95.8, "95%", fontsize=7, color="#ffc940", ha='right')
    ax_cap.set_title(f"유효 냉방능력 유지율  (6조합 비교, CMM={CMM_ref} m³/min)",
                     fontsize=10, color=DARK["text"])
    ax_cap.set_xlabel("Gap G [mm]", fontsize=9)
    ax_cap.set_ylabel("Cap. Retention [%]", fontsize=9)
    ax_cap.set_ylim(0, 108); ax_cap.grid(True)
    ax_cap.legend(fontsize=7.5, framealpha=0.4, loc='lower right')
    ax_cap.tick_params(labelsize=7)
    _refr_panel(ax_ref, ref)

    # ── Row 2: Q_evap, Q_net, ΔT_recir ──────────────────────
    for cb in combos_v3:
        kw = dict(color=cb['color'], lw=1.8, ls=cb['ls'], label=cb['label'])
        ax_qev.plot(gaps_v3, [r['Q_evap']   for r in cb['results']], **kw)
        ax_qn.plot(gaps_v3,  [r['Q_net']    for r in cb['results']], **kw)
        ax_dtr.plot(gaps_v3, [r['dT_recir'] for r in cb['results']], **kw)
    for ax, title, ylabel in [
        (ax_qev, "총 증발기 열전달량  Q_evap", "Q [W]"),
        (ax_qn,  "순 냉방능력  Q_net",          "Q [W]"),
        (ax_dtr, "재순환 ΔT_recir",             "ΔT [K]"),
    ]:
        ax.set_title(title, fontsize=9, color=DARK["text"])
        ax.set_xlabel("Gap G [mm]", fontsize=8); ax.set_ylabel(ylabel, fontsize=8)
        ax.grid(True); ax.legend(fontsize=6, framealpha=0.3); ax.tick_params(labelsize=7)

    # ── Row 3: 민감도 히트맵 + UA ────────────────────────────
    gaps_2d = np.linspace(5, 100, 20)
    cmm_2d  = np.linspace(2, 20, 10)
    cmap_rg = plt.colormaps['RdYlGn']

    for ax_hm, cb_idx in [(ax_hm0, 0), (ax_hm1, 1)]:
        cb = combos_v3[cb_idx]
        A_f = cb['evap_spec'].W * cb['evap_spec'].H
        Z  = np.zeros((len(cmm_2d), len(gaps_2d)))
        for j, cmm_v in enumerate(cmm_2d):
            vf = cmm_v / (60.0 * A_f)
            for k, g in enumerate(gaps_2d):
                gp_tmp = GapParams(CMM=cmm_v, T_amb=25.0,
                                   frame_material=cb['gp'].frame_material)
                gp_tmp.set_V_face(A_f)
                ua_ev_t = compute_UA(cb['evap_spec'], cb['evap_geo'], ref, vf, 'evap')
                r = simulate_gap(g, cb['evap_spec'], cb['cond_spec'],
                                 cb['evap_geo'], cb['cond_geo'],
                                 ua_ev_t, cb['ua_cd'], ref, gp_tmp)
                Z[j, k] = r['cap_ratio']
        im = ax_hm.imshow(Z, origin='lower', aspect='auto', cmap=cmap_rg,
                          vmin=40, vmax=100,
                          extent=[gaps_2d[0], gaps_2d[-1], cmm_2d[0], cmm_2d[-1]])
        cb_bar = plt.colorbar(im, ax=ax_hm, pad=0.02)
        cb_bar.set_label("[%]", fontsize=7, color=DARK["text"])
        cb_bar.ax.tick_params(labelsize=6, colors=DARK["dim"])
        ax_hm.contour(gaps_2d, cmm_2d, Z, levels=[80,90,95],
                      colors=['white'], linewidths=0.7, alpha=0.7)
        ax_hm.set_title(f"민감도: {cb['label']}", fontsize=8, color=DARK["text"])
        ax_hm.set_xlabel("Gap G [mm]", fontsize=8)
        ax_hm.set_ylabel("CMM [m³/min]", fontsize=8)
        ax_hm.tick_params(labelsize=7)

    ua_dict = {}
    for cb in combos_v3:
        ua_dict[f"evap\n{cb['label'].split('/')[0]}"] = cb['ua_ev']
    for cb in combos_v3:
        ua_dict[f"cond\n{cb['label'].split('/')[1]}"] = cb['ua_cd']
    seen, ua_dedup = set(), {}
    for lbl, ua_d in ua_dict.items():
        key = round(ua_d['UA'], 1)
        if key not in seen: seen.add(key); ua_dedup[lbl] = ua_d
    _ua_bar(ax_ua, ua_dedup)

    fig.suptitle(
        "모듈 A — 냉방성능 유지율  │  Gap × 조합(Staggered/Inline/MCHX) × CMM",
        fontsize=13, color=DARK["text"], fontweight='bold', y=0.975)
    fig.text(0.97, 0.975,
             f"{ref.refrigerant}  T_sat: {ref.T_sat_evap}°C/{ref.T_sat_cond}°C"
             f"  CMM={CMM_ref} m³/min",
             fontsize=8, color=DARK["dim"], ha='right', va='top')
    return fig


# ═══════════════════════════════════════════════════════════════════
#  모듈 B Figure — 비말동반 위험도
# ═══════════════════════════════════════════════════════════════════

def make_module_b_figure(combined_results, case_colors, gaps, ref, CMM_ref):
    """
    모듈 B 독립 Figure:
      Row 0 : Carryover Flux vs Gap (전체 케이스) | Sweet-spot 통합
      Row 1 : Flux vs Gap (log) | 케이스별 q_cond bar | 케이스별 V_onset bar
    """
    apply_style()
    fig = plt.figure(figsize=(22, 12))
    fig.patch.set_facecolor(DARK["bg"])

    gs = gridspec.GridSpec(2, 4, figure=fig,
                           hspace=0.40, wspace=0.32,
                           top=0.91, bottom=0.06, left=0.06, right=0.97)

    ax_ss   = fig.add_subplot(gs[0, :2])
    ax_flux = fig.add_subplot(gs[0, 2:])
    ax_flog = fig.add_subplot(gs[1, 0:2])
    ax_qbar = fig.add_subplot(gs[1, 2])
    ax_vbar = fig.add_subplot(gs[1, 3])

    # ── Sweet-spot 분석 ──────────────────────────────────────
    _sweet_spot_panel(ax_ss, combined_results, case_colors)
    ax_ss.set_title(f"Sweet-spot 분석 (Cap Ret. ↔ Carryover Flux, CMM={CMM_ref})",
                    fontsize=10, color=DARK["text"])

    # ── Flux vs Gap (선형) ───────────────────────────────────
    for res, cc in zip(combined_results, case_colors):
        ax_flux.plot(res['gaps'], np.maximum(res['flux_arr'], 0.01),
                     color=cc, lw=1.8, label=res['case']['name'][:16])
        if res['sweet_gap'] is not None:
            ax_flux.axvline(res['sweet_gap'], color=cc, lw=0.8, ls=':', alpha=0.5)
    for fv, lbl, clr in [(FLUX_CAUTION,"안전",C[3]),(FLUX_DANGER,"주의",C[2]),(FLUX_SEVERE,"위험",C[1])]:
        ax_flux.axhline(fv, color=clr, lw=1.0, ls=':', alpha=0.7)
        ax_flux.text(98, fv*1.08, lbl, ha='right', fontsize=7, color=clr)
    ax_flux.set_title("Carryover Flux vs Gap (전체 케이스)", fontsize=10, color=DARK["text"])
    ax_flux.set_xlabel("Gap G [mm]", fontsize=9)
    ax_flux.set_ylabel("Flux [mg/m²s]", fontsize=9)
    ax_flux.legend(fontsize=6.5, framealpha=0.3, loc='upper right')
    ax_flux.grid(True, alpha=0.3); ax_flux.tick_params(labelsize=7)

    # ── Flux vs Gap (log) ────────────────────────────────────
    _flux_gap_panel(ax_flog, combined_results, case_colors)
    ax_flog.set_title("Carryover Flux vs Gap (log scale)", fontsize=10, color=DARK["text"])

    # ── q_cond bar ───────────────────────────────────────────
    names = [r['case']['name'].split('—')[1].strip()[:12] for r in combined_results]
    qconds = [r['cr']['q_cond'] for r in combined_results]
    x_q = np.arange(len(names))
    bars_q = ax_qbar.barh(x_q, qconds, 0.6, color=case_colors, alpha=0.8)
    for bar, v in zip(bars_q, qconds):
        ax_qbar.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
                     f"{v:.2f}", va='center', fontsize=7, color=DARK["text"])
    ax_qbar.set_yticks(x_q); ax_qbar.set_yticklabels(names, fontsize=7)
    ax_qbar.set_xlabel("q_cond [g/s]", fontsize=8)
    ax_qbar.set_title("응결수 생성량", fontsize=9, color=DARK["text"])
    ax_qbar.grid(axis='x', alpha=0.3); ax_qbar.tick_params(labelsize=7)

    # ── V_onset bar ──────────────────────────────────────────
    vonsets = [r['V_on'] for r in combined_results]
    bars_v = ax_vbar.barh(x_q, vonsets, 0.6, color=case_colors, alpha=0.8)
    for bar, v in zip(bars_v, vonsets):
        ax_vbar.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height()/2,
                     f"{v:.2f}", va='center', fontsize=7, color=DARK["text"])
    ax_vbar.set_yticks(x_q); ax_vbar.set_yticklabels(names, fontsize=7)
    ax_vbar.set_xlabel("V_onset [m/s]", fontsize=8)
    ax_vbar.set_title("비말 이탈 개시 풍속", fontsize=9, color=DARK["text"])
    ax_vbar.grid(axis='x', alpha=0.3); ax_vbar.tick_params(labelsize=7)

    fig.suptitle(
        "모듈 B — 비말동반 위험도  │  Carryover Flux × Sweet-spot × 케이스 비교",
        fontsize=13, color=DARK["text"], fontweight='bold', y=0.975)
    fig.text(0.97, 0.975,
             f"{ref.refrigerant}  CMM={CMM_ref} m³/min",
             fontsize=8, color=DARK["dim"], ha='right', va='top')
    return fig


# ═══════════════════════════════════════════════════════════════════
#  모듈 C Figure — 압력강하
# ═══════════════════════════════════════════════════════════════════

def make_module_c_figure(dp_results: dict, gaps, inlet, ref, CMM_ref):
    """
    모듈 C 독립 Figure (2행×4열):
      Row 0 : Gap vs ΔP_system (조합비교) | ΔP stacked | 대표 Gap bar
      Row 1 : Gap ΔP 상세 | f-factor bar | 습면 보정 vs RH | 민감도 맵
    """
    from module_c import (f_factor_ft, f_factor_mchx,
                          wet_dp_correction, InletCondition)

    apply_style()
    fig = plt.figure(figsize=(22, 14))
    fig.patch.set_facecolor(DARK["bg"])

    gs = gridspec.GridSpec(2, 4, figure=fig,
                           hspace=0.40, wspace=0.34,
                           top=0.91, bottom=0.06, left=0.06, right=0.97)

    ax_main = fig.add_subplot(gs[0, :2])
    ax_stk  = fig.add_subplot(gs[0, 2])
    ax_bar2 = fig.add_subplot(gs[0, 3])
    ax_gapd = fig.add_subplot(gs[1, 0])
    ax_fbar = fig.add_subplot(gs[1, 1])
    ax_wet  = fig.add_subplot(gs[1, 2])
    ax_hm   = fig.add_subplot(gs[1, 3])

    colors_dp = [C[0], C[3], C[4]]
    labels = list(dp_results.keys())

    # ── Row 0 Left: Gap vs ΔP_system ────────────────────────
    for (label, sw), clr in zip(dp_results.items(), colors_dp):
        ax_main.plot(gaps, sw['dp_total'], color=clr, lw=2.2, label=label)
    ax_main.set_title(f"시스템 압력강하 ΔP vs Gap  (CMM={CMM_ref} m³/min)",
                      fontsize=10, color=DARK["text"])
    ax_main.set_xlabel("Gap G [mm]", fontsize=9)
    ax_main.set_ylabel("ΔP_system [Pa]", fontsize=9)
    ax_main.legend(fontsize=8, framealpha=0.4)
    ax_main.grid(True); ax_main.tick_params(labelsize=7)

    # ── Row 0 Center: ΔP stacked (첫 번째 조합) ─────────────
    sw0 = dp_results[labels[0]]
    ax_stk.fill_between(gaps, 0, sw0['dp_evap'],
                        alpha=0.5, color='#00c8ff', label='ΔP_evap')
    ax_stk.fill_between(gaps, sw0['dp_evap'],
                        sw0['dp_evap'] + sw0['dp_gap'],
                        alpha=0.5, color='#ffc940', label='ΔP_gap')
    ax_stk.fill_between(gaps, sw0['dp_evap'] + sw0['dp_gap'],
                        sw0['dp_total'],
                        alpha=0.5, color='#ff5e3a', label='ΔP_cond')
    ax_stk.plot(gaps, sw0['dp_total'], color=DARK["text"], lw=1.5, alpha=0.8)
    ax_stk.set_title(f"ΔP 성분 분해: {labels[0]}", fontsize=9, color=DARK["text"])
    ax_stk.set_xlabel("Gap G [mm]", fontsize=8)
    ax_stk.set_ylabel("ΔP [Pa]", fontsize=8)
    ax_stk.legend(fontsize=7, framealpha=0.3)
    ax_stk.grid(True, alpha=0.3); ax_stk.tick_params(labelsize=7)

    # ── Row 0 Right: 대표 Gap grouped bar ────────────────────
    gap_pts = [5, 20, 100]
    gap_idx = [np.argmin(np.abs(gaps - g)) for g in gap_pts]
    x_grp = np.arange(len(gap_pts))
    w_bar = 0.22
    for ic, ((label, sw), clr) in enumerate(zip(dp_results.items(), colors_dp)):
        vals = [sw['dp_total'][idx] for idx in gap_idx]
        bars_g = ax_bar2.bar(x_grp + ic * w_bar, vals, w_bar,
                             color=clr, alpha=0.8, label=label)
        for bar, v in zip(bars_g, vals):
            ax_bar2.text(bar.get_x() + bar.get_width()/2,
                        bar.get_height() + 0.3,
                        f"{v:.0f}", ha='center', fontsize=6.5,
                        color=clr, fontweight='bold')
    ax_bar2.set_xticks(x_grp + w_bar)
    ax_bar2.set_xticklabels([f"G={g}mm" for g in gap_pts], fontsize=8)
    ax_bar2.set_ylabel("ΔP [Pa]", fontsize=8)
    ax_bar2.set_title("대표 Gap × 조합 ΔP 비교", fontsize=9, color=DARK["text"])
    ax_bar2.legend(fontsize=7, framealpha=0.3)
    ax_bar2.grid(axis='y', alpha=0.3); ax_bar2.tick_params(labelsize=7)

    # ── Row 1 Left: Gap ΔP 상세 비교 ────────────────────────
    for (label, sw), clr in zip(dp_results.items(), colors_dp):
        ax_gapd.plot(gaps, sw['dp_gap'], color=clr, lw=2.0, label=label)
    ax_gapd.fill_between(gaps, 0, sw0['dp_gap'], alpha=0.15, color=C[2])
    ax_gapd.set_title("Gap 공간 압력손실 비교", fontsize=9, color=DARK["text"])
    ax_gapd.set_xlabel("Gap G [mm]", fontsize=8)
    ax_gapd.set_ylabel("ΔP_gap [Pa]", fontsize=8)
    ax_gapd.legend(fontsize=7, framealpha=0.3)
    ax_gapd.grid(True, alpha=0.3); ax_gapd.tick_params(labelsize=7)

    # ── Row 1 Center-Left: f-factor bar ──────────────────────
    f_labels_b, f_vals_b, f_colors_b = [], [], []
    for (label, sw), clr in zip(dp_results.items(), colors_dp):
        f_labels_b.append(f"evap\n{label.split('/')[0]}")
        f_vals_b.append(sw['f_evap_dry'])
        f_colors_b.append(clr)
    x_f = np.arange(len(f_vals_b))
    bars = ax_fbar.bar(x_f, f_vals_b, 0.55, color=f_colors_b, alpha=0.8)
    for bar, fv in zip(bars, f_vals_b):
        ax_fbar.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.001,
                     f"{fv:.4f}", ha='center', fontsize=7, color=DARK["text"])
    ax_fbar.set_xticks(x_f); ax_fbar.set_xticklabels(f_labels_b, fontsize=7)
    ax_fbar.set_ylabel("f-factor [-]", fontsize=8)
    ax_fbar.set_title("증발기 f-factor (건면)", fontsize=9, color=DARK["text"])
    ax_fbar.grid(axis='y', alpha=0.3); ax_fbar.tick_params(labelsize=7)

    # ── Row 1 Center-Right: 습면 보정 vs RH ──────────────────
    from dataclasses import dataclass as dc
    @dc
    class _MkFT:
        fin_pitch: float = 1.8e-3; tube_do: float = 9.52e-3; hx_type: str = "FT"
    @dc
    class _MkMCHX:
        louver_pitch: float = 1.2e-3; hx_type: str = "MCHX"

    rh_range = np.linspace(0.30, 0.95, 30)
    wet_ft   = [wet_dp_correction(_MkFT(), 27.0, rh, 5.0) for rh in rh_range]
    wet_mchx = [wet_dp_correction(_MkMCHX(), 27.0, rh, 5.0) for rh in rh_range]
    wet_ft30 = [wet_dp_correction(_MkFT(), 30.0, rh, 5.0) for rh in rh_range]
    ax_wet.plot(rh_range*100, wet_ft, color=C[0], lw=2.0, label='FT 27°C')
    ax_wet.plot(rh_range*100, wet_ft30, color=C[1], lw=1.5, ls='--', label='FT 30°C')
    ax_wet.plot(rh_range*100, wet_mchx, color=C[4], lw=2.0, label='MCHX 27°C')
    ax_wet.axhline(1.0, color=DARK["dim"], lw=0.8, ls=':')
    ax_wet.axvline(70, color='#ffc940', lw=0.8, ls=':', alpha=0.5)
    ax_wet.set_title("습면 보정 f_wet/f_dry vs RH", fontsize=9, color=DARK["text"])
    ax_wet.set_xlabel("상대습도 RH [%]", fontsize=8)
    ax_wet.set_ylabel("f_wet / f_dry [-]", fontsize=8)
    ax_wet.legend(fontsize=7, framealpha=0.3)
    ax_wet.grid(True, alpha=0.3); ax_wet.tick_params(labelsize=7)

    # ── Row 1 Right: 민감도 히트맵 (Gap × CMM) ──────────────
    cmm_range = np.linspace(2, 15, 20)
    gap_range = np.linspace(5, 100, 20)
    Z = np.zeros((len(cmm_range), len(gap_range)))
    for j, cmm_v in enumerate(cmm_range):
        for k, g in enumerate(gap_range):
            idx_g = np.argmin(np.abs(gaps - g))
            V_ratio = (cmm_v / CMM_ref)**2
            Z[j, k] = sw0['dp_total'][idx_g] * V_ratio

    cmap = plt.colormaps['magma_r']
    im = ax_hm.imshow(Z, origin='lower', aspect='auto', cmap=cmap,
                       extent=[gap_range[0], gap_range[-1],
                               cmm_range[0], cmm_range[-1]])
    cb_bar = plt.colorbar(im, ax=ax_hm, pad=0.02)
    cb_bar.set_label("ΔP [Pa]", fontsize=7, color=DARK["text"])
    cb_bar.ax.tick_params(labelsize=6, colors=DARK["dim"])
    cs = ax_hm.contour(gap_range, cmm_range, Z,
                        levels=[10, 30, 50, 100, 200],
                        colors=['white'], linewidths=0.6, alpha=0.7)
    ax_hm.clabel(cs, fontsize=6, colors='white', fmt='%.0f Pa')
    ax_hm.axhline(CMM_ref, color='#ffc940', lw=1.0, ls=':', alpha=0.7)
    ax_hm.text(97, CMM_ref+0.3, f"CMM={CMM_ref}", fontsize=6,
               color='#ffc940', ha='right')
    ax_hm.set_title(f"ΔP 민감도: Gap × CMM ({labels[0]})", fontsize=9,
                    color=DARK["text"])
    ax_hm.set_xlabel("Gap G [mm]", fontsize=8)
    ax_hm.set_ylabel("CMM [m³/min]", fontsize=8)
    ax_hm.tick_params(labelsize=7)

    fig.suptitle(
        "모듈 C — 공기측 압력강하  │  f-factor × ΔP_core × Gap 손실 × 습면 보정",
        fontsize=13, color=DARK["text"], fontweight='bold', y=0.975)
    fig.text(0.97, 0.975,
             f"{ref.refrigerant}  CMM={CMM_ref} m³/min"
             f"  T_in={inlet.T_in:.0f}°C  RH={inlet.RH_in*100:.0f}%",
             fontsize=8, color=DARK["dim"], ha='right', va='top')
    return fig



# ═══════════════════════════════════════════════════════════════
#  단일 케이스 가시화 — run_single.py용
# ═══════════════════════════════════════════════════════════════

def _styl(ax):
    ax.set_facecolor(DARK["panel"])
    ax.tick_params(colors=DARK["text"], labelsize=9)
    for s in ax.spines.values(): s.set_color(DARK["border"])
    ax.grid(True, color=DARK["grid"], alpha=0.4, lw=0.5)

_C = C  # alias


def make_single_fig_a(gaps, results, gp, ref, tag):
    """
    Module A — 냉방성능 vs Gap (2×3)
    (0,0) Cap Retention + Cap*
    (0,1) Q 분해 (Sensible/Latent stack + Q_net)
    (0,2) 손실 분해 (Q_rad, Q_cond, Q_mix, Q_carry stack)
    (1,0) ΔT_recir
    (1,1) SHR
    (1,2) 제습량
    """
    fig, axes = plt.subplots(2, 3, figsize=(22, 11), facecolor=DARK["bg"])
    fig.subplots_adjust(hspace=0.34, wspace=0.32, left=0.06, right=0.97, top=0.89, bottom=0.07)

    cap   = np.array([r['cap_ratio'] for r in results])
    cap_c = np.array([r.get('cap_corrected', r['cap_ratio']) for r in results])
    Qnet  = np.array([r['Q_net'] for r in results])
    Qsen  = np.array([r['Q_sensible'] for r in results])
    Qlat  = np.array([r['Q_latent'] for r in results])
    Qrad  = np.array([r['Q_rad'] for r in results])
    Qcond = np.array([r['Q_cond'] for r in results])
    Qmix  = np.array([r['Q_mix'] for r in results])
    Qcarry= np.array([r.get('Q_carry_penalty', 0) for r in results])
    dTr   = np.array([r['dT_recir'] for r in results])
    SHR   = np.array([r['SHR'] for r in results])
    deh   = np.array([r['dehumid_rate'] for r in results])

    # ── (0,0) Cap Retention + Cap* ──
    ax = axes[0,0]; _styl(ax)
    ax.fill_between(gaps, cap, alpha=0.12, color=_C[0])
    ax.plot(gaps, cap, color=_C[0], lw=2.5, label='Cap')
    if np.max(np.abs(cap - cap_c)) > 0.1:
        ax.plot(gaps, cap_c, color=_C[4], lw=2, ls='--', label='Cap* (carry)')
    ax.axhline(90, color=_C[3], ls='--', alpha=0.6, lw=1)
    ax.text(gaps[-1]*0.97, 91, '90%', color=_C[3], fontsize=8, ha='right')
    # 90% 도달 Gap 표시
    idx90 = np.where(cap >= 90)[0]
    if len(idx90) > 0:
        g90 = gaps[idx90[0]]
        ax.axvline(g90, color=_C[3], ls=':', alpha=0.4)
        ax.annotate(f'G≥90% = {g90:.0f}mm', (g90, 90), fontsize=8,
                    color=_C[3], xytext=(8, -15), textcoords='offset points')
    ax.set_ylabel('Cap Retention [%]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Capacity Retention', color=DARK["text"], fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='lower right', framealpha=0.3,
              labelcolor=DARK["text"], facecolor=DARK["panel"])

    # ── (0,1) Q 분해 ──
    ax = axes[0,1]; _styl(ax)
    ax.stackplot(gaps, Qsen, Qlat, colors=[_C[0]+'88', _C[4]+'88'],
                 labels=['Sensible', 'Latent'])
    ax.plot(gaps, Qnet, color=_C[3], lw=2.5, label='$Q_{net}$')
    ax.plot(gaps, Qsen+Qlat, color=DARK["text"], lw=1, ls=':', alpha=0.4, label='$Q_{total}$')
    ax.set_ylabel('Q [W]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Heat Transfer Breakdown', color=DARK["text"], fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='best', framealpha=0.3,
              labelcolor=DARK["text"], facecolor=DARK["panel"])

    # ── (0,2) 손실 분해 (stacked area) ──
    ax = axes[0,2]; _styl(ax)
    losses = [Qrad, Qcond, Qmix]
    labels = ['$Q_{rad}$', '$Q_{cond}$', '$Q_{mix}$']
    colors = [_C[1]+'88', _C[5]+'88', _C[2]+'88']
    if np.max(Qcarry) > 0.01:
        losses.append(Qcarry)
        labels.append('$Q_{carry}$')
        colors.append(_C[4]+'88')
    ax.stackplot(gaps, *losses, colors=colors, labels=labels)
    Qtot_loss = sum(losses)
    ax.plot(gaps, Qtot_loss, color=DARK["text"], lw=1.5, ls='--', alpha=0.6, label='Total loss')
    ax.set_ylabel('$Q_{loss}$ [W]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Loss Decomposition', color=DARK["text"], fontsize=11, fontweight='bold')
    ax.legend(fontsize=7, loc='upper right', framealpha=0.3,
              labelcolor=DARK["text"], facecolor=DARK["panel"])

    # ── (1,0) ΔT_recir ──
    ax = axes[1,0]; _styl(ax)
    ax.fill_between(gaps, dTr, alpha=0.12, color=_C[1])
    ax.plot(gaps, dTr, color=_C[1], lw=2.5)
    if np.max(dTr) > 0.01:
        ax.annotate(f'max={np.max(dTr):.1f}°C\n@ G={gaps[np.argmax(dTr)]:.0f}mm',
                    (gaps[np.argmax(dTr)], np.max(dTr)), fontsize=8, color=_C[1],
                    xytext=(15, -5), textcoords='offset points')
    ax.set_ylabel(r'$\Delta T_{recir}$ [°C]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Recirculation Temperature Rise', color=DARK["text"], fontsize=11, fontweight='bold')
    # gap_mode 표시
    ax.text(0.97, 0.95, f'mode: {gp.gap_mode}', transform=ax.transAxes,
            fontsize=8, color=DARK["dim"], ha='right', va='top')

    # ── (1,1) SHR ──
    ax = axes[1,1]; _styl(ax)
    ax.plot(gaps, SHR, color=_C[2], lw=2.5)
    ax.set_ylim(0, 1.05)
    ax.axhline(0.75, color=DARK["dim"], ls=':', alpha=0.3)
    ax.text(gaps[0]*1.05, 0.76, 'SHR=0.75', fontsize=7, color=DARK["dim"])
    mean_shr = np.mean(SHR[SHR > 0]) if np.any(SHR > 0) else 0
    ax.text(0.97, 0.95, f'avg SHR = {mean_shr:.2f}', transform=ax.transAxes,
            fontsize=9, color=_C[2], ha='right', va='top', fontweight='bold')
    ax.set_ylabel('SHR [-]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Sensible Heat Ratio', color=DARK["text"], fontsize=11, fontweight='bold')

    # ── (1,2) 제습량 ──
    ax = axes[1,2]; _styl(ax)
    ax.fill_between(gaps, deh, alpha=0.12, color=_C[4])
    ax.plot(gaps, deh, color=_C[4], lw=2.5)
    mean_deh = np.mean(deh[deh > 0]) if np.any(deh > 0) else 0
    ax.text(0.97, 0.95, f'avg = {mean_deh:.0f} g/h', transform=ax.transAxes,
            fontsize=9, color=_C[4], ha='right', va='top', fontweight='bold')
    ax.set_ylabel('Dehumidification [g/h]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Dehumidification Rate', color=DARK["text"], fontsize=11, fontweight='bold')

    fig.suptitle(f"Module A — Cooling Performance vs Gap\n{tag}",
                 fontsize=14, color=DARK["text"], fontweight='bold', y=0.97)
    fig.savefig("module_a.png", dpi=200, bbox_inches='tight', facecolor=DARK["bg"])
    print("    → module_a.png")
    return fig


def make_single_fig_b(gaps, result_b, viz_b, tag):
    """
    Module B — 비말동반 vs Gap (2×3)
    (0,0) Cap vs Flux (dual axis + sweet-spot)
    (0,1) Carryover Flux vs Gap (+ danger/moderate 기준선)
    (0,2) P_reach vs Gap
    (1,0) Weber 비율 vs V_face (현재 V_face 마커)
    (1,1) Risk Zone Map + Sweet-spot
    (1,2) Q_carry penalty vs Gap

    viz_b: dict with P_reach_arr, We_ratio, V_face, V_onset,
           q_cond, carry_penalty, We_arr_v (V sweep)
    """
    fig, axes = plt.subplots(2, 3, figsize=(22, 11), facecolor=DARK["bg"])
    fig.subplots_adjust(hspace=0.34, wspace=0.35, left=0.06, right=0.97, top=0.89, bottom=0.07)

    cap_arr  = result_b['cap_arr']
    flux_arr = result_b['flux_arr']
    sg       = result_b.get('sweet_gap')
    V_on     = result_b['V_on']
    # 4단계: FLUX_CAUTION / FLUX_DANGER / FLUX_SEVERE (from common)

    P_reach  = viz_b['P_reach_arr']
    We_ratio = viz_b['We_ratio']
    V_face   = viz_b['V_face']
    V_onset  = viz_b['V_onset']
    q_cond   = viz_b['q_cond']
    carry_p  = viz_b['carry_penalty']

    # 위험도 분류
    risk = ['SAFE' if f < FLUX_CAUTION else ('CAUTION' if f < FLUX_DANGER else ('DANGER' if f < FLUX_SEVERE else 'SEVERE'))
            for f in flux_arr]

    # ── (0,0) Cap vs Flux (dual axis) ──
    ax = axes[0,0]; _styl(ax)
    ax.plot(gaps, cap_arr, color=_C[0], lw=2.5, label='Cap [%]')
    ax.axhline(90, color=_C[3], ls='--', alpha=0.4)
    ax.set_ylabel('Cap Retention [%]', color=_C[0], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax2 = ax.twinx()
    ax2.plot(gaps, flux_arr, color=_C[1], lw=2, ls='--', label='Flux')
    ax2.axhline(FLUX_DANGER, color=_C[1], ls=':', alpha=0.3)
    ax2.set_ylabel('Carryover Flux', color=_C[1], fontsize=10)
    ax2.tick_params(colors=_C[1], labelsize=9)
    ax2.spines['right'].set_color(_C[1])
    if sg:
        ax.axvline(sg, color=DARK["text"], ls='--', lw=2, alpha=0.7)
        ax.text(sg+2, 50, f'G*={sg:.0f}mm', color=DARK["text"], fontsize=10, fontweight='bold')
    ax.set_title('Cap vs Flux (Sweet-spot)', color=DARK["text"], fontsize=11, fontweight='bold')
    ax.legend(fontsize=8, loc='center left', framealpha=0.3,
              labelcolor=DARK["text"], facecolor=DARK["panel"])
    ax2.legend(fontsize=8, loc='center right', framealpha=0.3,
               labelcolor=DARK["text"], facecolor=DARK["panel"])

    # ── (0,1) Flux vs Gap ──
    ax = axes[0,1]; _styl(ax)
    ax.fill_between(gaps, flux_arr, alpha=0.12, color=_C[5])
    ax.plot(gaps, flux_arr, color=_C[5], lw=2.5)
    ax.axhline(FLUX_CAUTION, color=_C[3], ls=':', lw=1, alpha=0.5)
    ax.axhline(FLUX_DANGER, color=_C[2], ls='--', lw=1.5, alpha=0.7)
    ax.axhline(FLUX_SEVERE, color=_C[1], ls='--', lw=1.5, alpha=0.5)
    ax.text(gaps[-1]*0.97, FLUX_DANGER*1.08, 'DANGER', color=_C[2], fontsize=8, ha='right')
    ax.text(gaps[-1]*0.97, FLUX_CAUTION*1.15, 'CAUTION', color=_C[3], fontsize=7, ha='right')
    ax.text(gaps[-1]*0.97, FLUX_SEVERE*1.08, 'SEVERE', color=_C[1], fontsize=7, ha='right')
    # max flux annotation
    if np.max(flux_arr) > 0.01:
        ax.annotate(f'max={np.max(flux_arr):.2f}', (gaps[np.argmax(flux_arr)], np.max(flux_arr)),
                    fontsize=8, color=_C[5], xytext=(5, 5), textcoords='offset points')
    ax.set_ylabel('Carryover Flux [-]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Carryover Flux vs Gap', color=DARK["text"], fontsize=11, fontweight='bold')

    # ── (0,2) P_reach vs Gap ──
    ax = axes[0,2]; _styl(ax)
    ax.plot(gaps, P_reach, color=_C[3], lw=2.5)
    ax.fill_between(gaps, P_reach, alpha=0.12, color=_C[3])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel('$P_{reach}$ (도달 확률) [-]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Droplet Reach Probability', color=DARK["text"], fontsize=11, fontweight='bold')
    ax.text(0.97, 0.95, f'Gap=5mm: {P_reach[0]:.2f}\nGap={gaps[-1]:.0f}mm: {P_reach[-1]:.2f}',
            transform=ax.transAxes, fontsize=9, color=_C[3], ha='right', va='top', fontweight='bold')

    # ── (1,0) Weber ratio vs V_face ──
    ax = axes[1,0]; _styl(ax)
    V_sweep = viz_b.get('V_sweep', np.linspace(0.5, 6, 50))
    We_sweep = viz_b.get('We_sweep', np.ones_like(V_sweep))
    ax.plot(V_sweep, We_sweep, color=_C[4], lw=2.5)
    ax.axhline(1.0, color=_C[2], ls='--', lw=1.5, alpha=0.7)
    ax.text(V_sweep[-1]*0.97, 1.08, '$We/We_{crit} = 1$', fontsize=8, color=_C[2], ha='right')
    ax.axvline(V_onset, color=_C[3], ls=':', lw=1.5, alpha=0.5)
    ax.text(V_onset+0.05, np.max(We_sweep)*0.9,
            f'$V_{{onset}}$={V_onset:.2f} m/s', fontsize=8, color=_C[3])
    # 현재 V_face 마커
    ax.scatter([V_face], [We_ratio], s=150, c=_C[1], marker='*', zorder=5, edgecolors=DARK["text"], lw=1.5)
    ax.annotate(f'V={V_face:.1f}\nWe/Wc={We_ratio:.2f}',
                (V_face, We_ratio), fontsize=8, color=_C[1],
                xytext=(10, 10), textcoords='offset points')
    ax.fill_between(V_sweep, 0, 1, alpha=0.05, color=_C[3])
    ax.text(V_onset*0.5, 0.5, 'SAFE', fontsize=10, color=_C[3], ha='center', alpha=0.5)
    ax.set_ylabel('$We / We_{crit}$ [-]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('$V_{face}$ [m/s]', color=DARK["text"])
    ax.set_title('Weber Number Ratio', color=DARK["text"], fontsize=11, fontweight='bold')

    # ── (1,1) Risk Zone Map ──
    ax = axes[1,1]; _styl(ax)
    cmap = {'SAFE': _C[3], 'CAUTION': _C[2], 'DANGER': _C[1], 'SEVERE': '#7b0000'}
    bw = (gaps[-1]-gaps[0])/len(gaps)*0.9
    for g, r in zip(gaps, risk):
        ax.bar(g, 1, width=bw, color=cmap[r], alpha=0.6)
    if sg:
        ax.axvline(sg, color=DARK["text"], ls='--', lw=2.5, alpha=0.9)
        ax.text(sg+2, 0.5, f'G*={sg:.0f}mm', color=DARK["text"], fontsize=11, fontweight='bold', va='center')
    ax.set_yticks([])
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Risk Zone Map', color=DARK["text"], fontsize=11, fontweight='bold')
    from matplotlib.patches import Patch
    ax.legend(handles=[Patch(color=_C[3], label='SAFE'), Patch(color=_C[2], label='CAUTION'),
                       Patch(color=_C[1], label='DANGER'), Patch(color='#7b0000', label='SEVERE')],
              fontsize=8, loc='upper right', framealpha=0.4,
              labelcolor=DARK["text"], facecolor=DARK["panel"])

    # ── (1,2) Q_carry penalty ──
    ax = axes[1,2]; _styl(ax)
    ax.fill_between(gaps, carry_p, alpha=0.12, color=_C[4])
    ax.plot(gaps, carry_p, color=_C[4], lw=2.5)
    if np.max(carry_p) > 0.01:
        ax.annotate(f'max={np.max(carry_p):.1f}W', (gaps[np.argmax(carry_p)], np.max(carry_p)),
                    fontsize=8, color=_C[4], xytext=(5, 5), textcoords='offset points')
    else:
        ax.text(0.5, 0.5, f'V_face < V_onset\n→ No carryover\n→ Q_carry = 0',
                transform=ax.transAxes, fontsize=10, color=DARK["dim"],
                ha='center', va='center')
    # q_cond 표시
    ax.text(0.03, 0.95, f'$q_{{cond}}$ = {q_cond:.2f} g/s\n$V_{{onset}}$ = {V_onset:.2f} m/s',
            transform=ax.transAxes, fontsize=8, color=DARK["dim"], va='top')
    ax.set_ylabel('$Q_{carry}$ [W]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('Condenser Carry Penalty', color=DARK["text"], fontsize=11, fontweight='bold')

    fig.suptitle(f"Module B — Carryover Risk vs Gap\n{tag}",
                 fontsize=14, color=DARK["text"], fontweight='bold', y=0.97)
    fig.savefig("module_b.png", dpi=200, bbox_inches='tight', facecolor=DARK["bg"])
    print("    → module_b.png")
    return fig


def make_single_fig_c(gaps, rc, inlet, tag):
    """
    Module C — 압력강하 vs Gap (2×3)
    (0,0) ΔP_total vs Gap
    (0,1) ΔP stacked area (evap/gap/cond)
    (0,2) ΔP_gap / ΔP_total 비율 [%]
    (1,0) ΔP_gap vs Gap (상세)
    (1,1) f_dry vs f_wet 비교 (bar)
    (1,2) 핵심 수치 요약 패널
    """
    fig, axes = plt.subplots(2, 3, figsize=(22, 11), facecolor=DARK["bg"])
    fig.subplots_adjust(hspace=0.34, wspace=0.32, left=0.06, right=0.97, top=0.89, bottom=0.07)

    dp_t = rc['dp_total']
    dp_e = rc['dp_evap']
    dp_g = rc['dp_gap']
    dp_c = rc['dp_cond']
    f_dry = rc['f_evap_dry']
    f_wet = rc['f_evap_wet']
    f_ratio = rc.get('f_wet_ratio', f_wet/(f_dry+1e-9))

    # ── (0,0) ΔP_total ──
    ax = axes[0,0]; _styl(ax)
    ax.fill_between(gaps, dp_t, alpha=0.12, color=_C[0])
    ax.plot(gaps, dp_t, color=_C[0], lw=2.5)
    ax.annotate(f'G=5: {dp_t[0]:.1f} Pa', (gaps[0], dp_t[0]),
                fontsize=8, color=_C[0], xytext=(10, 5), textcoords='offset points')
    ax.annotate(f'G={gaps[-1]:.0f}: {dp_t[-1]:.1f} Pa', (gaps[-1], dp_t[-1]),
                fontsize=8, color=_C[0], xytext=(-60, 5), textcoords='offset points')
    ax.set_ylabel(r'$\Delta P_{system}$ [Pa]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title('System Pressure Drop', color=DARK["text"], fontsize=12, fontweight='bold')

    # ── (0,1) ΔP stacked ──
    ax = axes[0,1]; _styl(ax)
    ax.stackplot(gaps, dp_e, dp_g, dp_c,
                 colors=[_C[0]+'88', _C[2]+'88', _C[1]+'88'],
                 labels=[r'$\Delta P_{evap}$', r'$\Delta P_{gap}$', r'$\Delta P_{cond}$'])
    ax.plot(gaps, dp_t, color=DARK["text"], lw=1.5, ls='--', alpha=0.6, label='Total')
    ax.set_ylabel(r'$\Delta P$ [Pa]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title(r'$\Delta P$ Component Stack', color=DARK["text"], fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='upper right', framealpha=0.3,
              labelcolor=DARK["text"], facecolor=DARK["panel"])

    # ── (0,2) Gap 비율 ──
    ax = axes[0,2]; _styl(ax)
    frac_g = dp_g / (dp_t + 1e-9) * 100
    frac_e = dp_e / (dp_t + 1e-9) * 100
    frac_c = dp_c / (dp_t + 1e-9) * 100
    ax.stackplot(gaps, frac_e, frac_g, frac_c,
                 colors=[_C[0]+'55', _C[2]+'55', _C[1]+'55'],
                 labels=['Evap', 'Gap', 'Cond'])
    ax.set_ylim(0, 100)
    ax.set_ylabel(r'$\Delta P$ fraction [%]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title(r'$\Delta P$ Component Fraction', color=DARK["text"], fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='center right', framealpha=0.3,
              labelcolor=DARK["text"], facecolor=DARK["panel"])
    # 소간격 Gap 비율 annotation
    ax.annotate(f'Gap: {frac_g[0]:.0f}%', (gaps[0], frac_e[0]+frac_g[0]/2),
                fontsize=8, color=_C[2], xytext=(10, 0), textcoords='offset points')

    # ── (1,0) ΔP_gap 상세 ──
    ax = axes[1,0]; _styl(ax)
    ax.fill_between(gaps, dp_g, alpha=0.12, color=_C[2])
    ax.plot(gaps, dp_g, color=_C[2], lw=2.5, label=r'$\Delta P_{gap}$')
    ax2 = ax.twinx()
    ax2.plot(gaps, frac_g, color=_C[5], lw=1.5, ls='--', label='Gap/Total %')
    ax2.set_ylabel('Gap / Total [%]', color=_C[5], fontsize=10)
    ax2.tick_params(colors=_C[5], labelsize=9)
    ax2.spines['right'].set_color(_C[5])
    ax.set_ylabel(r'$\Delta P_{gap}$ [Pa]', color=DARK["text"], fontsize=10)
    ax.set_xlabel('Gap [mm]', color=DARK["text"])
    ax.set_title(r'$\Delta P_{gap}$ Detail', color=DARK["text"], fontsize=12, fontweight='bold')
    ax.legend(fontsize=8, loc='upper left', framealpha=0.3,
              labelcolor=DARK["text"], facecolor=DARK["panel"])
    ax2.legend(fontsize=8, loc='upper right', framealpha=0.3,
               labelcolor=DARK["text"], facecolor=DARK["panel"])

    # ── (1,1) f_dry vs f_wet bar ──
    ax = axes[1,1]; _styl(ax)
    bars_x = [0, 1]
    bars_h = [f_dry, f_wet]
    bars_c = [_C[0], _C[1]]
    bars_l = ['$f_{dry}$', '$f_{wet}$']
    brs = ax.bar(bars_x, bars_h, width=0.5, color=bars_c, alpha=0.7,
                 edgecolor=DARK["text"], lw=0.5)
    for b, v, l in zip(brs, bars_h, bars_l):
        ax.text(b.get_x()+b.get_width()/2, v+0.001, f'{v:.4f}',
                ha='center', fontsize=10, color=DARK["text"], fontweight='bold')
        ax.text(b.get_x()+b.get_width()/2, -0.004, l,
                ha='center', fontsize=10, color=DARK["text"])
    # f_wet/f_dry ratio annotation
    ax.text(0.5, max(bars_h)*0.6,
            f'$f_{{wet}}/f_{{dry}}$ = {f_ratio:.3f}',
            fontsize=12, color=_C[2], fontweight='bold', ha='center',
            bbox=dict(boxstyle='round,pad=0.3', fc=DARK["panel"], ec=DARK["border"]))
    ax.set_xticks([])
    ax.set_ylabel('Fanning f [-]', color=DARK["text"], fontsize=10)
    ax.set_title('Friction Factor: Dry vs Wet', color=DARK["text"], fontsize=12, fontweight='bold')

    # ── (1,2) 수치 요약 패널 ──
    ax = axes[1,2]; _styl(ax)
    ax.axis('off')
    idx20 = np.argmin(np.abs(gaps - 20)) if len(gaps) > 0 else 0
    idx50 = np.argmin(np.abs(gaps - 50)) if len(gaps) > 0 else 0
    V = rc.get('V_face', 0)
    summary_lines = [
        (f'$V_{{face}}$ = {V:.2f} m/s', _C[0]),
        (f'$f_{{dry}}$ = {f_dry:.4f}   $f_{{wet}}$ = {f_wet:.4f}', DARK["text"]),
        (f'$f_{{wet}}/f_{{dry}}$ = {f_ratio:.3f}', _C[2]),
        ('', ''),
        (f'G=5mm:  $\\Delta P$ = {dp_t[0]:.1f} Pa  (gap {frac_g[0]:.0f}%)', DARK["text"]),
        (f'G=20mm: $\\Delta P$ = {dp_t[idx20]:.1f} Pa  (gap {frac_g[idx20]:.0f}%)', DARK["text"]),
        (f'G=50mm: $\\Delta P$ = {dp_t[idx50]:.1f} Pa  (gap {frac_g[idx50]:.0f}%)', DARK["text"]),
        ('', ''),
        (f'$\\Delta P_{{evap}}$ = {dp_e[0]:.1f} Pa (fixed)', _C[0]),
        (f'$\\Delta P_{{cond}}$ = {dp_c[0]:.1f} Pa (fixed)', _C[1]),
    ]
    for i, (txt, col) in enumerate(summary_lines):
        if txt:
            ax.text(0.05, 0.92 - i*0.09, txt, transform=ax.transAxes,
                    fontsize=10, color=col, va='top')
    ax.set_title('Key Numbers', color=DARK["text"], fontsize=12, fontweight='bold')

    fig.suptitle(f"Module C — Pressure Drop vs Gap\n{tag}",
                 fontsize=14, color=DARK["text"], fontweight='bold', y=0.97)
    fig.savefig("module_c.png", dpi=200, bbox_inches='tight', facecolor=DARK["bg"])
    print("    → module_c.png")
    return fig


def make_single_schematic(evap_spec, cond_spec, evap_geo, cond_geo,
                          ua_ev, ua_cd, gp, ref, result_a_ref, viz_b, tag):
    """
    Gap 물리 현상 개략도 — _combo_schematic 스타일 확장
    비말동반 대표 궤적 1개 포함 (해석결과 반영)

    Parameters
    ----------
    result_a_ref : simulate_gap 결과 (대표 Gap, e.g. G=20mm)
    viz_b : dict with V_face, V_onset, q_cond, P_reach_arr, carry_penalty, ...
    """
    from module_b import CarryoverSpec, MCHXCarryoverSpec, monte_carlo, eta_co, v_onset, we_ch, we_crit

    fig, ax = plt.subplots(1, 1, figsize=(14, 8), facecolor=DARK["bg"])
    ax.set_xlim(0, 14)
    ax.set_ylim(-1.2, 8.8)
    ax.set_facecolor(DARK["bg"])
    ax.axis('off')

    r = result_a_ref  # 대표 Gap 결과
    gap_mm = r['gap_mm']
    V_face = viz_b['V_face']
    V_on   = viz_b['V_onset']
    q_cond = viz_b['q_cond']
    carry_at_g = viz_b['carry_penalty'][np.argmin(np.abs(viz_b['gaps'] - gap_mm))] if 'gaps' in viz_b else 0

    # ══════════════════════════════════════════════════════════
    #  HX 블록
    # ══════════════════════════════════════════════════════════

    def _hx_block(x0, spec, ua_d, side):
        clr = C[0] if side == "evap" else C[1]
        tag_hx = "EVAP" if side == "evap" else "COND"
        typ = spec.hx_type
        rect = mpatches.FancyBboxPatch((x0, 1.0), 2.8, 5.0,
                   boxstyle="round,pad=0.08", fc=clr+"18", ec=clr, lw=2.0)
        ax.add_patch(rect)
        # 핀
        for y in np.linspace(1.3, 5.7, 9):
            ax.plot([x0+0.2, x0+2.6], [y, y], color=clr, lw=0.7, alpha=0.4)
        # 튜브
        for cy in np.linspace(2.0, 5.0, 3):
            for cx_off in np.linspace(0.6, 2.2, 3):
                ax.add_patch(plt.Circle((x0+cx_off, cy), 0.32, fc=clr+"33", ec=clr, lw=1.0))
        # 라벨
        ax.text(x0+1.4, 6.3, f"{tag_hx}\n{typ}", ha='center', fontsize=11,
                color=clr, fontweight='bold')
        # UA 정보
        ax.text(x0+1.4, 0.3,
                f"UA={ua_d['UA']:.0f}W/K\n"
                f"$h_o$={ua_d['h_o']:.0f}  $h_i$={ua_d['h_i']:.0f}",
                ha='center', fontsize=8, color=DARK["dim"])
        # 핀 타입
        ft = getattr(spec, 'fin_type', 'MCHX' if spec.hx_type == 'MCHX' else 'plain')
        layout = getattr(spec, 'tube_layout', 'mchx')[0].upper()
        ax.text(x0+1.4, 0.0, f"{ft} / {layout}",
                ha='center', fontsize=7, color=clr, alpha=0.6)

    _hx_block(0.3, evap_spec, ua_ev, "evap")
    _hx_block(10.9, cond_spec, ua_cd, "cond")

    evap_right = 3.1
    cond_left = 10.9

    # ══════════════════════════════════════════════════════════
    #  Gap 영역
    # ══════════════════════════════════════════════════════════

    # Gap 양방향 화살표
    ax.annotate('', xy=(cond_left-0.05, 3.5), xytext=(evap_right+0.05, 3.5),
                arrowprops=dict(arrowstyle='<->', color=C[2], lw=2.0))
    ax.text(7.0, 3.75, f"Gap G = {gap_mm:.0f} mm", ha='center', fontsize=11,
            color=C[2], fontweight='bold')

    # ── ① 재순환 (Recirculation) ──
    ax.annotate('', xy=(evap_right+0.05, 5.5), xytext=(cond_left-0.05, 5.5),
                arrowprops=dict(arrowstyle='->', color=C[5], lw=2.0,
                                connectionstyle='arc3,rad=0.45'))
    ax.text(7.0, 7.3, '① 재순환', ha='center', fontsize=10, color=C[5], fontweight='bold')
    dTr = r['dT_recir']
    ax.text(7.0, 6.85, f'$\\Delta T_{{recir}}$ = {dTr:.2f}°C  ({gp.gap_mode})',
            ha='center', fontsize=8, color=C[5], alpha=0.85)

    # ── ② 공기 주 흐름 ──
    ax.annotate('', xy=(cond_left-0.1, 2.5), xytext=(evap_right+0.1, 2.5),
                arrowprops=dict(arrowstyle='->', color='#5599cc', lw=2.0))
    ax.text(7.0, 2.2, f'공기 흐름 →  V={V_face:.2f} m/s  CMM={gp.CMM}',
            ha='center', fontsize=8, color='#5599cc')

    # ── ③ 복사 (wavy line) ──
    xr = np.linspace(evap_right+0.2, cond_left-0.2, 80)
    yr = 4.6 + 0.08 * np.sin(xr * 10)
    ax.plot(xr, yr, color=C[1], lw=1.2, ls='--', alpha=0.7)
    Q_rad = r['Q_rad']
    ax.text(7.0, 4.9, f'③ 복사  $Q_{{rad}}$ = {Q_rad:.1f}W',
            ha='center', fontsize=8.5, color=C[1], fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.15', fc=DARK["bg"]+'cc', ec=C[1], alpha=0.8))

    # ── ④ 전도 단락 (bottom bar) ──
    ax.plot([evap_right+0.1, cond_left-0.1], [1.2, 1.2], color=C[2], lw=3.0,
            solid_capstyle='round', alpha=0.7)
    ax.annotate('', xy=(cond_left-0.1, 1.2), xytext=(7.0, 1.2),
                arrowprops=dict(arrowstyle='->', color=C[2], lw=1.5))
    Q_cond = r['Q_cond']
    ax.text(7.0, 0.8, f'④ 전도단락  $Q_{{cond}}$ = {Q_cond:.1f}W',
            ha='center', fontsize=8.5, color=C[2], fontweight='bold')

    # ── ⑤ 비말동반 궤적 (대표 1개) ──
    # 실제 MC 궤적 계산
    if evap_spec.hx_type == 'MCHX':
        co_spec = MCHXCarryoverSpec(evap_spec, evap_geo)
    else:
        co_spec = CarryoverSpec(evap_spec, evap_geo)
    from module_b import track_batch, sample_bimodal
    V_onset_val = v_onset(co_spec) if evap_spec.hx_type == 'FT' else 2.5
    eta_val = eta_co(V_face, co_spec, 5e-4) if V_face > V_onset_val else 0

    if eta_val > 0:
        # 비말동반 활성 → 도달 궤적 표시
        rng = np.random.RandomState(42)
        d_arr = sample_bimodal(50, rng, 0.75)
        h_drain = getattr(evap_spec, 'h_drain', 0.025)
        y0_arr = rng.uniform(-evap_spec.H/2 + h_drain, evap_spec.H/2, 50)
        gap_m = gap_mm / 1000.0
        reached_mask, y_hit = track_batch(d_arr, y0_arr, gap_m, V_face, co_spec, 0.0)

        if np.any(reached_mask):
            idx_r = np.where(reached_mask)[0][0]
            n_pts = 25
            x_vis = np.linspace(evap_right+0.15, cond_left-0.15, n_pts)
            y0_vis = 1.5 + (y0_arr[idx_r] / evap_spec.H + 0.5) * 4.0
            np.random.seed(idx_r)
            y_vis = y0_vis + np.cumsum(np.random.normal(0, 0.03, n_pts)) - 0.015*np.arange(n_pts)
            ax.plot(x_vis, y_vis, color=C[0], lw=2.0, alpha=0.8)
            ax.scatter(x_vis[-1], y_vis[-1], s=60, c=C[1], marker='o',
                       edgecolors=DARK["text"], lw=1.0, zorder=5)
            ax.scatter(x_vis[0], y_vis[0], s=30, c=C[0], marker='o',
                       edgecolors=DARK["text"], lw=0.8, zorder=5)
            ax.annotate('도달!', (x_vis[-1], y_vis[-1]),
                        fontsize=7, color=C[1], fontweight='bold',
                        xytext=(5, 5), textcoords='offset points')
    else:
        # η_co = 0 → 액적 이탈 없음 → 응결수는 표면에 머무름
        # 낙하 궤적 (이탈 실패 시뮬레이션)
        n_pts = 15
        x_vis = np.linspace(evap_right+0.15, evap_right + 1.5, n_pts)
        y0_vis = 3.5
        np.random.seed(7)
        y_vis = y0_vis - 0.08*np.arange(n_pts) + np.cumsum(np.random.normal(0, 0.015, n_pts))
        ax.plot(x_vis, y_vis, color=C[3], lw=1.5, alpha=0.5, ls='--')
        ax.scatter(x_vis[-1], y_vis[-1], s=30, c=C[3], marker='v', zorder=5)
        ax.annotate('낙하 (이탈 불가)', (x_vis[-1], y_vis[-1]),
                    fontsize=7, color=C[3], xytext=(5, -10), textcoords='offset points')

    # 응결수 표시 (증발기 출구면)
    for yy in [2.0, 3.0, 4.0, 5.0]:
        ax.scatter(evap_right+0.05, yy, s=20, c='#3399cc', alpha=0.8, marker='o')

    # 비말동반 정보 박스
    P_r = viz_b['P_reach_arr'][np.argmin(np.abs(viz_b.get('gaps', np.array([gap_mm])) - gap_mm))] if 'P_reach_arr' in viz_b else 0
    carry_status = "ACTIVE" if eta_val > 0 else "SAFE (V < V_onset)"

    ax.text(7.0, -0.2,
            f'⑤ 비말동반  $q_{{cond}}$={q_cond:.2f}g/s  '
            f'$V_{{onset}}$={V_on:.2f}m/s  '
            f'$\\eta_{{co}}$={eta_val:.1e}  '
            f'$P_{{reach}}$={P_r:.2f}  '
            f'$Q_{{carry}}$={carry_at_g:.1f}W',
            ha='center', fontsize=7.5, color=C[0],
            bbox=dict(boxstyle='round,pad=0.2', fc=DARK["bg"]+'cc', ec=C[0], alpha=0.8))
    ax.text(7.0, -0.6, f'[{carry_status}]',
            ha='center', fontsize=8, color=C[3] if eta_val == 0 else C[1], fontweight='bold')

    # ══════════════════════════════════════════════════════════
    #  요약 정보 박스 (좌상단 — EVAP 위)
    # ══════════════════════════════════════════════════════════
    info_lines = [
        (f'{ref.refrigerant}  T_evap={ref.T_sat_evap:.0f}°C  T_cond={ref.T_sat_cond:.0f}°C', C[4]),
        (f'Cap = {r["cap_ratio"]:.1f}%   Q_net = {r["Q_net"]:.0f}W   SHR = {r["SHR"]:.2f}', C[0]),
        (f'Q_total = {r["Q_total"]:.0f}W  (sen={r["Q_sensible"]:.0f} + lat={r["Q_latent"]:.0f})   제습={r["dehumid_rate"]:.0f}g/h', DARK["text"]),
        (f'Losses: rad={Q_rad:.1f} + cond={Q_cond:.1f} + mix={r["Q_mix"]:.1f} + carry={carry_at_g:.1f}W', DARK["dim"]),
    ]
    for i, (txt, clr) in enumerate(info_lines):
        ax.text(0.2, 8.5 - i*0.35, txt, ha='left', fontsize=7.5,
                color=clr, fontweight='bold' if i < 2 else 'normal')

    fig.suptitle(f'Gap 물리 현상도  |  {tag}  |  G = {gap_mm:.0f}mm',
                 fontsize=13, color=DARK["text"], fontweight='bold', y=0.98)

    fig.savefig("gap_schematic.png", dpi=200, bbox_inches='tight', facecolor=DARK["bg"])
    print("    → gap_schematic.png")
    return fig
