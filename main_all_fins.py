"""
main_all_fins.py — 4종 핀 × 3종 모듈 = 12개 가시화 일괄 생성

[핀 타입]  Plain / Wavy / Louvered / Slit
[모듈]     A (냉방성능) / B (비말동반) / C (압력강하)
"""

import numpy as np
import warnings
warnings.filterwarnings("ignore")

from common import (
    FinTubeSpec, MCHXSpec, RefrigerantState,
    compute_ft_geometry, compute_mchx_geometry, compute_UA,
    make_mchx_evap, FLUX_DANGER, C,
)
from module_a import GapParams, simulate_gap, sweep
from module_b import analyze_combined
from module_c import InletCondition, sweep_dp
from visualize import make_module_a_figure, make_module_b_figure, make_module_c_figure


# ═══════════════════════════════════════════════════════════════════
#  공통 설정
# ═══════════════════════════════════════════════════════════════════

REFRIGERANT = "R410A"
T_SAT_EVAP  = 5.0
T_SAT_COND  = 45.0
CMM_REF     = 6.75
GAPS        = np.linspace(1, 100, 20)

# 4종 핀 정의
FIN_TYPES = {
    "plain": dict(
        fin_type="plain",
        label="Plain (평판핀)",
    ),
    "wavy": dict(
        fin_type="wavy",
        wavy_height=1.5e-3,
        wavy_angle=17.5,
        label="Wavy (파형핀)",
    ),
    "slit": dict(
        fin_type="slit",
        slit_num=6,
        slit_height=1.0e-3,
        label="Slit (랜스형 슬릿핀)",
    ),
    "louvered": dict(
        fin_type="louvered",
        ft_louver_pitch=1.7e-3,
        ft_louver_angle=28.0,
        label="Louvered (루버핀)",
    ),
}

# 공통 FT 기하 (핀 타입만 바꿈)
FT_BASE_EVAP = dict(
    W=0.30, H=0.25, D=0.045,
    fin_pitch=1.8e-3, fin_thickness=0.10e-3, fin_height=25.4e-3,
    tube_do=9.52e-3, tube_di=8.52e-3, tube_rows=2, tube_cols=8,
    tube_pitch_t=25.4e-3, tube_pitch_l=22.0e-3,
    tube_layout="staggered",
    fin_material="Al", tube_material="Cu", h_drain=0.025,
)
FT_BASE_COND = dict(
    W=0.35, H=0.28, D=0.050,
    fin_pitch=1.5e-3, fin_thickness=0.10e-3, fin_height=25.4e-3,
    tube_do=9.52e-3, tube_di=8.52e-3, tube_rows=2, tube_cols=9,
    tube_pitch_t=25.4e-3, tube_pitch_l=22.0e-3,
    tube_layout="staggered",
    fin_material="Al", tube_material="Cu", h_drain=0.025,
)


def make_cases(fin_key, fin_cfg):
    """핀 타입별 비말동반 분석 케이스 생성"""
    ft = fin_cfg['fin_type']
    lbl = fin_cfg['label']
    return [
        dict(name=f"Case 1 — {lbl} 기준",
             T_in=27.0, RH_in=0.70, CMM=13.5, T_wall=5.0,
             gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
             gap_mode="semi", seal_fraction=0.7,
             hx_type="FT", tube_layout="staggered"),
        dict(name=f"Case 2 — {lbl} 고풍속/고온고습",
             T_in=30.0, RH_in=0.85, CMM=18.0, T_wall=5.0,
             gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
             gap_mode="semi", seal_fraction=0.7,
             hx_type="FT", tube_layout="staggered"),
        dict(name=f"Case 3 — {lbl} 친수코팅",
             T_in=27.0, RH_in=0.70, CMM=15.75, T_wall=8.0,
             gap_mm=50, theta_deg=15.0, eta_coeff=1e-4, w_bridge=0.60,
             gap_mode="semi", seal_fraction=0.7,
             hx_type="FT", tube_layout="staggered"),
        dict(name=f"Case 4 — {lbl} Inline 기준",
             T_in=27.0, RH_in=0.70, CMM=13.5, T_wall=5.0,
             gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
             gap_mode="semi", seal_fraction=0.7,
             hx_type="FT", tube_layout="inline"),
    ]


# ═══════════════════════════════════════════════════════════════════
#  핀별 시뮬레이션 실행
# ═══════════════════════════════════════════════════════════════════

def run_single_fin(fin_key, fin_cfg, ref, gaps):
    """단일 핀 타입에 대해 모듈 A/B/C 전체 실행"""
    label = fin_cfg['label']
    print(f"\n{'='*70}")
    print(f"  핀 타입: {label} ({fin_key})")
    print(f"{'='*70}")

    # ── HX 스펙 생성 ────────────────────────────────────────
    fin_params = {k: v for k, v in fin_cfg.items() if k != 'label'}

    evap_s = FinTubeSpec(**{**FT_BASE_EVAP, **fin_params, 'tube_layout': 'staggered'})
    evap_i = FinTubeSpec(**{**FT_BASE_EVAP, **fin_params, 'tube_layout': 'inline'})
    cond_s = FinTubeSpec(**{**FT_BASE_COND, **fin_params, 'tube_layout': 'staggered'})
    cond_i = FinTubeSpec(**{**FT_BASE_COND, **fin_params, 'tube_layout': 'inline'})
    cond_mchx = MCHXSpec(W=0.35, H=0.25, D=0.018,
        ch_width=1.0e-3, ch_height=2.0e-3, ch_wall=0.5e-3, n_ports=14,
        fin_pitch=1.2e-3, fin_thickness=0.07e-3,
        louver_pitch=1.4e-3, louver_angle=25.0, n_slabs=1, slab_pitch=10.0e-3)

    # ── 기하 + UA ────────────────────────────────────────────
    geo_es  = compute_ft_geometry(evap_s)
    geo_ei  = compute_ft_geometry(evap_i)
    geo_cs  = compute_ft_geometry(cond_s)
    geo_ci  = compute_ft_geometry(cond_i)
    geo_cm  = compute_mchx_geometry(cond_mchx)

    Vref_ft = CMM_REF / (60.0 * evap_s.W * evap_s.H)
    ua_es  = compute_UA(evap_s, geo_es, ref, Vref_ft, 'evap')
    ua_ei  = compute_UA(evap_i, geo_ei, ref, Vref_ft, 'evap')
    ua_cs  = compute_UA(cond_s, geo_cs, ref, Vref_ft, 'cond')
    ua_ci  = compute_UA(cond_i, geo_ci, ref, Vref_ft, 'cond')
    ua_cm  = compute_UA(cond_mchx, geo_cm, ref, CMM_REF/(60*cond_mchx.W*cond_mchx.H), 'cond')

    print(f"  UA(evap S)={ua_es['UA']:.1f}  h_o={ua_es['h_o']:.1f}")
    print(f"  UA(evap I)={ua_ei['UA']:.1f}  h_o={ua_ei['h_o']:.1f}")

    # ══════════════════════════════════════════════════════════
    #  모듈 A: 6조합
    # ══════════════════════════════════════════════════════════
    print("  [A] 냉방성능...")
    combos = [
        dict(label=f"FT(S)/FT(S) open", color=C[0], ls="-",
             evap_spec=evap_s, evap_geo=geo_es, ua_ev=ua_es,
             cond_spec=cond_s, cond_geo=geo_cs, ua_cd=ua_cs,
             gp=GapParams(CMM=CMM_REF, RH_in=0.50, mode="forced", gap_mode="open")),
        dict(label=f"FT(S)/FT(S) semi", color=C[0], ls="--",
             evap_spec=evap_s, evap_geo=geo_es, ua_ev=ua_es,
             cond_spec=cond_s, cond_geo=geo_cs, ua_cd=ua_cs,
             gp=GapParams(CMM=CMM_REF, RH_in=0.50, mode="forced", gap_mode="semi", seal_fraction=0.7)),
        dict(label=f"FT(S)/FT(S) sealed", color=C[0], ls=":",
             evap_spec=evap_s, evap_geo=geo_es, ua_ev=ua_es,
             cond_spec=cond_s, cond_geo=geo_cs, ua_cd=ua_cs,
             gp=GapParams(CMM=CMM_REF, RH_in=0.50, mode="forced", gap_mode="sealed")),
        dict(label=f"FT(I)/FT(I) semi", color=C[3], ls="-.",
             evap_spec=evap_i, evap_geo=geo_ei, ua_ev=ua_ei,
             cond_spec=cond_i, cond_geo=geo_ci, ua_cd=ua_ci,
             gp=GapParams(CMM=CMM_REF, RH_in=0.50, mode="forced", gap_mode="semi", seal_fraction=0.7)),
        dict(label=f"FT(S)/MCHX semi", color=C[1], ls="--",
             evap_spec=evap_s, evap_geo=geo_es, ua_ev=ua_es,
             cond_spec=cond_mchx, cond_geo=geo_cm, ua_cd=ua_cm,
             gp=GapParams(CMM=CMM_REF, RH_in=0.50, mode="forced", gap_mode="semi", seal_fraction=0.7)),
        dict(label=f"MCHX/MCHX semi", color=C[4], ls="-.",
             evap_spec=MCHXSpec(W=0.30, H=0.20, D=0.020,
                 ch_width=0.8e-3, ch_height=1.5e-3, ch_wall=0.4e-3, n_ports=12,
                 fin_pitch=1.0e-3, fin_thickness=0.06e-3,
                 louver_pitch=1.2e-3, louver_angle=27.0, n_slabs=1, slab_pitch=8.0e-3),
             evap_geo=compute_mchx_geometry(MCHXSpec(W=0.30, H=0.20, D=0.020,
                 ch_width=0.8e-3, ch_height=1.5e-3, ch_wall=0.4e-3, n_ports=12,
                 fin_pitch=1.0e-3, fin_thickness=0.06e-3,
                 louver_pitch=1.2e-3, louver_angle=27.0, n_slabs=1, slab_pitch=8.0e-3)),
             ua_ev=compute_UA(MCHXSpec(W=0.30, H=0.20, D=0.020,
                 ch_width=0.8e-3, ch_height=1.5e-3, ch_wall=0.4e-3, n_ports=12,
                 fin_pitch=1.0e-3, fin_thickness=0.06e-3,
                 louver_pitch=1.2e-3, louver_angle=27.0, n_slabs=1, slab_pitch=8.0e-3),
                 compute_mchx_geometry(MCHXSpec(W=0.30, H=0.20, D=0.020,
                 ch_width=0.8e-3, ch_height=1.5e-3, ch_wall=0.4e-3, n_ports=12,
                 fin_pitch=1.0e-3, fin_thickness=0.06e-3,
                 louver_pitch=1.2e-3, louver_angle=27.0, n_slabs=1, slab_pitch=8.0e-3)),
                 ref, CMM_REF/(60*0.30*0.20), 'evap'),
             cond_spec=cond_mchx, cond_geo=geo_cm, ua_cd=ua_cm,
             gp=GapParams(CMM=CMM_REF, RH_in=0.50, mode="forced", gap_mode="semi", seal_fraction=0.7)),
    ]
    for cb in combos:
        cb['results'] = sweep(gaps, cb['evap_spec'], cb['cond_spec'],
                              cb['evap_geo'], cb['cond_geo'],
                              cb['ua_ev'], cb['ua_cd'], ref, cb['gp'])

    # ══════════════════════════════════════════════════════════
    #  모듈 B: 4케이스
    # ══════════════════════════════════════════════════════════
    print("  [B] 비말동반...")
    cases = make_cases(fin_key, fin_cfg)
    case_colors = [C[0], C[1], C[2], C[3]]
    combined = []
    for case, cc in zip(cases, case_colors):
        if case.get('tube_layout','staggered') == 'inline':
            res = analyze_combined(case, gaps, evap_i, geo_ei, cond_i, geo_ci,
                                   ua_ei, ua_ci, ref)
        else:
            res = analyze_combined(case, gaps, evap_s, geo_es, cond_s, geo_cs,
                                   ua_es, ua_cs, ref)
        combined.append(res)

    # ── A+B 통합: 비말동반 → 응축기 잠열 페널티 ──
    from module_b import compute_carry_penalty, apply_carry_penalty
    for res in combined:
        case = res['case']
        if case.get('tube_layout','staggered') == 'inline':
            ev_sp, ev_geo = evap_i, geo_ei
        else:
            ev_sp, ev_geo = evap_s, geo_es
        penalty = compute_carry_penalty(res, ev_sp, ev_geo, gaps, T_SAT_COND)
        res['carry_penalty'] = penalty

    # combo 결과에 대표 페널티 적용 (Case 1 기준)
    if combined:
        for cb in combos:
            apply_carry_penalty(cb['results'], combined[0].get('carry_penalty', np.zeros(len(gaps))))

    # ══════════════════════════════════════════════════════════
    #  모듈 C: 압력강하
    # ══════════════════════════════════════════════════════════
    print("  [C] 압력강하...")
    inlet = InletCondition(T_in=27.0, RH_in=0.70, CMM=CMM_REF, T_wall_evap=T_SAT_EVAP)
    dp_results = {
        f"FT(S)/FT(S)": sweep_dp(gaps, evap_s, cond_s, geo_es, geo_cs, inlet),
        f"FT(I)/FT(I)": sweep_dp(gaps, evap_i, cond_i, geo_ei, geo_ci, inlet),
    }

    # ══════════════════════════════════════════════════════════
    #  가시화
    # ══════════════════════════════════════════════════════════
    print("  [V] Figure 생성...")
    prefix = fin_key

    fig_a = make_module_a_figure(combos, gaps, ref, CMM_REF)
    fig_a.suptitle(
        f"모듈 A — {label}  │  냉방성능 유지율 × Gap × 조합",
        fontsize=13, color="white", fontweight='bold', y=0.975)
    out_a = f"{prefix}_module_a.png"
    fig_a.savefig(out_a, dpi=150, bbox_inches='tight', facecolor=fig_a.get_facecolor())
    import matplotlib.pyplot as plt; plt.close(fig_a)

    fig_b = make_module_b_figure(combined, case_colors, gaps, ref, CMM_REF)
    fig_b.suptitle(
        f"모듈 B — {label}  │  비말동반 위험도 × Carryover Flux",
        fontsize=13, color="white", fontweight='bold', y=0.975)
    out_b = f"{prefix}_module_b.png"
    fig_b.savefig(out_b, dpi=150, bbox_inches='tight', facecolor=fig_b.get_facecolor())
    plt.close(fig_b)

    fig_c = make_module_c_figure(dp_results, gaps, inlet, ref, CMM_REF)
    fig_c.suptitle(
        f"모듈 C — {label}  │  공기측 압력강하 × f-factor × Gap",
        fontsize=13, color="white", fontweight='bold', y=0.975)
    out_c = f"{prefix}_module_c.png"
    fig_c.savefig(out_c, dpi=150, bbox_inches='tight', facecolor=fig_c.get_facecolor())
    plt.close(fig_c)

    print(f"  → {out_a}, {out_b}, {out_c}")

    # ── 콘솔 요약 ────────────────────────────────────────────
    print(f"\n  Cap Retention (semi, G=20/50/100mm):")
    for cb in combos[:3]:
        caps = [r['cap_ratio'] for r in cb['results']]
        caps_c = [r.get('cap_corrected', r['cap_ratio']) for r in cb['results']]
        for g_target in [20, 50, 100]:
            idx = np.argmin(np.abs(gaps - g_target))
            diff = caps[idx] - caps_c[idx]
            carry_note = f" (Cap*={caps_c[idx]:.1f}%, Δ={diff:.1f})" if diff > 0.05 else ""
            print(f"    {cb['label']:<24} G={g_target:>3}mm → {caps[idx]:>5.1f}%{carry_note}")

    dp0 = dp_results[f"FT(S)/FT(S)"]
    idx20 = np.argmin(np.abs(gaps - 20))
    print(f"\n  ΔP (Staggered, G=20mm): {dp0['dp_total'][idx20]:.1f} Pa"
          f"  f_dry={dp0['f_evap_dry']:.4f}  f_wet={dp0['f_evap_wet']:.4f}")


# ═══════════════════════════════════════════════════════════════════
#  메인
# ═══════════════════════════════════════════════════════════════════

def main():
    print("냉매 물성 계산 중...")
    ref = RefrigerantState(refrigerant=REFRIGERANT,
                           T_sat_evap=T_SAT_EVAP, T_sat_cond=T_SAT_COND)
    gaps = GAPS

    for fin_key, fin_cfg in FIN_TYPES.items():
        run_single_fin(fin_key, fin_cfg, ref, gaps)

    print(f"\n\n{'='*70}")
    print(f"  전체 완료: {len(FIN_TYPES)} 핀 × 3 모듈 = {len(FIN_TYPES)*3} Figure 생성")
    print(f"{'='*70}")


if __name__ == "__main__":
    main()
