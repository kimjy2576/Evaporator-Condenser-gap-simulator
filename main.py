"""
main.py — 메인 실행 파일
전처리(입력 변수) → 해석(모듈 A + B + C) → 후처리(가시화)

[실행 흐름]
  1. 전처리: HX 스펙 정의, 냉매 조건, 입구 조건(CMM), 케이스 설정
  2. 해석:   모듈 A (냉방성능) + 모듈 B (비말동반) + 모듈 C (압력강하)
  3. 후처리: 모듈별 독립 Figure 생성 및 콘솔 요약
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
#  1. 전처리 — 입력 변수 정의
# ═══════════════════════════════════════════════════════════════════

# ── 냉매 조건 ─────────────────────────────────────────────────────
REFRIGERANT  = "R410A"
T_SAT_EVAP   = 5.0    # 증발 포화온도 [°C]
T_SAT_COND   = 45.0   # 응축 포화온도 [°C]

# ── 풍량 기준 ─────────────────────────────────────────────────────
CMM_REF = 6.75         # 기준 풍량 [m³/min]

# ── Gap sweep 범위 ───────────────────────────────────────────────
GAPS = np.linspace(1, 100, 30)

# ── 비말동반 분석 케이스 ──────────────────────────────────────────
CASES = [
    # FT Staggered
    dict(name="Case 1 — FT(S) 기준",
         T_in=27.0, RH_in=0.70, CMM=13.5, T_wall=5.0,
         gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
         gap_mode="semi", seal_fraction=0.7,
         hx_type="FT", tube_layout="staggered"),
    dict(name="Case 2 — FT(S) 고풍속/고온고습",
         T_in=30.0, RH_in=0.85, CMM=18.0, T_wall=5.0,
         gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
         gap_mode="semi", seal_fraction=0.7,
         hx_type="FT", tube_layout="staggered"),
    dict(name="Case 3 — FT(S) 친수코팅/θ=15°",
         T_in=27.0, RH_in=0.70, CMM=15.75, T_wall=8.0,
         gap_mm=50, theta_deg=15.0, eta_coeff=1e-4, w_bridge=0.60,
         gap_mode="semi", seal_fraction=0.7,
         hx_type="FT", tube_layout="staggered"),
    # FT Inline
    dict(name="Case 4 — FT(I) 기준",
         T_in=27.0, RH_in=0.70, CMM=13.5, T_wall=5.0,
         gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
         gap_mode="semi", seal_fraction=0.7,
         hx_type="FT", tube_layout="inline"),
    dict(name="Case 5 — FT(I) 고풍속/고온고습",
         T_in=30.0, RH_in=0.85, CMM=18.0, T_wall=5.0,
         gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
         gap_mode="semi", seal_fraction=0.7,
         hx_type="FT", tube_layout="inline"),
    # MCHX
    dict(name="Case 6 — MCHX 기준 (Lp=1.2mm)",
         T_in=27.0, RH_in=0.70, CMM=13.5, T_wall=5.0,
         gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
         gap_mode="semi", seal_fraction=0.7,
         hx_type="MCHX"),
    dict(name="Case 7 — MCHX 고밀도핀 (Lp=0.9mm)",
         T_in=30.0, RH_in=0.85, CMM=18.0, T_wall=5.0,
         gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.80,
         gap_mode="semi", seal_fraction=0.7,
         hx_type="MCHX", louver_pitch_override=0.9e-3),
    dict(name="Case 8 — MCHX 친수코팅",
         T_in=27.0, RH_in=0.70, CMM=15.75, T_wall=5.0,
         gap_mm=20, theta_deg=0.0, eta_coeff=1e-4, w_bridge=0.85,
         gap_mode="semi", seal_fraction=0.7,
         hx_type="MCHX"),
]


# ═══════════════════════════════════════════════════════════════════
#  2. 해석 — 모듈 A + B + C
# ═══════════════════════════════════════════════════════════════════

def run_analysis():
    gaps = GAPS

    # ── 전처리: 냉매 물성 ────────────────────────────────────
    print("전처리: 냉매 물성 계산 중 (CoolProp)...")
    ref = RefrigerantState(refrigerant=REFRIGERANT,
                           T_sat_evap=T_SAT_EVAP, T_sat_cond=T_SAT_COND)

    # ── 전처리: HX 스펙 정의 ─────────────────────────────────
    evap_ft = FinTubeSpec(
        W=0.30, H=0.25, D=0.045,
        fin_pitch=1.8e-3, fin_thickness=0.10e-3, fin_height=25.4e-3,
        tube_do=9.52e-3, tube_di=8.52e-3, tube_rows=2, tube_cols=8,
        tube_pitch_t=25.4e-3, tube_pitch_l=22.0e-3,
        tube_layout="staggered",
        fin_material="Al", tube_material="Cu", h_drain=0.025)
    evap_ft_inline = FinTubeSpec(
        W=0.30, H=0.25, D=0.045,
        fin_pitch=1.8e-3, fin_thickness=0.10e-3, fin_height=25.4e-3,
        tube_do=9.52e-3, tube_di=8.52e-3, tube_rows=2, tube_cols=8,
        tube_pitch_t=25.4e-3, tube_pitch_l=22.0e-3,
        tube_layout="inline",
        fin_material="Al", tube_material="Cu", h_drain=0.025)
    evap_mchx = MCHXSpec(
        W=0.30, H=0.20, D=0.020,
        ch_width=0.8e-3, ch_height=1.5e-3, ch_wall=0.4e-3, n_ports=12,
        fin_pitch=1.0e-3, fin_thickness=0.06e-3,
        louver_pitch=1.2e-3, louver_angle=27.0, n_slabs=1, slab_pitch=8.0e-3)
    cond_ft = FinTubeSpec(
        W=0.35, H=0.28, D=0.050,
        fin_pitch=1.5e-3, fin_thickness=0.10e-3, fin_height=25.4e-3,
        tube_do=9.52e-3, tube_di=8.52e-3, tube_rows=2, tube_cols=9,
        tube_pitch_t=25.4e-3, tube_pitch_l=22.0e-3,
        tube_layout="staggered",
        fin_material="Al", tube_material="Cu", h_drain=0.025)
    cond_ft_inline = FinTubeSpec(
        W=0.35, H=0.28, D=0.050,
        fin_pitch=1.5e-3, fin_thickness=0.10e-3, fin_height=25.4e-3,
        tube_do=9.52e-3, tube_di=8.52e-3, tube_rows=2, tube_cols=9,
        tube_pitch_t=25.4e-3, tube_pitch_l=22.0e-3,
        tube_layout="inline",
        fin_material="Al", tube_material="Cu", h_drain=0.025)
    cond_mchx = MCHXSpec(
        W=0.35, H=0.25, D=0.018,
        ch_width=1.0e-3, ch_height=2.0e-3, ch_wall=0.5e-3, n_ports=14,
        fin_pitch=1.2e-3, fin_thickness=0.07e-3,
        louver_pitch=1.4e-3, louver_angle=25.0, n_slabs=1, slab_pitch=10.0e-3)

    # ── 전처리: 기하 + UA ────────────────────────────────────
    geo_ev_ft     = compute_ft_geometry(evap_ft)
    geo_ev_ft_inl = compute_ft_geometry(evap_ft_inline)
    geo_ev_mchx   = compute_mchx_geometry(evap_mchx)
    geo_cd_ft     = compute_ft_geometry(cond_ft)
    geo_cd_ft_inl = compute_ft_geometry(cond_ft_inline)
    geo_cd_mchx   = compute_mchx_geometry(cond_mchx)

    Vref_ft   = CMM_REF / (60.0 * evap_ft.W * evap_ft.H)
    Vref_ft_i = CMM_REF / (60.0 * evap_ft_inline.W * evap_ft_inline.H)
    Vref_mchx = CMM_REF / (60.0 * evap_mchx.W * evap_mchx.H)

    ua_ev_ft     = compute_UA(evap_ft,        geo_ev_ft,     ref, Vref_ft,   'evap')
    ua_ev_ft_inl = compute_UA(evap_ft_inline, geo_ev_ft_inl, ref, Vref_ft_i, 'evap')
    ua_ev_mchx   = compute_UA(evap_mchx,      geo_ev_mchx,  ref, Vref_mchx, 'evap')
    ua_cd_ft     = compute_UA(cond_ft,        geo_cd_ft,     ref, Vref_ft,   'cond')
    ua_cd_ft_inl = compute_UA(cond_ft_inline, geo_cd_ft_inl, ref, Vref_ft_i, 'cond')
    ua_cd_mchx   = compute_UA(cond_mchx,      geo_cd_mchx,  ref, Vref_mchx, 'cond')

    print(f"  UA(evap FT(S))={ua_ev_ft['UA']:.1f}  UA(evap FT(I))={ua_ev_ft_inl['UA']:.1f}"
          f"  UA(evap MCHX)={ua_ev_mchx['UA']:.1f}")

    # ══════════════════════════════════════════════════════════
    #  모듈 A: 냉방성능 Gap sweep (6조합)
    # ══════════════════════════════════════════════════════════
    print("\\n[모듈 A] 냉방성능 계산 중...")
    combo_defs = [
        dict(label="FT(S)/FT(S) open", color=C[0], ls="-",
             evap_spec=evap_ft, evap_geo=geo_ev_ft, ua_ev=ua_ev_ft,
             cond_spec=cond_ft, cond_geo=geo_cd_ft, ua_cd=ua_cd_ft,
             gp=GapParams(CMM=CMM_REF, mode="forced", gap_mode="open")),
        dict(label="FT(S)/FT(S) semi", color=C[0], ls="--",
             evap_spec=evap_ft, evap_geo=geo_ev_ft, ua_ev=ua_ev_ft,
             cond_spec=cond_ft, cond_geo=geo_cd_ft, ua_cd=ua_cd_ft,
             gp=GapParams(CMM=CMM_REF, mode="forced", gap_mode="semi", seal_fraction=0.7)),
        dict(label="FT(S)/FT(S) sealed", color=C[0], ls=":",
             evap_spec=evap_ft, evap_geo=geo_ev_ft, ua_ev=ua_ev_ft,
             cond_spec=cond_ft, cond_geo=geo_cd_ft, ua_cd=ua_cd_ft,
             gp=GapParams(CMM=CMM_REF, mode="forced", gap_mode="sealed")),
        dict(label="FT(I)/FT(I) semi", color=C[3], ls="-.",
             evap_spec=evap_ft_inline, evap_geo=geo_ev_ft_inl, ua_ev=ua_ev_ft_inl,
             cond_spec=cond_ft_inline, cond_geo=geo_cd_ft_inl, ua_cd=ua_cd_ft_inl,
             gp=GapParams(CMM=CMM_REF, mode="forced", gap_mode="semi", seal_fraction=0.7)),
        dict(label="MCHX/MCHX semi", color=C[4], ls="-.",
             evap_spec=evap_mchx, evap_geo=geo_ev_mchx, ua_ev=ua_ev_mchx,
             cond_spec=cond_mchx, cond_geo=geo_cd_mchx, ua_cd=ua_cd_mchx,
             gp=GapParams(CMM=CMM_REF, mode="forced", gap_mode="semi", seal_fraction=0.7)),
        dict(label="FT(S)/MCHX semi", color=C[1], ls="--",
             evap_spec=evap_ft, evap_geo=geo_ev_ft, ua_ev=ua_ev_ft,
             cond_spec=cond_mchx, cond_geo=geo_cd_mchx, ua_cd=ua_cd_mchx,
             gp=GapParams(CMM=CMM_REF, mode="forced", gap_mode="semi", seal_fraction=0.7)),
    ]
    for cb in combo_defs:
        cb['results'] = sweep(gaps, cb['evap_spec'], cb['cond_spec'],
                              cb['evap_geo'], cb['cond_geo'],
                              cb['ua_ev'], cb['ua_cd'], ref, cb['gp'])

    # ══════════════════════════════════════════════════════════
    #  모듈 B: 비말동반 위험도 (8케이스)
    # ══════════════════════════════════════════════════════════
    print("\\n[모듈 B] 비말동반 위험도 계산 중...")
    case_colors = [C[0], C[1], C[2], C[0], C[1], C[3], C[4], C[5]]
    combined_results = []

    for case, cc in zip(CASES, case_colors):
        print(f"  {case['name']}")
        if case.get('hx_type', 'FT') == 'MCHX':
            evap_spec_case = make_mchx_evap(evap_mchx, case)
            evap_geo_case  = (compute_mchx_geometry(evap_spec_case)
                              if 'louver_pitch_override' in case else geo_ev_mchx)
            ua_evap_case   = (compute_UA(evap_spec_case, evap_geo_case, ref, Vref_mchx, 'evap')
                              if 'louver_pitch_override' in case else ua_ev_mchx)
            res = analyze_combined(case, gaps, evap_spec_case, evap_geo_case,
                                   cond_ft, geo_cd_ft, ua_evap_case, ua_cd_ft, ref)
        elif case.get('tube_layout', 'staggered') == 'inline':
            res = analyze_combined(case, gaps, evap_ft_inline, geo_ev_ft_inl,
                                   cond_ft_inline, geo_cd_ft_inl,
                                   ua_ev_ft_inl, ua_cd_ft_inl, ref)
        else:
            res = analyze_combined(case, gaps, evap_ft, geo_ev_ft,
                                   cond_ft, geo_cd_ft, ua_ev_ft, ua_cd_ft, ref)
        combined_results.append(res)

    # ── A+B 통합: 비말동반 → 응축기 잠열 페널티 ──────────────
    print("\\n[A+B] 응축기 잠열 페널티 계산 중...")
    from module_b import compute_carry_penalty, apply_carry_penalty
    for res in combined_results:
        case = res['case']
        hx_type = case.get('hx_type', 'FT')
        if hx_type == 'MCHX':
            ev_sp, ev_geo = evap_mchx, geo_ev_mchx
        elif case.get('tube_layout', 'staggered') == 'inline':
            ev_sp, ev_geo = evap_ft_inline, geo_ev_ft_inl
        else:
            ev_sp, ev_geo = evap_ft, geo_ev_ft
        penalty = compute_carry_penalty(res, ev_sp, ev_geo, gaps, T_SAT_COND)
        res['carry_penalty'] = penalty
        print(f"    {case['name']:<36} Q_carry(G=20)={penalty[np.argmin(np.abs(gaps-20))]:.1f}W")

    # 각 combo_def의 Module A 결과에도 대표 페널티 적용
    # (combo_defs는 단일 케이스가 아니므로 Case 1 기준으로 대표값 사용)
    if combined_results:
        for cb in combo_defs:
            apply_carry_penalty(cb['results'], combined_results[0].get('carry_penalty', np.zeros(len(gaps))))

    # ══════════════════════════════════════════════════════════
    #  모듈 C: 압력강하 (고정 CMM)
    # ══════════════════════════════════════════════════════════
    print("\\n[모듈 C] 압력강하 계산 중...")
    inlet = InletCondition(T_in=27.0, RH_in=0.70, CMM=CMM_REF, T_wall_evap=T_SAT_EVAP)
    dp_results = {}
    dp_combos = [
        ("FT(S)/FT(S)", evap_ft, cond_ft, geo_ev_ft, geo_cd_ft),
        ("FT(I)/FT(I)", evap_ft_inline, cond_ft_inline, geo_ev_ft_inl, geo_cd_ft_inl),
        ("MCHX/MCHX",   evap_mchx, cond_mchx, geo_ev_mchx, geo_cd_mchx),
    ]
    for label, es, cs, ge, gc in dp_combos:
        dp_results[label] = sweep_dp(gaps, es, cs, ge, gc, inlet)
        print(f"  {label}: ΔP(G=20)={dp_results[label]['dp_total'][np.argmin(np.abs(gaps-20))]:.1f} Pa")

    # ══════════════════════════════════════════════════════════
    #  3. 후처리 — 가시화 + 콘솔 요약
    # ══════════════════════════════════════════════════════════
    print("\\n[후처리] Figure 생성 중...")
    fig_a = make_module_a_figure(combo_defs, gaps, ref, CMM_REF)
    fig_a.savefig("module_a_cooling_performance.png", dpi=150,
                  bbox_inches='tight', facecolor=fig_a.get_facecolor())
    print("  → module_a_cooling_performance.png")

    fig_b = make_module_b_figure(combined_results, case_colors, gaps, ref, CMM_REF)
    fig_b.savefig("module_b_carryover_risk.png", dpi=150,
                  bbox_inches='tight', facecolor=fig_b.get_facecolor())
    print("  → module_b_carryover_risk.png")

    fig_c = make_module_c_figure(dp_results, gaps, inlet, ref, CMM_REF)
    fig_c.savefig("module_c_pressure_drop.png", dpi=150,
                  bbox_inches='tight', facecolor=fig_c.get_facecolor())
    print("  → module_c_pressure_drop.png")

    # ── 콘솔 요약 ────────────────────────────────────────────
    print(f"\\n{'='*80}")
    print(f"  통합 분석 결과 요약  |  {REFRIGERANT}  T_evap={T_SAT_EVAP}°C  T_cond={T_SAT_COND}°C")
    print(f"{'='*80}")

    print(f"\\n  [모듈 A] Gap별 Cap Retention (CMM={CMM_REF} m³/min)")
    print(f"  {'Gap':>5} | " + " | ".join(f"{cb['label']:>20}" for cb in combo_defs))
    print(f"  {'-'*140}")
    for g in [5, 10, 20, 30, 50, 75, 100]:
        row = f"  {g:>5} | "
        for cb in combo_defs:
            r = simulate_gap(g, cb['evap_spec'], cb['cond_spec'],
                             cb['evap_geo'], cb['cond_geo'],
                             cb['ua_ev'], cb['ua_cd'], ref, cb['gp'])
            row += f"{r['cap_ratio']:>19.1f}% | "
        print(row)

    print(f"\\n  [모듈 B] Sweet-spot 요약")
    for res in combined_results:
        c = res['case']; g = res['sweet_gap']
        hx = res.get('evap_hx_type', 'FT')
        cp20 = res.get('carry_penalty', np.zeros(len(gaps)))[np.argmin(np.abs(gaps-20))]
        if g is not None:
            idx = np.argmin(np.abs(res['gaps'] - g))
            print(f"    {c['name']:<36} G*={g:.0f}mm  Cap={res['cap_arr'][idx]:.1f}%"
                  f"  Flux={res['flux_arr'][idx]:.2f}  Q_carry={cp20:.1f}W")
        else:
            print(f"    {c['name']:<36} N/A (조건 미충족)  Q_carry={cp20:.1f}W")

    print(f"\\n  [모듈 C] 시스템 ΔP (G=20mm)")
    for label, dp in dp_results.items():
        idx20 = np.argmin(np.abs(gaps - 20))
        print(f"    {label:<16} ΔP_total={dp['dp_total'][idx20]:.1f} Pa"
              f"  (evap={dp['dp_evap'][idx20]:.1f} + gap={dp['dp_gap'][idx20]:.1f}"
              f" + cond={dp['dp_cond'][idx20]:.1f})")

    print(f"\\n{'='*80}")
    return combo_defs, combined_results, dp_results, ref


if __name__ == "__main__":
    run_analysis()
