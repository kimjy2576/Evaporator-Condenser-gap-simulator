"""
run_single.py — 단일 설계 조건 Gap 해석 프로그램

gap_sim_config.py에서 생성한 JSON 또는 Python 파일을 읽어
선택된 1개 케이스에 대해 Gap에 따른 Module A/B/C 영향을 분석

사용법:
  python run_single.py                         # 기본 설정 사용
  python run_single.py gap_config.json         # JSON 설정 파일
  python run_single.py config_out.py           # Python 설정 파일 (향후)

출력: module_a.png, module_b.png, module_c.png + 콘솔 요약
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import os, sys, json, warnings
warnings.filterwarnings("ignore")

from common import (
    FinTubeSpec, MCHXSpec, RefrigerantState,
    compute_ft_geometry, compute_mchx_geometry, compute_UA, C,
    DARK, FLUX_CAUTION, FLUX_DANGER, FLUX_SEVERE,
)
from module_a import GapParams, simulate_gap, sweep
from module_b import analyze_combined
from module_c import InletCondition, sweep_dp

_NOTO = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
if os.path.exists(_NOTO):
    from matplotlib import font_manager as fm
    fm.fontManager.addfont(_NOTO)
    matplotlib.rcParams['font.family'] = fm.FontProperties(fname=_NOTO).get_name()
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['mathtext.fontset'] = 'cm'

# ── 색상 ──
BG    = "#0a0e17"
from visualize import make_single_fig_a, make_single_fig_b, make_single_fig_c


# ═══════════════════════════════════════════════════════════════
#  1. 설정 파일 읽기
# ═══════════════════════════════════════════════════════════════

def load_config(path=None):
    """JSON 또는 Python 설정 파일 로드. 없으면 기본값 사용."""
    if path and path.endswith('.json') and os.path.exists(path):
        print(f"  Config: {path}")
        with open(path) as f:
            cfg = json.load(f)
        return cfg

    elif path and path.endswith('.py') and os.path.exists(path):
        # Python config에서 직접 import
        print(f"  Config: {path}")
        import importlib.util
        spec = importlib.util.spec_from_file_location("cfg", path)
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        # 모듈 속성 → dict 변환
        cfg = {}
        if hasattr(mod, 'evap_spec'):
            cfg['evap'] = _spec_to_dict(mod.evap_spec)
        if hasattr(mod, 'cond_spec'):
            cfg['cond'] = _spec_to_dict(mod.cond_spec)
        if hasattr(mod, 'ref'):
            cfg['ref'] = dict(refrigerant=mod.ref.refrigerant,
                              T_sat_evap=mod.ref.T_sat_evap,
                              T_sat_cond=mod.ref.T_sat_cond)
        if hasattr(mod, 'gp'):
            gp = mod.gp
            cfg['gap'] = dict(T_amb=gp.T_amb, RH_in=gp.RH_in, CMM=gp.CMM,
                              gap_mode=gp.gap_mode, seal_fraction=gp.seal_fraction,
                              mode=gp.mode, frame_material=gp.frame_material,
                              A_frame=gp.A_frame*1e4)
        if hasattr(mod, 'GAPS'):
            cfg['sim'] = dict(gap_min=float(mod.GAPS[0]),
                              gap_max=float(mod.GAPS[-1]),
                              gap_points=len(mod.GAPS))
        return cfg

    else:
        print("  Config: defaults")
        return None


def _spec_to_dict(sp):
    """FinTubeSpec → dict (mm 단위)"""
    return dict(
        W=sp.W, H=sp.H, D=sp.D,
        fin_pitch=sp.fin_pitch*1e3, fin_thickness=sp.fin_thickness*1e3,
        fin_height=sp.fin_height*1e3,
        tube_do=sp.tube_do*1e3, tube_di=sp.tube_di*1e3,
        tube_rows=sp.tube_rows, tube_cols=sp.tube_cols,
        tube_pitch_t=sp.tube_pitch_t*1e3, tube_pitch_l=sp.tube_pitch_l*1e3,
        fin_type=sp.fin_type, corr_mode=getattr(sp,'corr_mode','etype'),
        tube_layout=sp.tube_layout,
        wavy_height=getattr(sp,'wavy_height',1.5e-3)*1e3,
        wavy_angle=getattr(sp,'wavy_angle',17.5),
        ft_louver_pitch=getattr(sp,'ft_louver_pitch',1.7e-3)*1e3,
        ft_louver_angle=getattr(sp,'ft_louver_angle',28.0),
        slit_num=getattr(sp,'slit_num',6),
        slit_height=getattr(sp,'slit_height',1.0e-3)*1e3,
        fin_material=sp.fin_material, tube_material=sp.tube_material,
        h_drain=sp.h_drain*1e3,
    )


def build_spec(d):
    """dict (mm 단위) → FinTubeSpec. fpi>0이면 fin_pitch 자동 계산.
    fin_height는 tube_pitch_t로 자동 설정 (레거시 호환)."""
    fpi = d.get('fpi', 0)
    Pt = d.get('tube_pitch_t', 25.4)
    return FinTubeSpec(
        W=d['W'], H=d['H'], D=d['D'],
        fin_pitch=d.get('fin_pitch', 1.8)*1e-3,
        fpi=fpi, fin_thickness=d['fin_thickness']*1e-3,
        fin_height=Pt*1e-3,  # = tube_pitch_t (자동)
        tube_do=d['tube_do']*1e-3, tube_di=d['tube_di']*1e-3,
        tube_rows=int(d['tube_rows']), tube_cols=int(d['tube_cols']),
        tube_pitch_t=d['tube_pitch_t']*1e-3, tube_pitch_l=d['tube_pitch_l']*1e-3,
        fin_type=d.get('fin_type','plain'),
        corr_mode=d.get('corr_mode','etype'),
        tube_layout=d.get('tube_layout','staggered'),
        wavy_height=d.get('wavy_height',1.5)*1e-3,
        wavy_angle=d.get('wavy_angle',17.5),
        ft_louver_pitch=d.get('ft_louver_pitch',1.7)*1e-3,
        ft_louver_angle=d.get('ft_louver_angle',28.0),
        slit_num=int(d.get('slit_num',6)),
        slit_height=d.get('slit_height',1.0)*1e-3,
        fin_material=d.get('fin_material','Al'),
        tube_material=d.get('tube_material','Cu'),
        h_drain=d.get('h_drain',25.0)*1e-3,
    )


# ═══════════════════════════════════════════════════════════════
#  2. 해석 실행
# ═══════════════════════════════════════════════════════════════

def run_analysis(cfg=None):
    """단일 설계 조건 Gap 해석"""

    # 기본값
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

    # ── 스펙 빌드 ──
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

    # ── 전처리 ──
    geo_e = compute_ft_geometry(evap)
    geo_c = compute_ft_geometry(cond)
    V_face = gp.CMM / (60.0 * evap.W * evap.H)
    ua_e = compute_UA(evap, geo_e, ref, V_face, 'evap')
    ua_c = compute_UA(cond, geo_c, ref, V_face, 'cond')

    tag = (f"{evap.fin_type.capitalize()} {evap.tube_layout[0].upper()} | "
           f"{ref.refrigerant} | T_evap={ref.T_sat_evap}°C | "
           f"CMM={gp.CMM} | {gp.gap_mode}")

    print(f"\n{'='*70}")
    print(f"  Single-Case Gap Analysis")
    print(f"  {tag}")
    print(f"{'='*70}")
    print(f"  Evap: {evap.fin_type} {evap.tube_layout} {evap.tube_rows}R×{evap.tube_cols}C"
          f"  UA={ua_e['UA']:.1f} W/K  h_o={ua_e['h_o']:.1f} W/m²K")
    print(f"  Cond: {cond.fin_type} {cond.tube_layout} {cond.tube_rows}R×{cond.tube_cols}C"
          f"  UA={ua_c['UA']:.1f} W/K  h_o={ua_c['h_o']:.1f} W/m²K")
    print(f"  Gap:  {sim_d['gap_min']}~{sim_d['gap_max']}mm ({sim_d['gap_points']}pts)"
          f"  mode={gp.gap_mode}  RH={gp.RH_in:.0%}  T_amb={gp.T_amb}°C")

    # ══════════════════════════════════════════════════════════
    #  Module A: 냉방성능 Gap sweep
    # ══════════════════════════════════════════════════════════
    print("\n  [A] Cooling performance...")
    results_a = sweep(gaps, evap, cond, geo_e, geo_c, ua_e, ua_c, ref, gp)

    # ══════════════════════════════════════════════════════════
    #  Module B: 비말동반
    # ══════════════════════════════════════════════════════════
    print("  [B] Carryover risk...")
    case_b = dict(
        name=tag, T_in=gap_d['T_amb'], RH_in=gap_d['RH_in'],
        CMM=gap_d['CMM'], T_wall=ref.T_sat_evap,
        gap_mm=20, theta_deg=0.0, eta_coeff=5e-4, w_bridge=0.75,
        gap_mode=gap_d['gap_mode'], seal_fraction=gap_d['seal_fraction'],
        hx_type="FT", tube_layout=evap.tube_layout,
    )
    result_b = analyze_combined(case_b, gaps, evap, geo_e, cond, geo_c,
                                ua_e, ua_c, ref)

    # ══════════════════════════════════════════════════════════
    #  A+B 통합: 비말동반 → 응축기 잠열 페널티
    # ══════════════════════════════════════════════════════════
    from module_b import compute_carry_penalty, apply_carry_penalty
    carry_penalty = compute_carry_penalty(result_b, evap, geo_e, gaps, ref.T_sat_cond)
    apply_carry_penalty(results_a, carry_penalty)
    print(f"    Carry penalty @ G=20: {carry_penalty[np.argmin(np.abs(gaps-20))]:.1f}W")

    # ══════════════════════════════════════════════════════════
    #  Module C: 압력강하
    # ══════════════════════════════════════════════════════════
    print("  [C] Pressure drop...")
    inlet = InletCondition(T_in=gap_d['T_amb'], RH_in=gap_d['RH_in'],
                           CMM=gap_d['CMM'], T_wall_evap=ref.T_sat_evap)
    result_c = sweep_dp(gaps, evap, cond, geo_e, geo_c, inlet)

    # ══════════════════════════════════════════════════════════
    #  3. 가시화
    # ══════════════════════════════════════════════════════════
    print("\n  Generating figures...")

    # Module B 가시화 추가 데이터 계산
    from module_b import CarryoverSpec, we_ch, we_crit, v_onset, monte_carlo, eta_co
    co_spec = CarryoverSpec(evap, geo_e)
    We_ch_val = we_ch(V_face, co_spec)
    We_crit_val = we_crit(co_spec)
    V_onset_val = v_onset(co_spec)

    # P_reach vs Gap
    P_reach_arr = np.array([
        monte_carlo(g, V_face, co_spec, 0, 0.75, N=40, seed=42)['P_reach']
        for g in gaps])

    # We ratio vs V_face sweep
    V_sweep = np.linspace(0.3, max(V_face*2.5, 6.0), 60)
    We_sweep = np.array([we_ch(v, co_spec) / We_crit_val for v in V_sweep])

    viz_b = dict(
        P_reach_arr=P_reach_arr,
        We_ratio=We_ch_val / (We_crit_val + 1e-9),
        V_face=V_face,
        V_onset=V_onset_val,
        q_cond=result_b['cr']['q_cond'],
        carry_penalty=carry_penalty,
        V_sweep=V_sweep,
        We_sweep=We_sweep,
        gaps=gaps,
    )

    make_single_fig_a(gaps, results_a, gp, ref, tag)
    make_single_fig_b(gaps, result_b, viz_b, tag)
    make_single_fig_c(gaps, result_c, inlet, tag)

    # Gap 물리 현상 개략도 (대표 Gap = 20mm 또는 중간값)
    from visualize import make_single_schematic
    g_repr = 20 if 20 <= gaps[-1] else gaps[len(gaps)//2]
    idx_repr = np.argmin(np.abs(gaps - g_repr))
    make_single_schematic(evap, cond, geo_e, geo_c, ua_e, ua_c,
                          gp, ref, results_a[idx_repr], viz_b, tag)

    # ══════════════════════════════════════════════════════════
    #  4. CSV 결과 파일 생성
    # ══════════════════════════════════════════════════════════
    import csv
    print("\n  Exporting CSV...")

    # ── Module A: 냉방성능 ──
    with open('module_a_results.csv', 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow([
            'Gap_mm',
            # 성능 지표
            'Cap_Retention_%', 'Cap_Corrected_%',
            'Q_net_W', 'Q_total_W', 'Q_useful_W',
            'Q_sensible_W', 'Q_latent_W', 'Q_evap_W',
            'SHR',
            'Dehumid_rate_g_per_h',
            # 손실 분해
            'Q_rad_W', 'Q_cond_W', 'Q_mix_W', 'Q_carry_W',
            'Q_recir_waste_W',
            # 온도/습도
            'T_in_eff_C', 'T_evap_out_C', 'T_dp_C',
            'dT_recir_C',
            'W_in_kg_per_kg', 'W_out_kg_per_kg',
            # 열교환 파라미터
            'UA_evap_W_per_K', 'NTU', 'eps',
            'eta_mix',
            # 운전 조건
            'V_face_m_per_s', 'CMM', 'gap_mode',
            # Segment 정보
            'wet_rows', 'Nr', 'q_cond_g_per_s',
        ])
        for i, r in enumerate(results_a):
            w.writerow([
                f"{gaps[i]:.1f}",
                f"{r['cap_ratio']:.2f}", f"{r.get('cap_corrected', r['cap_ratio']):.2f}",
                f"{r['Q_net']:.1f}", f"{r['Q_total']:.1f}", f"{r['Q_useful']:.1f}",
                f"{r['Q_sensible']:.1f}", f"{r['Q_latent']:.1f}", f"{r['Q_evap']:.1f}",
                f"{r['SHR']:.4f}",
                f"{r['dehumid_rate']:.1f}",
                f"{r['Q_rad']:.2f}", f"{r['Q_cond']:.2f}", f"{r['Q_mix']:.2f}",
                f"{r.get('Q_carry_penalty', 0):.2f}",
                f"{r['Q_recir_waste']:.2f}",
                f"{r['T_in_eff']:.2f}", f"{r['T_evap_out']:.2f}", f"{r['T_dp']:.2f}",
                f"{r['dT_recir']:.3f}",
                f"{r['W_in']:.6f}", f"{r['W_out']:.6f}",
                f"{r['UA_evap']:.2f}", f"{r['NTU']:.4f}", f"{r['eps']:.4f}",
                f"{r['eta_mix']:.5f}",
                f"{r['V_face']:.3f}", f"{r['CMM']:.2f}", r['gap_mode'],
                f"{r.get('wet_rows', '')}", f"{r.get('Nr', '')}", f"{r.get('q_cond', 0):.3f}",
            ])
    print("    → module_a_results.csv")

    # ── Module B: 비말동반 ──
    with open('module_b_results.csv', 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow([
            'Gap_mm',
            # 비말동반 지표
            'Carryover_Flux',
            'P_reach',
            'Risk_Level',
            # Module A 연계
            'Cap_Retention_%',
            # 응축기 페널티
            'Q_carry_penalty_W',
        ])
        for i, g in enumerate(gaps):
            flux = result_b['flux_arr'][i]
            risk = 'SAFE' if flux < FLUX_CAUTION else ('CAUTION' if flux < FLUX_DANGER else ('DANGER' if flux < FLUX_SEVERE else 'SEVERE'))
            w.writerow([
                f"{g:.1f}",
                f"{flux:.4f}",
                f"{viz_b['P_reach_arr'][i]:.4f}",
                risk,
                f"{result_b['cap_arr'][i]:.2f}",
                f"{carry_penalty[i]:.2f}",
            ])
        # 하단에 고정 파라미터 요약 행
        w.writerow([])
        w.writerow(['# Fixed Parameters'])
        w.writerow(['q_cond_g_per_s', f"{viz_b['q_cond']:.3f}"])
        w.writerow(['V_face_m_per_s', f"{viz_b['V_face']:.3f}"])
        w.writerow(['V_onset_m_per_s', f"{viz_b['V_onset']:.3f}"])
        w.writerow(['We_ratio', f"{viz_b['We_ratio']:.4f}"])
        w.writerow(['eta_co', f"{eta_co(V_face, co_spec, 5e-4) if V_face > V_onset_val else 0:.6f}"])
        w.writerow(['Sweet_gap_mm', f"{result_b['sweet_gap']}" if result_b['sweet_gap'] else 'N/A'])
        w.writerow(['cr_Q_total_W', f"{result_b['cr']['Q_total']:.1f}"])
        w.writerow(['cr_Q_latent_W', f"{result_b['cr']['Q_lat']:.1f}"])
        w.writerow(['cr_SHR', f"{result_b['cr']['SHR']:.3f}"])
        w.writerow(['cr_wet_fraction', f"{result_b['cr']['wet_fraction']:.2f}"])
    print("    → module_b_results.csv")

    # ── Module C: 압력강하 ──
    with open('module_c_results.csv', 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        w.writerow([
            'Gap_mm',
            # 압력강하 분해
            'dP_total_Pa', 'dP_evap_Pa', 'dP_gap_Pa', 'dP_cond_Pa',
            # 비율
            'dP_gap_fraction_%', 'dP_evap_fraction_%', 'dP_cond_fraction_%',
        ])
        for i, g in enumerate(gaps):
            dp_t = result_c['dp_total'][i]
            dp_e = result_c['dp_evap'][i]
            dp_g = result_c['dp_gap'][i]
            dp_cc = result_c['dp_cond'][i]
            w.writerow([
                f"{g:.1f}",
                f"{dp_t:.2f}", f"{dp_e:.2f}", f"{dp_g:.2f}", f"{dp_cc:.2f}",
                f"{dp_g/(dp_t+1e-9)*100:.1f}",
                f"{dp_e/(dp_t+1e-9)*100:.1f}",
                f"{dp_cc/(dp_t+1e-9)*100:.1f}",
            ])
        # 고정 파라미터
        w.writerow([])
        w.writerow(['# Fixed Parameters'])
        w.writerow(['V_face_m_per_s', f"{result_c['V_face']:.3f}"])
        w.writerow(['f_evap_dry', f"{result_c['f_evap_dry']:.5f}"])
        w.writerow(['f_evap_wet', f"{result_c['f_evap_wet']:.5f}"])
        w.writerow(['f_wet_ratio', f"{result_c.get('f_wet_ratio', result_c['f_evap_wet']/(result_c['f_evap_dry']+1e-9)):.4f}"])
        w.writerow(['f_cond', f"{result_c['f_cond']:.5f}"])
        w.writerow(['evap_layout', result_c.get('evap_layout', '')])
        w.writerow(['cond_layout', result_c.get('cond_layout', '')])
    print("    → module_c_results.csv")

    # ── 통합 요약 CSV (전 모듈 주요 지표 한 파일) ──
    with open('gap_analysis_summary.csv', 'w', encoding='utf-8', newline='') as f:
        w = csv.writer(f, lineterminator='\n')
        # 헤더
        w.writerow([
            'Gap_mm',
            # A
            'Cap_%', 'Cap*_%', 'Q_net_W', 'Q_sensible_W', 'Q_latent_W',
            'SHR', 'Dehumid_g_per_h', 'dT_recir_C',
            'Q_rad_W', 'Q_cond_W', 'Q_mix_W', 'Q_carry_W',
            # B
            'Flux', 'P_reach', 'Risk',
            # C
            'dP_total_Pa', 'dP_evap_Pa', 'dP_gap_Pa', 'dP_cond_Pa', 'dP_gap_%',
        ])
        for i, g in enumerate(gaps):
            r = results_a[i]
            flux = result_b['flux_arr'][i]
            risk = 'SAFE' if flux < FLUX_CAUTION else ('CAUTION' if flux < FLUX_DANGER else ('DANGER' if flux < FLUX_SEVERE else 'SEVERE'))
            w.writerow([
                f"{g:.1f}",
                f"{r['cap_ratio']:.2f}", f"{r.get('cap_corrected', r['cap_ratio']):.2f}",
                f"{r['Q_net']:.1f}", f"{r['Q_sensible']:.1f}", f"{r['Q_latent']:.1f}",
                f"{r['SHR']:.4f}", f"{r['dehumid_rate']:.1f}", f"{r['dT_recir']:.3f}",
                f"{r['Q_rad']:.2f}", f"{r['Q_cond']:.2f}", f"{r['Q_mix']:.2f}",
                f"{r.get('Q_carry_penalty', 0):.2f}",
                f"{flux:.4f}", f"{viz_b['P_reach_arr'][i]:.4f}", risk,
                f"{result_c['dp_total'][i]:.2f}", f"{result_c['dp_evap'][i]:.2f}",
                f"{result_c['dp_gap'][i]:.2f}", f"{result_c['dp_cond'][i]:.2f}",
                f"{result_c['dp_gap'][i]/(result_c['dp_total'][i]+1e-9)*100:.1f}",
            ])
        # 설계 조건 요약
        w.writerow([])
        w.writerow(['# Design Conditions'])
        w.writerow(['Refrigerant', ref.refrigerant])
        w.writerow(['T_sat_evap_C', ref.T_sat_evap])
        w.writerow(['T_sat_cond_C', ref.T_sat_cond])
        w.writerow(['T_amb_C', gap_d['T_amb']])
        w.writerow(['RH_in', gap_d['RH_in']])
        w.writerow(['CMM', gap_d['CMM']])
        w.writerow(['gap_mode', gap_d['gap_mode']])
        w.writerow(['seal_fraction', gap_d['seal_fraction']])
        w.writerow(['Evap_fin_type', evap.fin_type])
        w.writerow(['Evap_layout', evap.tube_layout])
        w.writerow(['Evap_FPI', f"{25.4/(evap.fin_pitch*1e3):.1f}" if evap.fin_pitch > 0 else 'N/A'])
        w.writerow(['Evap_UA_W_per_K', f"{ua_e['UA']:.2f}"])
        w.writerow(['Cond_fin_type', cond.fin_type])
        w.writerow(['Cond_UA_W_per_K', f"{ua_c['UA']:.2f}"])
        w.writerow(['V_face_m_per_s', f"{V_face:.3f}"])
        w.writerow(['V_onset_m_per_s', f"{V_onset_val:.3f}"])
        w.writerow(['f_dry', f"{result_c['f_evap_dry']:.5f}"])
        w.writerow(['f_wet', f"{result_c['f_evap_wet']:.5f}"])
        w.writerow(['f_wet_ratio', f"{result_c.get('f_wet_ratio', 0):.4f}"])
    print("    → gap_analysis_summary.csv")

    # ── 콘솔 요약 ──
    print(f"\n  {'Gap':>5} | {'Cap':>6} {'Cap*':>6} | {'Q_net':>7} {'Q_lat':>6} {'SHR':>5} | "
          f"{'Flux':>6} {'Q_crr':>6} {'Risk':>8} | {'dP_tot':>7}")
    print(f"  {'-'*82}")
    for g in [5, 10, 20, 30, 50, 75, 100]:
        if g > gaps[-1]:
            break
        idx = np.argmin(np.abs(gaps - g))
        r = results_a[idx]
        flux = result_b['flux_arr'][idx]
        risk = 'SAFE' if flux < FLUX_CAUTION else ('CAUTION' if flux < FLUX_DANGER else ('DANGER' if flux < FLUX_SEVERE else 'SEVERE'))
        dp = result_c['dp_total'][idx]
        qc = r.get('Q_carry_penalty', 0)
        cap_c = r.get('cap_corrected', r['cap_ratio'])
        print(f"  {g:>5} | {r['cap_ratio']:>5.1f}% {cap_c:>5.1f}% | {r['Q_net']:>6.0f}W {r['Q_latent']:>5.0f}W"
              f" {r['SHR']:>4.2f} | {flux:>5.3f} {qc:>5.1f}W {risk:>8} | {dp:>6.1f}Pa")

    print(f"\n  → module_a.png, module_b.png, module_c.png")
    return results_a, result_b, result_c


#  Entry point
# ═══════════════════════════════════════════════════════════════

if __name__ == '__main__':
    path = sys.argv[1] if len(sys.argv) > 1 else None
    cfg = load_config(path)
    run_analysis(cfg)
