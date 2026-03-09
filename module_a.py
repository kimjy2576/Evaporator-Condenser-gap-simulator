"""module_a.py — 모듈 A: Gap 열유동 모델 (재순환/혼합/복사/전도 → Cap Retention)"""

import numpy as np
from dataclasses import dataclass
from common import *

@dataclass
class GapParams:
    """
    Gap 열유동 파라미터

    gap_mode:
      'open'   — 개방 구조 (기존 v1 모델)
                  외부 재순환 경로 존재: 응축기 출구 → leak → 증발기 입구
                  건조기 린트필터 개방부, 벽걸이 에어컨 등
      'sealed' — 밀폐 덕트 (v2 신규)
                  외부 재순환 없음, Gap 내부 혼합만 존재
                  ① 급확대 와류 (recirculation eddy)
                  ② 열적 부력 (Richardson 수 기반)
                  ③ 제트 간 혼합
      'semi'   — 반밀폐 (건조기 실제 조건)
                  외부 재순환 (감쇠) + 내부 혼합 동시 적용
                  seal_fraction: 밀폐 비율 (0=완전 개방, 1=완전 밀폐)
    """
    T_amb:          float = 25.0
    RH_in:          float = 0.50    # 입구 상대습도 (건조기: 0.3~0.9)
    CMM:            float = 6.75    # 풍량 [m³/min]
    V_face:         float = 0.0     # 자동 계산됨
    gap_mode:       str   = "open"  # 'open' | 'sealed' | 'semi'
    seal_fraction:  float = 0.7     # 반밀폐 시 밀폐 비율 (0~1)
    mode:           str   = "forced"   # "forced" | "natural"
    baffle:         bool  = False
    frame_material: str   = "Al"
    A_frame:        float = 4e-4

    def set_V_face(self, A_face: float):
        """정면적으로부터 V_face 환산"""
        self.V_face = self.CMM / (60.0 * max(A_face, 1e-6))


def _recirculation_dT_open(gap_mm: float, p: GapParams,
                           T_cond_surf: float, T_amb: float) -> float:
    """
    [개방 구조] 외부 재순환 온도 상승량 (기존 v1 모델)

    응축기 출구 고온 공기 → 외부 leak path → 증발기 입구 재유입
    ΔT_recir = α · exp(−β · G)
    """
    gap = gap_mm / 1000.0
    dT  = T_cond_surf - T_amb
    if p.mode == "natural":
        alpha, beta = dT * 0.65, 18.0
    else:
        vf    = max(0.3, 1.0 - (p.V_face - 0.5) / 6.0)
        alpha = dT * 0.55 * vf
        beta  = 25.0 + p.V_face * 4.0
    return alpha * np.exp(-beta * gap) * (0.05 if p.baffle else 1.0)


def _mixing_effectiveness(gap_mm: float, p: GapParams, evap_spec,
                          T_cond_surf: float, T_evap_surf: float,
                          cond_spec=None) -> dict:
    """
    [밀폐 덕트] Gap 내부 혼합 효율 모델

    4가지 메커니즘의 합산:

    [메커니즘 ①: 급확대 재순환 와류] — 증발기 출구 → Gap
    증발기 코어 출구에서 급확대 시 유동 박리 → 코너 재순환 와류.
    σ_evap 작을수록 확대비 크고 와류 강함.

    [메커니즘 ②: 열적 부력 — Richardson 수]
    Gap 양벽 온도차(ΔT ≈ 40K)에 의한 자연대류 셀.

    [메커니즘 ③: 제트 간 와류]
    핀 채널별 제트 확산 + 상호 간섭.

    [메커니즘 ④: 급축소 상류 와류] — Gap → 응축기 입구 (v2 신규)
    Gap 공기가 응축기 코어로 진입할 때 급축소 발생.
    축소 직전(상류) 코너에 정체 영역 + 재순환 와류 형성.
    σ_cond 작을수록 축소비 크고 상류 와류 강함.

    물리:
      - 확대 와류(①)는 증발기 면에서 발생 → 냉기 쪽 혼합
      - 축소 와류(④)는 응축기 면에서 발생 → 고온 쪽 혼합
      - 두 와류 영역이 Gap 양끝에서 서로를 향해 확장
      - Gap이 작으면 두 영역이 중첩 → 혼합 급증

    η_cont ∝ (1 - σ_cond) — σ_cond 작을수록 축소 심화
      FT (σ≈0.55): 강한 축소 → η_cont 큼
      MCHX (σ≈0.76): 완만한 축소 → η_cont 작음

    Returns
    -------
    dict : eta_mix, eta_exp, eta_buoy, eta_jet, eta_cont, Ri, L_recirc
    """
    gap = max(gap_mm / 1000.0, 0.001)
    sigma_evap = 0.55
    if hasattr(evap_spec, 'fin_pitch') and hasattr(evap_spec, 'tube_do'):
        sigma_evap = (evap_spec.fin_pitch - evap_spec.fin_thickness) / evap_spec.fin_pitch * \
                (evap_spec.fin_height - evap_spec.tube_do) / evap_spec.fin_height
    elif hasattr(evap_spec, 'fin_pitch'):
        # MCHX: sigma from geo
        sigma_evap = getattr(evap_spec, 'sigma', 0.55)

    H_gap = evap_spec.H
    V_gap = p.V_face

    # ── ① 급확대 재순환 와류 (증발기 출구 → Gap) ────────────
    if hasattr(evap_spec, 'fin_pitch'):
        ch_h_evap = evap_spec.fin_pitch - evap_spec.fin_thickness
    else:
        ch_h_evap = 1.0e-3
    step_h_exp   = ch_h_evap * (1.0 - sigma_evap) * 0.5
    D_h_step_exp = max(2.0 * step_h_exp, 0.5e-3)
    L_recirc_exp = 6.0 * D_h_step_exp

    if gap < L_recirc_exp:
        eta_exp = 0.08 * (1.0 - gap / L_recirc_exp)**0.5
    else:
        eta_exp = 0.01 * np.exp(-(gap - L_recirc_exp) / (3.0 * L_recirc_exp))

    # ── ② 열적 부력 (Richardson 수) ──────────────────────────
    dT = abs(T_cond_surf - T_evap_surf)
    T_mean_K = (T_cond_surf + T_evap_surf) / 2.0 + 273.15
    beta_th = 1.0 / T_mean_K

    V_eff = max(V_gap, 0.1)
    Ri = 9.81 * beta_th * dT * gap / V_eff**2

    if Ri < 0.1:
        eta_buoy = 0.005 * Ri / 0.1
    elif Ri < 1.0:
        eta_buoy = 0.005 + 0.065 * (Ri - 0.1) / 0.9
    else:
        eta_buoy = min(0.15, 0.07 + 0.03 * np.log(Ri))

    if V_eff < 0.5:
        eta_buoy *= (1.0 + 2.0 * (0.5 - V_eff))

    # ── ③ 제트 간 와류 ───────────────────────────────────────
    # [v3 보정] Gap, V_face, 핀 피치 의존성 추가
    # 기존: η_jet = 0.02 고정 → 결과의 85% 지배 (비물리적)
    #
    # 물리:
    #   핀 채널에서 나오는 개별 제트가 Gap 공간에서 확산하며
    #   인접 제트와 간섭 → 소규모 와류 혼합
    #
    # 의존 변수:
    #   Fp (핀 피치): 작을수록 제트 간격 ↓ → 간섭 ↑ → η ↑
    #   V_face: 높을수록 제트 운동량 ↑ → 확산 빠름 → 혼합 ↑ (but also Gap 통과 빠름)
    #   Gap: 작을수록 제트가 응축기 면에 충돌 전 확산 미완 → 혼합 ↓
    #         클수록 제트가 완전 확산 → 혼합 ↑ → 어느 시점 이후 포화
    #
    # η_jet = η_jet_base × f(Gap) × f(V) × f(Fp)
    #   η_jet_base = 0.005 (제트 확산 기본값)
    #   f(Gap) = 1 - exp(-gap/L_mix), L_mix = 10mm (제트 확산 특성 길이)
    #   f(V) = (V/V_ref)^0.3, V_ref = 1.5 m/s
    #   f(Fp) = (Fp_ref/Fp)^0.5, Fp_ref = 2mm (피치 작을수록 간섭 ↑)

    Fp = evap_spec.fin_pitch if hasattr(evap_spec, 'fin_pitch') else 1.8e-3
    Fp_ref = 2.0e-3
    V_ref  = 1.5
    L_mix  = 0.010  # 제트 확산 특성 길이 [m]

    f_gap = 1.0 - np.exp(-gap / L_mix)
    f_vel = (max(V_gap, 0.3) / V_ref)**0.3
    f_fp  = (Fp_ref / max(Fp, 0.5e-3))**0.5

    eta_jet = 0.005 * f_gap * f_vel * f_fp
    eta_jet = np.clip(eta_jet, 0.001, 0.03)

    # ── ④ 급축소 상류 와류 (Gap → 응축기 입구) ───────────────
    eta_cont = 0.0
    L_recirc_cont = 0.0
    sigma_cond = 0.55  # 기본값

    if cond_spec is not None:
        # 응축기 σ_cond 계산
        if hasattr(cond_spec, 'fin_pitch') and hasattr(cond_spec, 'fin_height'):
            # FT 응축기
            sigma_cond = ((cond_spec.fin_pitch - cond_spec.fin_thickness)
                         / cond_spec.fin_pitch
                         * (cond_spec.fin_height - cond_spec.tube_do)
                         / cond_spec.fin_height)
            ch_h_cond = cond_spec.fin_pitch - cond_spec.fin_thickness
        elif hasattr(cond_spec, 'slab_pitch'):
            # MCHX 응축기
            sigma_cond = ((cond_spec.fin_pitch - cond_spec.fin_thickness)
                         / cond_spec.fin_pitch
                         * (cond_spec.slab_pitch - cond_spec.ch_height)
                         / cond_spec.slab_pitch)
            ch_h_cond = cond_spec.fin_pitch - cond_spec.fin_thickness
        else:
            ch_h_cond = 1.0e-3

        # 축소 step 높이 & 상류 재순환 길이
        # 축소 와류는 확대 와류보다 약함 (L_cont ≈ 3×D_h vs L_exp ≈ 6×D_h)
        step_h_cont   = ch_h_cond * (1.0 - sigma_cond) * 0.5
        D_h_step_cont = max(2.0 * step_h_cont, 0.5e-3)
        L_recirc_cont = 3.0 * D_h_step_cont   # 축소 재순환: ~확대의 절반

        if gap < L_recirc_cont:
            # 축소 와류가 Gap 상당 부분 침범
            eta_cont = 0.06 * (1.0 - sigma_cond) * (1.0 - gap / L_recirc_cont)**0.5
        else:
            # 축소 와류가 응축기 면 근처에만 국한
            eta_cont = 0.005 * (1.0 - sigma_cond) * np.exp(
                -(gap - L_recirc_cont) / (3.0 * L_recirc_cont))

        # ── 확대+축소 와류 중첩 효과 ────────────────────────
        # Gap이 두 재순환 영역의 합보다 작으면 와류가 중첩 → 혼합 증폭
        L_total = L_recirc_exp + L_recirc_cont
        if gap < L_total:
            overlap = 1.0 - gap / L_total
            # 중첩 시 20% 추가 혼합 (비선형 상호작용)
            eta_cont += 0.02 * overlap * (1.0 - sigma_cond)

    # ── 총 혼합 효율 ─────────────────────────────────────────
    eta_mix = min(eta_exp + eta_buoy + eta_jet + eta_cont, 0.30)

    return dict(
        eta_mix=eta_mix,
        eta_exp=eta_exp,
        eta_buoy=eta_buoy,
        eta_jet=eta_jet,
        eta_cont=eta_cont,
        sigma_cond=sigma_cond,
        Ri=Ri,
        L_recirc=L_recirc_exp * 1000,       # mm
        L_recirc_cont=L_recirc_cont * 1000,  # mm
    )


def recirculation_dT(gap_mm: float, p: GapParams,
                     T_cond_surf: float, T_amb: float,
                     evap_spec=None, cond_spec=None,
                     T_evap_surf: float = 5.0) -> dict:
    """
    Gap 모드에 따른 열적 페널티 통합 계산

    Returns
    -------
    dict : dT_recir, eta_mix, Q_mix_penalty, mix_detail, T_in_eff
    """
    rho_air = 1.18; cp_air = 1006.0

    gm = p.gap_mode.lower()

    if gm == 'open':
        dTr = _recirculation_dT_open(gap_mm, p, T_cond_surf, T_amb)
        return dict(dT_recir=dTr, eta_mix=0.0, Q_mix_penalty=0.0,
                    mix_detail=None, T_in_eff=T_amb + dTr)

    elif gm == 'sealed':
        dTr = 0.0
        if evap_spec is not None:
            mx = _mixing_effectiveness(gap_mm, p, evap_spec,
                                       T_cond_surf, T_evap_surf,
                                       cond_spec=cond_spec)
            eta_mix = mx['eta_mix']
        else:
            eta_mix = 0.03
            mx = None
        Q_mix = 0.0
        return dict(dT_recir=0.0, eta_mix=eta_mix, Q_mix_penalty=Q_mix,
                    mix_detail=mx, T_in_eff=T_amb)

    else:  # 'semi'
        sf = np.clip(p.seal_fraction, 0.0, 1.0)
        dTr = _recirculation_dT_open(gap_mm, p, T_cond_surf, T_amb) * (1.0 - sf)

        if evap_spec is not None:
            mx = _mixing_effectiveness(gap_mm, p, evap_spec,
                                       T_cond_surf, T_evap_surf,
                                       cond_spec=cond_spec)
            eta_mix = mx['eta_mix'] * sf
        else:
            eta_mix = 0.03 * sf
            mx = None

        Q_mix = 0.0
        return dict(dT_recir=dTr, eta_mix=eta_mix, Q_mix_penalty=Q_mix,
                    mix_detail=mx, T_in_eff=T_amb + dTr)


def view_factor(gap_mm: float, W: float) -> float:
    return 1.0 / (1.0 + 2.0 * gap_mm / 1000.0 / W)


def radiation_flux(gap_mm: float, T_cond: float, T_evap: float, W: float) -> float:
    Tc, Te = T_cond + 273.15, T_evap + 273.15
    eps_eff = 1.0 / (2.0 / EPS_AL - 1.0)
    return max(0.0, SIGMA_SB * eps_eff * view_factor(gap_mm, W) * (Tc**4 - Te**4))


def conduction_shortcircuit(gap_mm: float, p: GapParams,
                             T_cond: float, T_evap: float) -> float:
    """
    프레임 전도 단락 열전달 [W]

    [v3 보정] 3가지 물리적 제한 추가:

    ① 접촉저항 (R_contact)
       프레임-HX 코어 접합부에서 불완전 접촉에 의한 열저항.
       일반 금속-금속 접촉: R_c ≈ 1e-4 ~ 1e-3 m²K/W
       Al 프레임-Cu 튜브 (볼트/브레이징): R_c ≈ 5e-4 m²K/W 가정

    ② 프레임 핀효율 (η_frame)
       프레임은 HX 코어에서 Gap을 가로질러 뻗어있으며,
       측면이 Gap 내 공기에 노출되어 열을 빼앗김.
       → 유효 ΔT가 프레임 중앙에서 감소

       m = √(h_conv × P / (k × A_c))
       P: 프레임 둘레, A_c: 프레임 단면적
       η_frame = tanh(m×L) / (m×L)

    ③ 물리적 상한
       전도 열전달이 증발기 냉방능력을 초과할 수 없음.
       Q_cond ≤ 0.3 × Q_evap_max (보수적 상한)
    """
    k = FRAME_K[p.frame_material]
    L = max(gap_mm / 1000.0, 0.002)
    A_frame = p.A_frame

    # ① 접촉저항 (양쪽 접합부)
    R_contact = 5e-4   # [m²K/W] Al-Cu 접합 대표값
    R_cond    = L / k   # 프레임 전도 저항 [m²K/W]
    R_total   = R_cond + 2.0 * R_contact  # 양쪽 접촉저항

    # ② 프레임 핀효율
    # 프레임 단면: ~20mm × 2mm (대표값)
    w_frame = 0.020   # 프레임 폭 [m]
    t_frame = 0.002   # 프레임 두께 [m]
    P_frame = 2.0 * (w_frame + t_frame)  # 둘레
    h_conv  = 10.0    # Gap 내 자연대류+약한 강제대류 [W/m²K]
    m_frame = np.sqrt(h_conv * P_frame / (k * A_frame + 1e-9))
    mL = m_frame * L
    eta_frame = np.tanh(mL) / (mL + 1e-9) if mL > 0.01 else 1.0

    # 전도 열전달 (보정 적용)
    dT = T_cond - T_evap
    Q_cond = eta_frame * A_frame * dT / R_total

    # ③ 물리적 상한: Q_cond ≤ 100W (증발기 능력의 ~30%)
    Q_cond = min(Q_cond, 100.0)

    return max(Q_cond, 0.0)


def simulate_gap(gap_mm, evap_spec, cond_spec, evap_geo, cond_geo,
                 ua_evap, ua_cond, ref: RefrigerantState, gp: GapParams,
                 m_ref=0.00458, x_in=0.22, evap_corr='auto',
                 flow='counter', N_seg=5, _coil_ref_cache=None) -> dict:
    """
    단일 간격에서 전체 열유동 계산 (모듈 A) — Level 2 통합

    Level 2 Tube-Segment 모델 사용:
    - 세그먼트별 건도 추적 + T_wall 반복 수렴
    - 상관식 자동 선택 (evap_corr='auto')
    """
    from common import compute_coil_v3

    T_evap_surf = ref.T_sat_evap
    T_cond_surf = ref.T_sat_cond
    W_eff = max(evap_spec.W, cond_spec.W)
    A_face = evap_spec.W * evap_spec.H
    gp.set_V_face(A_face)

    # 재순환 + 혼합 통합 계산
    gap_result = recirculation_dT(gap_mm, gp, T_cond_surf, gp.T_amb,
                                  evap_spec=evap_spec, cond_spec=cond_spec,
                                  T_evap_surf=T_evap_surf)
    dTr       = gap_result['dT_recir']
    eta_mix   = gap_result['eta_mix']
    Q_mix     = gap_result['Q_mix_penalty']
    T_in_eff  = gap_result['T_in_eff']

    RH_in = getattr(gp, 'RH_in', 0.50)
    P_atm = 101325.0
    rho_air = P_atm / (287.05 * (T_in_eff + 273.15))
    cp_air = 1006.0
    m_dot_air = rho_air * gp.V_face * A_face

    # ── Level 2 코일 모델 ──
    coil = compute_coil_v3(evap_spec, evap_geo, ref,
                           T_in_eff, RH_in, gp.V_face,
                           m_ref, x_in, 'evap', N_seg, flow, evap_corr)
    if _coil_ref_cache is not None:
        coil_ref = _coil_ref_cache
    else:
        coil_ref = compute_coil_v3(evap_spec, evap_geo, ref,
                                   gp.T_amb, RH_in, gp.V_face,
                                   m_ref, x_in, 'evap', N_seg, flow, evap_corr)

    Q_total    = coil['Q_total']
    Q_sensible = coil['Q_sen']
    Q_latent   = coil['Q_lat']
    SHR        = coil['SHR']
    T_evap_out = coil['T_out']
    W_in       = coil['W_in']
    W_out      = coil.get('W_out', W_in)
    T_dp       = coil['T_dp']
    h_in       = coil['h_in']
    h_out      = coil.get('h_out', h_in - Q_total / (m_dot_air + 1e-9))
    h_amb      = coil_ref['h_in']
    Q_evap     = Q_total

    eps = coil.get('eps', Q_total / (m_dot_air * cp_air * max(abs(T_in_eff - T_evap_surf), 0.1) + 1e-9))
    NTU = -np.log(1 - min(eps, 0.999)) if eps < 0.999 else 7.0

    Q_useful = max(m_dot_air * (h_amb - h_out), 0.0)
    Q_recir_waste = max(Q_total - Q_useful, 0.0)
    Q_ref = max(coil_ref['Q_total'], 1.0)
    dehumid_rate = coil.get('dehumid_rate', 0)

    q_rad  = radiation_flux(gap_mm, T_cond_surf, T_evap_surf, W_eff)
    Q_rad  = q_rad * evap_spec.W * evap_spec.H
    Q_cond_loss = conduction_shortcircuit(gap_mm, gp, T_cond_surf, T_evap_surf)

    Q_mix = eta_mix * Q_useful
    Q_net = max(Q_useful - Q_rad - Q_cond_loss - Q_mix, 0.0)
    cap_ratio = min(100.0, Q_net / Q_ref * 100.0)

    return dict(gap_mm=gap_mm, dT_recir=dTr, T_in_eff=T_in_eff,
                q_rad=q_rad, Q_rad=Q_rad, Q_cond=Q_cond_loss,
                Q_evap=Q_evap, Q_total=Q_total,
                Q_sensible=Q_sensible, Q_latent=Q_latent, SHR=SHR,
                Q_useful=Q_useful, Q_recir_waste=Q_recir_waste,
                Q_mix=Q_mix, Q_net=Q_net,
                cap_ratio=cap_ratio,
                dehumid_rate=dehumid_rate,
                T_evap_out=T_evap_out, T_dp=T_dp,
                W_in=W_in, W_out=W_out,
                UA_evap=ua_evap['UA'], NTU=NTU, eps=eps,
                V_face=gp.V_face, CMM=gp.CMM,
                eta_mix=eta_mix, gap_mode=gp.gap_mode,
                mix_detail=gap_result.get('mix_detail'),
                wet_rows=0, Nr=coil.get('Nr', 1),
                q_cond=coil.get('q_cond', 0),
                phase_summary=coil.get('phase_summary', {}),
                x_out=coil.get('x_out', 0),
                T_ref_out=coil.get('T_ref_out', 0))


def sweep(gaps, evap_spec, cond_spec, evap_geo, cond_geo,
          ua_evap, ua_cond, ref, gp,
          m_ref=0.00458, x_in=0.22, evap_corr='auto',
          flow='counter', N_seg=5) -> list:
    # ★ 기준 조건 (Gap→∞) 1회만 계산 → 캐싱
    from common import compute_coil_v3
    RH_in = getattr(gp, 'RH_in', 0.50)
    A_face = evap_spec.W * evap_spec.H
    V_face = gp.CMM / (60.0 * A_face)
    coil_ref = compute_coil_v3(evap_spec, evap_geo, ref,
                               gp.T_amb, RH_in, V_face,
                               m_ref, x_in, 'evap', N_seg, flow, evap_corr)
    return [simulate_gap(g, evap_spec, cond_spec, evap_geo, cond_geo,
                         ua_evap, ua_cond, ref, gp,
                         m_ref, x_in, evap_corr, flow, N_seg,
                         _coil_ref_cache=coil_ref) for g in gaps]

