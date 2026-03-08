"""module_b.py — 모듈 B: 비말동반 위험도 (습면열전달 → Weber → MC → Flux)"""

import numpy as np
from common import *
from module_a import GapParams, simulate_gap

def _psat(T):
    Tk = T + 273.15
    return np.exp(-5.8002206e3/Tk + 1.3914993 - 4.864e-2*Tk
                  + 4.1765e-5*Tk**2 - 1.4452e-8*Tk**3 + 6.5460*np.log(Tk))
def _hum(T, RH, P=P_ATM):
    Pw = RH * _psat(T); return 0.62198 * Pw / (P - Pw)
def _enth(T, W): return CP_AIR*T + W*(2501000 + CP_VAP*T)
def _dp(T, RH):
    a, b = 17.625, 243.04
    g = np.log(max(RH, 1e-9)) + a*T/(b+T); return b*g/(a-g)
def _dWsdT(T):
    dT = 0.05; return (_hum(T+dT, 1.0) - _hum(T-dT, 1.0)) / (2*dT)


def compute_condensate(spec, T_in, RH_in, T_wall, V_face):
    """
    Threlkeld ε-NTU 습면 열전달
    T_wall ≥ T_dp 이면 건면 → q_cond = 0

    [v2 변경] tube_layout에 따라 j-factor 분기
      Staggered: 기존 상관식 (높은 j → 높은 h_o → 많은 응결수)
      Inline:    j_inline ≈ 0.78 × j_stag (wake 차폐 → h_o 감소 → 응결수 감소)
    """
    T_dp = _dp(T_in, RH_in)
    W_in = _hum(T_in, RH_in)
    h_in = _enth(T_in, W_in)

    rho = P_ATM*(1+W_in)/(287*(1+W_in*461.5/287)*(T_in+273.15))
    mu  = MU_AIR*((T_in+273.15)/298.15)**0.7
    k_a = K_AIR*((T_in+273.15)/298.15)**0.82
    cp_m = CP_AIR + W_in*CP_VAP; Pr = cp_m*mu/k_a
    G_c  = rho*V_face/spec.sigma; Re = G_c*spec.Do/mu
    Nr, Fp, Do, Pl, Pt = spec.Nr, spec.Fp, spec.Do, spec.Pl, spec.Pt

    # j-factor: fin_type(5종) × tube_layout(2종) 분기
    # [v4] common.py 직접 상관식에서 산출된 E_j 비율 적용
    layout = getattr(spec, 'tube_layout', 'staggered').lower()
    fin_type = getattr(spec, 'fin_type', 'plain').lower()

    j_stag = max(0.003, 0.086*Re**(-0.361-0.042*Nr/np.log(Re+1e-6)
            +0.158*np.log(Nr*(Fp/Do)**0.41))
            *Nr**(-0.170)*(Fp/Do)**(-1.224)*(Fp/Pl)**(-0.508)*(Pt/Pl)**0.257)

    # common.py 직접 상관식 기반 강화 비율 계산
    from common import (_j_plain_staggered, _j_wavy_staggered,
                         _j_louvered_staggered, _j_slit_staggered,
                         FinTubeSpec)
    _tmp = FinTubeSpec(
        fin_pitch=Fp, fin_thickness=getattr(spec, 'df', 0.1e-3),
        fin_height=Pt, tube_do=Do, tube_rows=Nr,
        tube_pitch_t=Pt, tube_pitch_l=Pl, fin_type=fin_type,
        wavy_height=getattr(spec, 'wavy_height', 1.5e-3),
        wavy_angle=getattr(spec, 'wavy_angle', 17.5),
        ft_louver_pitch=getattr(spec, 'ft_louver_pitch', 1.7e-3),
        ft_louver_angle=getattr(spec, 'ft_louver_angle', 28.0),
        slit_num=getattr(spec, 'slit_num', 6),
        slit_height=getattr(spec, 'slit_height', 1.0e-3),
    )
    j_plain_ref = _j_plain_staggered(Re, _tmp)
    _fn = {'wavy': _j_wavy_staggered, 'louvered': _j_louvered_staggered,
           'louver': _j_louvered_staggered, 'slit': _j_slit_staggered,
           }.get(fin_type)
    if _fn:
        E_j = _fn(Re, _tmp) / (j_plain_ref + 1e-9)
        j_stag *= max(E_j, 1.0)

    if layout in ('inline', 'parallel') and Nr >= 2:
        F_layout = 0.82 if Nr == 2 else 0.75
        j = j_stag * F_layout
    else:
        j = j_stag

    h_o = j * G_c * cp_m / Pr**(2/3)

    if T_wall >= T_dp:
        UA = h_o * spec.A_total; m_dot = G_c * spec.Ac
        eps = 1 - np.exp(-UA/m_dot)
        Q_sen = eps * m_dot * CP_AIR * (T_in - T_wall)
        return dict(q_cond=0.0, Q_total=Q_sen, Q_lat=0.0, Q_sen=Q_sen,
                    SHR=1.0, T_dp=T_dp, T_wall=T_wall,
                    wet_fraction=0.0, m_dot=m_dot)

    h_fg = 2501000 - 2381*T_wall; T_fin = (T_in+T_wall)/2
    xi   = _dWsdT(T_fin)*h_fg/(LE**(2/3)*CP_AIR)
    m_f  = np.sqrt(2*h_o/(K_AL*spec.df))*np.sqrt(1+max(xi,0))
    r1   = Do/2; r2e = max(1.27*Pt*np.sqrt(Pl/Pt-0.3), r1*1.01)
    phi  = (r2e/r1-1)*(1+0.35*np.log(r2e/r1)); arg = m_f*r1*phi
    eta_f = np.tanh(arg)/(arg+1e-9)
    eta_o = 1 - spec.A_fin/spec.A_total*(1-eta_f)
    cp_s  = CP_AIR + h_fg*_dWsdT(T_wall)
    UA_s  = cp_s/CP_AIR*eta_o*h_o*spec.A_total
    m_dot = G_c*spec.Ac; NTU = UA_s/m_dot; eps = 1 - np.exp(-NTU)
    W_sw  = _hum(T_wall, 1.0); h_sw = _enth(T_wall, W_sw)
    Q     = max(eps*m_dot*(h_in-h_sw), 0)
    W_out = max(W_in - eps*max(W_in-W_sw, 0), W_sw)
    Q_lat = m_dot*(W_in-W_out)*h_fg; Q_sen = Q-Q_lat
    return dict(q_cond=m_dot*(W_in-W_out)*1000,
                Q_total=Q, Q_lat=Q_lat, Q_sen=Q_sen,
                SHR=Q_sen/(Q+1e-9), T_dp=T_dp,
                T_wall=T_wall, wet_fraction=1.0, m_dot=m_dot)


def we_ch(V_face, spec):   return RHO_A*(V_face/spec.sigma)**2*spec.Fp/SIGMA_W
def we_crit(spec):          return 0.019/np.sqrt(spec.Fp)
def v_onset(spec):
    Wc = we_crit(spec); return np.sqrt(Wc*SIGMA_W/(RHO_A*spec.Fp))*spec.sigma

def eta_co(V_face, spec, coeff):
    Wch = we_ch(V_face, spec); Wc = we_crit(spec)
    if Wch <= Wc: return 0.0
    return min(0.30, coeff*(Wch/Wc-1.0)**1.2)


def sample_bimodal(N, rng, w_bridge=0.75):
    Nb = int(N*w_bridge); Ns = N-Nb
    ds = rng.lognormal(np.log(D50_SMALL),  SIG_SMALL,  Ns)
    db = rng.lognormal(np.log(D50_BRIDGE), SIG_BRIDGE, Nb)
    d  = np.concatenate([np.clip(ds,5e-6,1e-3), np.clip(db,0.5e-3,8e-3)])
    return d[rng.permutation(N)]


def _Cd(Re):
    Re = np.maximum(Re, 1e-9)
    return np.where(Re<=1.0, 24/Re,
           np.where(Re<=1000, 24/Re*(1+0.15*Re**0.687), 0.44))


def track_batch(d_arr, y0_arr, gap, V_face, spec, theta_deg=0.0, dt=8e-5, t_max=0.5):
    """
    RK4 액적 궤적 일괄 적분
    공기 흐름 방향(x): 증발기 출구 → 응축기 입구 (= 주 흐름 방향)
    """
    N = len(d_arr); H = spec.H
    y_lo = -H/2+spec.h_drain; y_hi = H/2
    s = np.zeros((N,4)); s[:,1] = y0_arr
    alive   = np.ones(N, dtype=bool)
    reached = np.zeros(N, dtype=bool); y_hit = np.full(N, np.nan)
    m_arr = RHO_W*np.pi*d_arr**3/6; Ap_arr = np.pi*d_arr**2/4
    V_ch  = V_face / spec.sigma

    th = np.radians(theta_deg)
    gx =  G_GRAV*np.sin(th)
    gy = -G_GRAV*np.cos(th)

    def _accel(s_):
        x_,y_,vx_,vy_ = s_[:,0],s_[:,1],s_[:,2],s_[:,3]
        tt   = np.clip(x_/max(gap,1e-9), 0, 1)
        Vax  = V_ch*(1-0.30*tt); Vay = 0.03*V_ch*np.sin(np.pi*tt)
        dvx  = Vax-vx_; dvy = Vay-vy_
        Vrel = np.maximum(np.sqrt(dvx**2+dvy**2), 1e-12)
        Fd   = 0.5*RHO_A*Vrel**2*_Cd(RHO_A*Vrel*d_arr/MU_A)*Ap_arr
        ax   = Fd*dvx/(m_arr*Vrel) - gx
        ay   = Fd*dvy/(m_arr*Vrel) + gy
        return np.column_stack([vx_,vy_,ax,ay])

    t = 0.0
    while t < t_max and alive.any():
        k1 = _accel(s); k2 = _accel(s+dt/2*k1)
        k3 = _accel(s+dt/2*k2); k4 = _accel(s+dt*k3)
        s_new = s + dt/6*(k1+2*k2+2*k3+k4)
        x_n = s_new[:,0]; y_n = s_new[:,1]
        hit  = alive & (x_n >= gap)
        miss = alive & ((y_n<y_lo)|(y_n>y_hi)|(x_n<-0.005))
        reached[hit] = True; y_hit[hit] = y_n[hit]
        alive[hit|miss] = False; s[alive] = s_new[alive]; t += dt

    return reached, y_hit


def monte_carlo(gap_mm, V_face, spec, theta_deg=0.0, w_bridge=0.75, N=300, seed=42):
    rng = np.random.RandomState(seed)
    gap = gap_mm / 1000
    d_arr  = sample_bimodal(N, rng, w_bridge)
    y0_arr = rng.uniform(-spec.H/2+spec.h_drain, spec.H/2, N)
    reached, y_hit = track_batch(d_arr, y0_arr, gap, V_face, spec, theta_deg)
    is_small = (d_arr < 0.5e-3)
    return dict(
        P_reach =reached.mean(),
        P_small =reached[is_small].mean()  if is_small.sum()>0  else 0.0,
        P_bridge=reached[~is_small].mean() if (~is_small).sum()>0 else 0.0,
        d_arr=d_arr, y0_arr=y0_arr, reached=reached, is_small=is_small
    )



# ═══════════════════════════════════════════════════════════════════
#  비말동반 위험도 모델 — 모듈 B MCHX 확장
#  ① Chang & Wang (1997) IJHMT 40(3):533-544  — 건면 j-factor
#  ② McLaughlin & Webb (2000) SAE 2000-01-0574 — 습면 브리징 보정
# ═══════════════════════════════════════════════════════════════════

def _j_chang_wang_1997(Re_Lp, Lp, Fp, Fh, Fd, Ll, tp, delta_f, theta_deg):
    """
    Chang & Wang (1997) MCHX 루버 핀 건면 j-factor

    j = Re_Lp^J1 × (θ/90)^0.27 × (Fp/Lp)^-0.14 × (Fh/Lp)^-0.29
        × (Fd/Lp)^-0.23 × (Ll/Lp)^0.66 × (tp/Lp)^-0.02 × (δf/Lp)^-0.48
    J1 = -0.49 × (θ/90)^0.27

    적용범위: 100 ≤ Re_Lp ≤ 3000, 91개 샘플, 평균편차 5.94%
    """
    th   = np.clip(theta_deg, 15.0, 45.0)
    tn   = th / 90.0
    J1   = -0.49 * tn**0.27
    j    = (Re_Lp**J1 * tn**0.27
            * (Fp/Lp)**(-0.14) * (Fh/Lp)**(-0.29)
            * (Fd/Lp)**(-0.23) * (Ll/Lp)**0.66
            * (tp/Lp)**(-0.02) * (delta_f/Lp)**(-0.48))
    return max(j, 0.003)


def _wet_penalty_mclaughlin_webb(Lp, theta_deg, is_hydrophilic=False):
    """
    McLaughlin & Webb (2000) 브리징 기반 습면 h_o 보정계수

    Lp_crit(θ) = LP_CRIT_MCHX × sin(30°)/sin(θ)
    Lp < Lp_crit → 심화 브리징: F_wet ≈ 0.50
    Lp ≥ Lp_crit → 경미한 습면: F_wet ≈ 0.85
    친수 코팅:  ×1.25 (McLaughlin & Webb: h_o,wet +25%)
    """
    theta_corr = np.sin(np.radians(30.0)) / np.sin(np.radians(max(theta_deg, 10.0)))
    Lp_crit    = LP_CRIT_MCHX * theta_corr

    if Lp < Lp_crit:
        frac  = Lp / Lp_crit
        F_wet = F_WET_MILD * frac + F_WET_BRIDGING * (1.0 - frac)
    else:
        F_wet = F_WET_MILD

    if is_hydrophilic:
        F_wet *= F_WET_HYDRO
    return min(F_wet, 1.0)


def compute_condensate_mchx(spec: 'MCHXCarryoverSpec',
                              T_in, RH_in, T_wall, V_face,
                              is_hydrophilic=False):
    """
    MCHX 습면 열전달 — Threlkeld ε-NTU

    [FT compute_condensate와의 주요 차이]
      j-factor  : Chang & Wang (1997) 루버 핀 상관식 (Re_Lp 기반)
      핀 효율   : 직사각형 핀 근사 — η_f = tanh(m·Fh/2)/(m·Fh/2)
      습면 보정  : McLaughlin & Webb (2000) 브리징 패널티 (F_wet)
    """
    T_dp = _dp(T_in, RH_in)
    W_in = _hum(T_in, RH_in)
    h_in = _enth(T_in, W_in)

    rho  = P_ATM*(1+W_in)/(287*(1+W_in*461.5/287)*(T_in+273.15))
    mu   = MU_AIR*((T_in+273.15)/298.15)**0.7
    k_a  = K_AIR*((T_in+273.15)/298.15)**0.82
    cp_m = CP_AIR + W_in*CP_VAP
    Pr   = cp_m*mu/k_a

    G_c   = rho*V_face/spec.sigma
    Re_Lp = np.clip(G_c*spec.Lp/mu, 100.0, 3000.0)

    # Chang & Wang (1997) 건면 j-factor
    j_dry   = _j_chang_wang_1997(Re_Lp, spec.Lp, spec.Fp, spec.Fh,
                                   spec.Fd, spec.Ll, spec.tp, spec.df,
                                   spec.theta_deg)
    h_o_dry = j_dry * G_c * cp_m / Pr**(2/3)

    # McLaughlin & Webb (2000) 습면 보정
    F_wet = (_wet_penalty_mclaughlin_webb(spec.Lp, spec.theta_deg, is_hydrophilic)
             if T_wall < T_dp else 1.0)
    h_o   = h_o_dry * F_wet

    # 직사각형 MCHX 핀 효율 계산 공통 함수
    def _rect_eta_o(h_o_val, extra_xi=0.0):
        m_f = np.sqrt(2*h_o_val*(1+max(extra_xi,0)) / (K_AL*spec.df+1e-15))
        arg = m_f * spec.Fh / 2
        ef  = np.tanh(arg)/(arg+1e-9)
        return 1 - spec.A_fin/spec.A_total*(1-ef)

    m_dot = G_c * spec.Ac

    # ── 건면 처리 ──────────────────────────────────────────────
    if T_wall >= T_dp:
        eta_o = _rect_eta_o(h_o)
        UA    = eta_o * h_o * spec.A_total
        eps   = 1 - np.exp(-UA/(m_dot*CP_AIR+1e-9))
        Q_sen = eps * m_dot * CP_AIR * (T_in - T_wall)
        return dict(q_cond=0.0, Q_total=Q_sen, Q_lat=0.0, Q_sen=Q_sen,
                    SHR=1.0, T_dp=T_dp, T_wall=T_wall,
                    wet_fraction=0.0, m_dot=m_dot,
                    j_dry=j_dry, F_wet=F_wet, Re_Lp=Re_Lp)

    # ── 습면 처리 (Threlkeld ε-NTU) ────────────────────────────
    h_fg = 2501000 - 2381*T_wall
    T_fin = (T_in + T_wall) / 2
    xi    = _dWsdT(T_fin) * h_fg / (LE**(2/3)*CP_AIR)
    eta_o = _rect_eta_o(h_o, xi)

    cp_s  = CP_AIR + h_fg*_dWsdT(T_wall)
    UA_s  = cp_s/CP_AIR * eta_o * h_o * spec.A_total
    NTU   = UA_s/(m_dot+1e-9)
    eps   = 1 - np.exp(-NTU)

    W_sw  = _hum(T_wall, 1.0)
    h_sw  = _enth(T_wall, W_sw)
    Q     = max(eps*m_dot*(h_in-h_sw), 0)
    W_out = max(W_in - eps*max(W_in-W_sw, 0), W_sw)
    Q_lat = m_dot*(W_in-W_out)*h_fg
    Q_sen = Q - Q_lat

    return dict(q_cond=m_dot*(W_in-W_out)*1000,
                Q_total=Q, Q_lat=Q_lat, Q_sen=Q_sen,
                SHR=Q_sen/(Q+1e-9), T_dp=T_dp,
                T_wall=T_wall, wet_fraction=1.0, m_dot=m_dot,
                j_dry=j_dry, F_wet=F_wet, Re_Lp=Re_Lp)


# ── MCHX Weber 수 기반 이탈 조건 (Lp 특성 길이) ──────────────────

def we_ch_mchx(V_face, spec):
    """MCHX 채널 Weber 수 — 특성 길이 Lp (루버 피치)"""
    return RHO_A*(V_face/spec.sigma)**2 * spec.Lp / SIGMA_W

def we_crit_mchx(spec):
    """MCHX 임계 Weber 수 (FT의 Fp → Lp 치환)"""
    return 0.019 / np.sqrt(spec.Lp)

def v_onset_mchx(spec):
    """MCHX 비말 이탈 개시 면속도 [m/s]"""
    Wc = we_crit_mchx(spec)
    return np.sqrt(Wc*SIGMA_W/(RHO_A*spec.Lp))*spec.sigma

def eta_co_mchx(V_face, spec, coeff):
    """
    MCHX 비말동반 효율 η_co

    [물리적 메커니즘 — McLaughlin & Webb 2000]
    ① 공기역학적 이탈 (aerodynamic): We_Lp > We_c → FT와 동일 메커니즘
       특성길이 Lp (루버 피치) 사용
    ② 브리징 구동 배출 (bridging-driven): Lp < Lp_crit 시 루버 채널 폐색 →
       채널 출구 강제 액체 배출 → V_onset 미만에서도 carryover 발생
       (FT에는 없는 MCHX 고유 메커니즘)

    친수 코팅: MCHX에서는 브리징 악화 → eta_co 증가 (FT와 반대 거동)
    단, eta_coeff가 낮으면 compute_condensate_mchx에서 is_hydrophilic=True
    → F_wet 증가 → q_cond 증가 → flux 최종 증가
    """
    th_corr  = np.sin(np.radians(30.0)) / np.sin(np.radians(max(spec.theta_deg, 10.0)))
    Lp_crit  = LP_CRIT_MCHX * th_corr

    # ② 브리징 구동 배출 (Lp < Lp_crit 에서만 존재)
    bridge_base = 0.0
    if spec.Lp < Lp_crit:
        # 브리징 강도: Lp가 Lp_crit에서 멀수록 강함
        bridge_intensity = 1.0 - spec.Lp / Lp_crit   # 0~1
        # 속도 의존성: 낮은 V에서도 존재하나 V 비례로 강화
        v_factor = min(1.0, V_face / 2.0)
        bridge_base = 0.12 * bridge_intensity * v_factor * (coeff / 5e-4)**0.3

    # ① 공기역학적 이탈 (We_Lp 기반)
    Wch = we_ch_mchx(V_face, spec)
    Wc  = we_crit_mchx(spec)

    if Wch > Wc:
        eta_aero = min(0.35, coeff * (Wch / Wc - 1.0)**1.2)
        if spec.Lp < Lp_crit:
            # 브리징 가속 효과
            bridge_factor = 1.0 + 2.0 * (1.0 - spec.Lp / Lp_crit)
            eta_aero = min(0.45, eta_aero * bridge_factor)
    else:
        eta_aero = 0.0

    return min(0.45, bridge_base + eta_aero)



def risk_level(f):
    if f < FLUX_CAUTION: return "안전"
    if f < FLUX_DANGER:  return "주의"
    if f < FLUX_SEVERE:  return "위험"
    return "매우위험"

def risk_score(f): return min(1.0, f/(FLUX_SEVERE*2))


# ═══════════════════════════════════════════════════════════════════
#  통합 케이스 분석 — Gap sweep 시 모듈 A + B 동시 계산
# ═══════════════════════════════════════════════════════════════════

def analyze_combined(case: dict, gaps: np.ndarray,
                     evap_spec, geo_evap: dict,
                     cond_spec, geo_cond: dict,
                     ua_evap: dict, ua_cond: dict,
                     ref: RefrigerantState) -> dict:
    """
    단일 케이스에 대해 Gap sweep → 모듈 A(냉방성능) + 모듈 B(비말동반) 동시 계산

    [통합 핵심]
    - 모듈 A: 재순환/복사/전도 열침입 → Cap Retention [%]
    - 모듈 B: 습면 응결수 → Weber 이탈 → MC 궤적 → Flux [mg/m²s]
      · FT 증발기: FT compute_condensate + Fp-기반 Weber 수
      · MCHX 증발기: Chang&Wang(1997) j-factor + McLaughlin&Webb(2000) 브리징 보정
                      Lp-기반 Weber 수 + 브리징 패널티
    - sweet_gap: 두 기준(Cap≥95%, Flux≤FLUX_DANGER)을 동시 만족하는 최소 Gap
    """
    CMM = case['CMM']
    A_face_evap = evap_spec.W * evap_spec.H
    V_face = CMM / (60.0 * A_face_evap)
    gap_mode = case.get('gap_mode', 'semi')
    seal_frac = case.get('seal_fraction', 0.7)
    gp = GapParams(CMM=CMM, T_amb=25.0, mode="forced",
                   gap_mode=gap_mode, seal_fraction=seal_frac)
    gp.set_V_face(A_face_evap)
    is_hydrophilic = (case.get('eta_coeff', 5e-4) <= 2e-4)  # 낮은 η_coeff = 친수

    # ── 공통 열전달 엔진으로 q_cond 산출 ──────────────────────
    if hasattr(evap_spec, 'hx_type') and evap_spec.hx_type == 'MCHX':
        # MCHX: 기존 방식 유지 (공통 엔진 미지원)
        co_spec = MCHXCarryoverSpec(evap_spec, geo_evap)
        cr      = compute_condensate_mchx(co_spec, case['T_in'], case['RH_in'],
                                           case['T_wall'], V_face, is_hydrophilic)
        V_on    = v_onset_mchx(co_spec)
        _eta_fn = eta_co_mchx
    else:
        # FT: 공통 Tube-Row Segmented Threlkeld 엔진 사용
        from common import compute_coil_performance_segmented
        coil = compute_coil_performance_segmented(evap_spec, geo_evap, ua_evap,
                                                   case['T_in'], case['RH_in'],
                                                   case['T_wall'], V_face)
        # cr dict 호환 형식으로 매핑
        cr = dict(q_cond=coil['q_cond'], Q_total=coil['Q_total'],
                  Q_lat=coil['Q_lat'], Q_sen=coil['Q_sen'],
                  SHR=coil['SHR'], T_dp=coil['T_dp'],
                  T_wall=case['T_wall'],
                  wet_fraction=1.0 if coil['is_wet'] else 0.0,
                  m_dot=coil['m_dot'])

        co_spec = CarryoverSpec(evap_spec, geo_evap)
        V_on    = v_onset(co_spec)
        _eta_fn = eta_co

    cap_arr   = []
    flux_arr  = []   # Mid η_co 기준

    for g in gaps:
        # 모듈 A: Cap Retention
        r_gap = simulate_gap(g, evap_spec, cond_spec, geo_evap, geo_cond,
                              ua_evap, ua_cond, ref, gp)
        cap_arr.append(r_gap['cap_ratio'])

        # 모듈 B: Carryover Flux (Mid η)
        mc  = monte_carlo(g, V_face, co_spec, case['theta_deg'],
                          case['w_bridge'], N=40, seed=42)
        eta = _eta_fn(V_face, co_spec, case['eta_coeff'])
        A   = co_spec.W * co_spec.H
        flux_arr.append(cr['q_cond'] * eta * 1000 * mc['P_reach'] / A)

    cap_arr  = np.array(cap_arr)
    flux_arr = np.array(flux_arr)

    # sweet spot: Cap ≥ 95% AND Flux ≤ FLUX_DANGER
    sweet_gap = None
    for g, cap, fl in zip(gaps, cap_arr, flux_arr):
        if cap >= 95.0 and fl <= FLUX_DANGER:
            sweet_gap = g
            break

    return dict(
        case=case, cr=cr, V_on=V_on,
        gaps=gaps, cap_arr=cap_arr, flux_arr=flux_arr,
        sweet_gap=sweet_gap,
        evap_hx_type=getattr(evap_spec, 'hx_type', 'FT'),
    )


def compute_carry_penalty(result_b, evap_spec, evap_geo, gaps,
                          T_cond_surf=45.0):
    """
    비말동반 → 응축기 잠열 페널티 계산

    물리: 비말동반 물방울이 응축기 표면에서 재증발
          → 증발잠열 h_fg를 응축기 방열에서 빼앗김

    Parameters
    ----------
    result_b : analyze_combined 반환값
    evap_spec : 증발기 스펙
    evap_geo : 증발기 기하
    gaps : gap 배열 [mm]
    T_cond_surf : 응축기 표면 온도 [°C]

    Returns
    -------
    np.ndarray : Q_carry_penalty [W] (gap별)
    """
    h_fg = 2501000 - 2381 * T_cond_surf  # J/kg
    n = len(gaps)
    penalty = np.zeros(n)

    q_cond = result_b['cr']['q_cond']
    if q_cond <= 0:
        return penalty

    case = result_b['case']
    V_face = case['CMM'] / (60.0 * evap_spec.W * evap_spec.H)
    eta_coeff = case.get('eta_coeff', 5e-4)

    hx_type = getattr(evap_spec, 'hx_type', 'FT')
    if hx_type == 'MCHX':
        co_spec = MCHXCarryoverSpec(evap_spec, evap_geo)
        eta = eta_co_mchx(V_face, co_spec, eta_coeff)
    else:
        co_spec = CarryoverSpec(evap_spec, evap_geo)
        eta = eta_co(V_face, co_spec, eta_coeff)

    if eta <= 0:
        return penalty

    for i, g in enumerate(gaps):
        mc = monte_carlo(g, V_face, co_spec,
                         case.get('theta_deg', 0),
                         case.get('w_bridge', 0.75),
                         N=40, seed=42)
        m_carry = q_cond * 1e-3 * eta * mc['P_reach']  # kg/s
        penalty[i] = m_carry * h_fg  # W

    return penalty


def apply_carry_penalty(results_a, carry_penalty):
    """
    Module A 결과에 비말동반 응축기 페널티 반영

    results_a : list of simulate_gap 결과 dict
    carry_penalty : np.ndarray [W]

    추가 필드:
      Q_carry_penalty : 응축기 잠열 페널티 [W]
      cap_corrected : 보정된 Cap Retention [%]
    """
    for i, r in enumerate(results_a):
        qc = carry_penalty[i] if i < len(carry_penalty) else 0
        r['Q_carry_penalty'] = qc
        if r['cap_ratio'] > 0 and qc > 0:
            # Q_ref 역산 → 페널티 반영
            Q_ref = r['Q_net'] / (r['cap_ratio'] / 100 + 1e-9)
            r['cap_corrected'] = max(0, (r['Q_net'] - qc) / (Q_ref + 1e-9) * 100)
        else:
            r['cap_corrected'] = r['cap_ratio']

