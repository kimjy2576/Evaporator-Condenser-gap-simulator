"""common.py — 공용 모듈: 상수, HX 스펙, 냉매, 기하, HTC, UA"""

import numpy as np
from dataclasses import dataclass, field
import CoolProp.CoolProp as CP

# ─── 물리 상수 ─────────────────────────────────────────
SIGMA_SB = 5.67e-8   # Stefan-Boltzmann [W/m²K⁴]
EPS_AL   = 0.08      # 알루미늄 방사율
K_AL     = 160.0     # Al 열전도도 [W/mK]
K_CU     = 385.0     # Cu 열전도도 [W/mK]
FRAME_K  = {"Al": 160.0, "Steel": 50.0, "SUS": 15.0, "PP": 0.2}

# 비말동반 물리 상수
P_ATM   = 101325; CP_AIR = 1006; CP_VAP = 1860
K_AIR   = 0.0263; MU_AIR = 1.84e-5; LE = 0.845
RHO_W   = 998.0;  MU_A   = 1.84e-5; RHO_A = 1.20
SIGMA_W = 0.0728; G_GRAV = 9.81

D50_SMALL  = 80e-6;  SIG_SMALL  = 0.50
D50_BRIDGE = 3e-3;   SIG_BRIDGE = 0.30

# ─── 스타일 ───────────────────────────────────────────────────────
DARK = dict(bg="#ffffff", panel="#f5f5f8", border="#c8ccd4",
            text="#1a1a2e", dim="#6a7080", grid="#e0e2e8")
C    = ["#0077cc", "#d63031", "#e67e22", "#00875a", "#8e44ad", "#e17055"]

RISK_CLR = {"안전": "#2a9d5c", "주의": "#ffc940", "위험": "#ff7c43", "매우위험": "#e03030"}
FLUX_CAUTION = 1.0; FLUX_DANGER = 4.0; FLUX_SEVERE = 12.0
ETA_LBLS  = ["Low(친수)", "Mid(기준)", "High(소수성)"]
ETA_COEFFS = [1e-4, 5e-4, 2e-3]
ETA_C     = [C[3], C[0], C[1]]

# ─── MCHX 비말동반 상수 (McLaughlin & Webb, 2000) ─────────────────
# SAE Technical Paper 2000-01-0574
# "Wet Air Side Performance of Louver Fin Automotive Evaporators"
LP_CRIT_MCHX   = 1.2e-3    # 임계 루버 피치 [m] (1.1~1.3 mm, θ=30° 기준 중간값)
F_WET_BRIDGING  = 0.50      # 브리징 심화 시 습면 h_o 저하율 (dry 대비 50% 수준)
F_WET_MILD      = 0.85      # 일반 습면 h_o 저하율 (dry 대비 15% 저하)
F_WET_HYDRO     = 1.25      # 친수 코팅 시 h_o 향상율 (McLaughlin & Webb: +25%)


# ─── HX 스펙 ───────────────────────────────────────────

@dataclass
class FinTubeSpec:
    """
    Fin-Tube 열교환기 형상 스펙 (gap_sim + carryover 공용)

    fin_type: 'plain' | 'wavy' | 'slit' | 'louvered'
      plain    — 평판핀 (기본), Wang (2000) IJHMT 43(15)
      wavy     — 파형핀, Wang et al. (1997) IJHMT 40(4):773-784
      slit     — 랜스형 슬릿핀 (핀 연결 유지, 표면에 슬릿 가공), Wang (2001)
      louvered — 루버핀, Wang et al. (1999) IJHMT 42(1):1-17

    tube_layout: 'staggered' | 'inline'
    """
    W: float = 0.30
    H: float = 0.25
    D: float = 0.045

    fin_pitch:     float = 1.8e-3     # [m] (FPI 입력 시 자동 계산됨)
    fpi:           float = 0          # Fins Per Inch (>0이면 fin_pitch 자동 계산)
    fin_thickness: float = 0.10e-3
    fin_height:    float = 25.4e-3   # [레거시] = tube_pitch_t 와 동일. 내부에서 자동 설정됨

    tube_do:       float = 9.52e-3
    tube_di:       float = 8.52e-3
    tube_rows:     int   = 2
    tube_cols:     int   = 8
    tube_pitch_t:  float = 25.4e-3
    tube_pitch_l:  float = 22.0e-3

    fin_type:      str   = "plain"     # 'plain'|'wavy'|'slit'|'louvered'
    corr_mode:     str   = "etype"     # 'etype' (E_type 보정) | 'direct' (직접 상관식)
    tube_layout:   str   = "staggered" # 'staggered' | 'inline'

    # ── Wavy 핀 전용 ────────────────────────────────────────
    wavy_height:   float = 1.5e-3     # 파형 진폭 (peak-to-valley) [m]
    wavy_angle:    float = 17.5       # 파형 경사각 [°] (일반 15~20°)

    # ── Louvered 핀 전용 (FT용, MCHX louver와 별개) ─────────
    ft_louver_pitch: float = 1.7e-3   # 루버 피치 [m] (FT 루버: 보통 1.2~2.5mm)
    ft_louver_angle: float = 28.0     # 루버 각도 [°] (보통 20~35°)

    # ── Slit (랜스형) 핀 전용 ────────────────────────────────
    slit_num:      int   = 6          # 핀당 슬릿 수 (보통 4~10)
    slit_height:   float = 1.0e-3     # 슬릿 높이 [m] (보통 0.8~1.5mm)


    fin_material:  str   = "Al"
    tube_material: str   = "Cu"
    h_drain:       float = 0.025

    hx_type: str = field(default="FT", init=False)

    def __post_init__(self):
        self.hx_type = "FT"
        # FPI → fin_pitch 자동 변환
        if self.fpi > 0:
            self.fin_pitch = 25.4e-3 / self.fpi  # FPI → m


@dataclass
class MCHXSpec:
    """Micro-Channel 열교환기 형상 스펙 (gap_sim 전용)"""
    W: float = 0.30
    H: float = 0.20
    D: float = 0.020

    ch_width:  float = 0.8e-3
    ch_height: float = 1.5e-3
    ch_wall:   float = 0.4e-3
    n_ports:   int   = 12

    fin_pitch:     float = 1.0e-3
    fin_thickness: float = 0.06e-3
    louver_pitch:  float = 1.2e-3
    louver_angle:  float = 27.0

    n_slabs:    int   = 1
    slab_pitch: float = 8.0e-3

    fin_material:  str = "Al"
    tube_material: str = "Al"

    hx_type: str = field(default="MCHX", init=False)

    def __post_init__(self):
        self.hx_type = "MCHX"


class CarryoverSpec:
    """
    어댑터: FinTubeSpec + geo dict → carryover 모델 인터페이스 변환
    공기 주 흐름: 증발기 → Gap → 응축기
    """
    def __init__(self, ft_spec: FinTubeSpec, geo: dict):
        # FinTubeSpec 필드를 carryover 모델 변수명으로 매핑
        self.W       = ft_spec.W
        self.H       = ft_spec.H
        self.D       = ft_spec.D
        self.Fp      = ft_spec.fin_pitch
        self.df      = ft_spec.fin_thickness
        self.Do      = ft_spec.tube_do
        self.Di      = ft_spec.tube_di
        self.Pt      = ft_spec.tube_pitch_t
        self.Pl      = ft_spec.tube_pitch_l
        self.Nr      = ft_spec.tube_rows
        self.Nt      = ft_spec.tube_cols
        self.h_drain = ft_spec.h_drain
        self.tube_layout = getattr(ft_spec, 'tube_layout', 'staggered')
        self.fin_type    = getattr(ft_spec, 'fin_type', 'plain')
        # 핀 타입별 추가 파라미터 (j/f 강화 인자 계산용)
        self.wavy_height    = getattr(ft_spec, 'wavy_height', 1.5e-3)
        self.wavy_angle     = getattr(ft_spec, 'wavy_angle', 17.5)
        self.ft_louver_pitch = getattr(ft_spec, 'ft_louver_pitch', 1.7e-3)
        self.ft_louver_angle = getattr(ft_spec, 'ft_louver_angle', 28.0)
        self.slit_num       = getattr(ft_spec, 'slit_num', 6)
        self.slit_height    = getattr(ft_spec, 'slit_height', 1.0e-3)
        # 기하 계산 결과
        self.sigma   = geo['sigma']
        self.A_total = geo['A_total']
        self.A_fin   = geo['A_fin']
        self.Ac      = geo['Ac']
        self.hx_type = 'FT'


class MCHXCarryoverSpec:
    """
    어댑터: MCHXSpec + geo dict → Module B 비말동반 모델 인터페이스 변환

    [적용 상관식]
      건면 j-factor : Chang & Wang (1997) — IJHMT 40(3):533-544
        j = Re_Lp^J1 × (θ/90)^0.27 × (Fp/Lp)^-0.14 × (Fh/Lp)^-0.29
            × (Fd/Lp)^-0.23 × (Ll/Lp)^0.66 × (tp/Lp)^-0.02 × (δf/Lp)^-0.48
        J1 = -0.49 × (θ/90)^0.27
        적용범위: 100 ≤ Re_Lp ≤ 3000

      습면 보정    : McLaughlin & Webb (2000) — SAE 2000-01-0574
        Lp < Lp_crit (1.1~1.3 mm): h_o,wet ≈ 0.50 × h_o,dry  (심화 브리징)
        Lp ≥ Lp_crit              : h_o,wet ≈ 0.85 × h_o,dry  (경미한 습면)
        친수 코팅 적용 시 h_o,wet +25% (McLaughlin & Webb)

      비말 이탈 조건: Weber 수 기반 (Lp 특성길이, FT의 Fp 대신 사용)
        We_Lp = ρ_a × V_ch² × Lp / σ_w
        We_c,MCHX = 0.019 / √Lp  (FT 형태 유지, Lp 치환)
        브리징 패널티: Lp < Lp_crit 시 η_co 추가 증가
    """
    def __init__(self, mchx_spec: MCHXSpec, geo: dict):
        self.W    = mchx_spec.W
        self.H    = mchx_spec.H
        self.D    = mchx_spec.D

        # 루버 핀 기하 (Chang & Wang 1997 파라미터)
        self.Lp   = mchx_spec.louver_pitch        # 루버 피치 [m]
        self.Fp   = mchx_spec.fin_pitch            # 핀 피치 [m]
        self.Fh   = mchx_spec.slab_pitch           # 핀 높이 = 슬랩 간격 [m]
        self.Fd   = mchx_spec.D                    # 유동 깊이 = 튜브 깊이 [m]
        self.Ll   = 0.70 * mchx_spec.D             # 루버 길이 ≈ 0.70 × Fd [m]
        self.tp   = mchx_spec.slab_pitch           # 튜브 피치 [m]
        self.df   = mchx_spec.fin_thickness        # 핀 두께 [m]
        self.theta_deg = mchx_spec.louver_angle    # 루버 각도 [°]

        self.h_drain = 0.010                       # drain pan 높이 [m]

        # 기하 계산 결과
        self.sigma   = geo['sigma']
        self.A_total = geo['A_total']
        self.A_fin   = geo.get('A_louver', geo['A_total'] * 0.85)
        self.Ac      = geo['Ac']
        self.hx_type = 'MCHX'



def make_mchx_evap(base_spec: 'MCHXSpec', case: dict) -> 'MCHXSpec':
    """케이스별 louver_pitch 오버라이드가 있는 경우 새 스펙 생성"""
    if 'louver_pitch_override' not in case:
        return base_spec
    from dataclasses import replace
    return replace(base_spec, louver_pitch=case['louver_pitch_override'])


# ─── 냉매 물성 ─────────────────────────────────────────

@dataclass
class RefrigerantState:
    refrigerant: str  = "R410A"
    T_sat_evap:  float = 5.0
    T_sat_cond:  float = 45.0

    P_evap:    float = field(default=0.0, init=False)
    P_cond:    float = field(default=0.0, init=False)
    h_evap_l:  float = field(default=0.0, init=False)
    h_evap_v:  float = field(default=0.0, init=False)
    h_cond_l:  float = field(default=0.0, init=False)
    h_cond_v:  float = field(default=0.0, init=False)
    rho_l_evap: float = field(default=0.0, init=False)
    rho_v_evap: float = field(default=0.0, init=False)
    mu_l_evap:  float = field(default=0.0, init=False)
    k_l_evap:   float = field(default=0.0, init=False)
    Pr_l_evap:  float = field(default=0.0, init=False)
    h_fg_evap:  float = field(default=0.0, init=False)

    def __post_init__(self): self.compute()

    def compute(self):
        rf = self.refrigerant
        Te = self.T_sat_evap + 273.15
        Tc = self.T_sat_cond + 273.15
        self.P_evap    = CP.PropsSI('P','T',Te,'Q',0, rf)
        self.P_cond    = CP.PropsSI('P','T',Tc,'Q',0, rf)
        self.h_evap_l  = CP.PropsSI('H','T',Te,'Q',0, rf)
        self.h_evap_v  = CP.PropsSI('H','T',Te,'Q',1, rf)
        self.h_cond_l  = CP.PropsSI('H','T',Tc,'Q',0, rf)
        self.h_cond_v  = CP.PropsSI('H','T',Tc,'Q',1, rf)
        self.h_fg_evap = self.h_evap_v - self.h_evap_l
        self.rho_l_evap = CP.PropsSI('D','T',Te,'Q',0, rf)
        self.rho_v_evap = CP.PropsSI('D','T',Te,'Q',1, rf)
        self.mu_l_evap  = CP.PropsSI('V','T',Te,'Q',0, rf)
        self.k_l_evap   = CP.PropsSI('L','T',Te,'Q',0, rf)
        self.Pr_l_evap  = CP.PropsSI('Prandtl','T',Te,'Q',0, rf)


# ─── 기하 계산 ─────────────────────────────────────────

def compute_ft_geometry(spec: FinTubeSpec) -> dict:
    """Fin-Tube: 공기측 면적, 핀 효율, UA 계산

    [핀 방향] plate fin-tube HX:
    - 튜브: W 방향(가로)으로 관통
    - 핀: H×D 면 (높이×깊이), W 방향으로 Fp 간격 적층
    - Nf = W / Fp (튜브 길이 방향 핀 수)
    - CoilDesigner 검증: A_fin=0.7325 m² 일치 (±1.5%)
    """
    Nf = int(spec.W / spec.fin_pitch)    # 핀 개수 (W 방향 적층)
    Nt = spec.tube_cols                   # 튜브 열 (높이 방향)
    Nr = spec.tube_rows                   # 튜브 행 (깊이 방향)
    N_tubes = Nr * Nt                     # 총 튜브 수
    Dc = spec.tube_do + 2 * spec.fin_thickness  # collar diameter

    # 최소 유로 면적비
    sigma = (spec.fin_pitch - spec.fin_thickness) / spec.fin_pitch * \
            (spec.tube_pitch_t - Dc) / spec.tube_pitch_t
    Afr = spec.W * spec.H
    Ac  = sigma * Afr

    # ── 핀 면적 (plate fin: H × D 면, W 방향 적층) ──
    # 1개 핀: 양면 × (H×D 판 면적 - 튜브 관통 구멍, collar 포함)
    A_fin_one = 2 * (spec.H * spec.D - N_tubes * np.pi * Dc**2 / 4)
    A_fin = max(A_fin_one, 0) * Nf

    # ── 튜브 노출 면적 (핀 사이 간격) ──
    # 각 튜브: π×Dc × W(튜브 길이) × (1 - δ_fin/Fp)(핀이 점유하지 않는 비율)
    A_tube = np.pi * Dc * spec.W * (1 - spec.fin_thickness / spec.fin_pitch) * N_tubes
    A_total = A_fin + A_tube

    # 핀 효율 (초기 추정, compute_UA에서 재계산)
    k_fin = K_AL if spec.fin_material == "Al" else K_CU
    h_o   = 55.0
    m_fin = np.sqrt(2 * h_o / (k_fin * spec.fin_thickness))
    # Schmidt 등가 반경법 (원형 핀 근사)
    r_i = spec.tube_do / 2
    # Schmidt 등가 반경 (staggered plate fin)
    Xm = spec.tube_pitch_t / 2
    XL = np.sqrt((spec.tube_pitch_t / 2)**2 + spec.tube_pitch_l**2) / 2
    r_eq = max(1.27 * Xm * np.sqrt(XL / Xm - 0.3), r_i * 1.01)
    phi = (r_eq / r_i - 1) * (1 + 0.35 * np.log(r_eq / r_i))
    eta_fin = np.tanh(m_fin * r_i * phi) / (m_fin * r_i * phi + 1e-9)
    eta_o   = 1 - (A_fin / A_total) * (1 - eta_fin)

    # ── 냉매측 면적 (튜브 내면) ──
    # 각 튜브 길이 = W (폭 방향 관통)
    A_i = np.pi * spec.tube_di * spec.W * N_tubes

    # ── 벽면 열저항 ──
    k_tube = K_CU if spec.tube_material == "Cu" else K_AL
    R_wall = np.log(spec.tube_do / spec.tube_di) / \
             (2 * np.pi * k_tube * spec.W * N_tubes)

    # ── 수력직경 ──
    Lc = Nr * spec.tube_pitch_l  # 코어 깊이 (공기흐름 방향)
    Dh = 4 * Ac * Lc / (A_total + 1e-9)  # 수력직경 (Wang 2000)

    return dict(
        Ac=Ac, Afr=Afr, A_fin=A_fin, A_tube=A_tube,
        A_total=A_total, A_i=A_i,
        eta_o=eta_o, eta_fin=eta_fin,
        R_wall=R_wall, sigma=sigma,
        N_tubes=N_tubes, depth=spec.D,
        Dc=Dc, Dh=Dh, Lc=Lc,
    )


def compute_mchx_geometry(spec: MCHXSpec) -> dict:
    """MCHX: 루버 핀 공기측 면적, 채널 냉매측 면적 계산"""
    pitch_ch = spec.ch_width + spec.ch_wall
    N_ch = int(spec.W / pitch_ch)
    n_fins = int(spec.D / spec.fin_pitch)

    sigma_mchx = (spec.fin_pitch - spec.fin_thickness) / spec.fin_pitch * \
                 (spec.slab_pitch - spec.ch_height) / spec.slab_pitch
    Afr = spec.W * spec.H
    Ac  = sigma_mchx * Afr

    A_louver_one = 2 * spec.D * (spec.slab_pitch - spec.ch_height) * \
                   (spec.fin_pitch - spec.fin_thickness) / spec.fin_pitch
    A_louver = A_louver_one * (int(spec.H / spec.slab_pitch)) * n_fins

    perimeter_ch = 2 * (spec.ch_width + spec.ch_height)
    A_i = perimeter_ch * spec.H * N_ch * spec.n_slabs

    k_fin = K_AL
    h_o   = 80.0
    m_fin = np.sqrt(2 * h_o / (k_fin * spec.fin_thickness))
    L_fin = (spec.slab_pitch - spec.ch_height) / 2
    eta_fin = np.tanh(m_fin * L_fin) / (m_fin * L_fin + 1e-9)
    eta_o   = 1 - (A_louver / (A_louver + spec.W * spec.H * 0.05)) * (1 - eta_fin)

    Dh = 4 * spec.ch_width * spec.ch_height / (2 * (spec.ch_width + spec.ch_height))

    return dict(
        Ac=Ac, Afr=Afr,
        A_louver=A_louver, A_i=A_i, A_total=A_louver,
        eta_o=eta_o, eta_fin=eta_fin,
        R_wall=0.0002, sigma=sigma_mchx,
        N_ch=N_ch, Dh=Dh, depth=spec.D,
    )


# ─── FT j-factor: 4종 핀 × 2종 배열 ──────────────────────────────
#
# [참고문헌]
#   Plain:    Wang, Chi & Chang (2000) IJHMT 43(15):2693 — Part I
#   Wavy:     Wang et al. (1997) IJHMT 40(4):773-784
#   Louvered: Wang, Lee & Chang (1999) IJHMT 42(1):1-17
#   Slit:     Wang et al. (2001) / Nakayama & Xu (1983) 기반
# ─────────────────────────────────────────────────────────────────

def _j_plain_staggered(Re, spec, geo=None):
    """Plain fin, Staggered — Wang, Chi & Chang (2000) IJHMT 43(15):2693-2700
    Table 3 원논문 정확 구현 (Dc, Dh 사용)"""
    Fp = spec.fin_pitch
    Dc = geo['Dc'] if geo and 'Dc' in geo else spec.tube_do + 2*spec.fin_thickness
    Dh = geo['Dh'] if geo and 'Dh' in geo else 2e-3  # fallback
    Nr = spec.tube_rows
    Pt = spec.tube_pitch_t
    Pl = spec.tube_pitch_l
    Re = np.clip(Re, 300, 20000)
    ln_Re = np.log(Re + 1e-9)

    if Nr == 1:
        # Nr=1 case (Wang Table 3)
        p1 = -0.991 - 0.1055*(Pt/Pl)**3.1 * np.log(Re)
        p2 = -0.7344 + 2.1059*Nr**0.55 / ln_Re
        j = 0.108 * Re**(-0.29) * (Pt/Pl)**p1 * (Fp/Dc)**p2 * Nr**(-0.031)
    elif Re < 1000:
        # Nr>=2, Re<1000
        p3 = -0.361 - 0.042*Nr/ln_Re + 0.158*np.log(Nr*(Fp/Dc)**0.41)
        p4 = -1.224 - 0.076*(Pl/Dh)**1.42 / ln_Re
        p5 = -0.083 + 0.058*Nr/ln_Re
        p6 = -5.735 + 1.21*np.log(Re/Nr)
        j = 0.086 * Re**p3 * Nr**p4 * (Fp/Dc)**p5 * (Fp/Dh)**p6 * (Fp/Pt)**(-0.93)
    else:
        # Nr>=2, Re>=1000
        p3 = -0.361 - 0.042*Nr/ln_Re + 0.158*np.log(Nr*(Fp/Dc)**0.41)
        p4 = -1.224 - 0.076*(Pl/Dh)**1.42 / ln_Re
        p5 = -0.083 + 0.058*Nr/ln_Re
        p6 = -5.735 + 1.21*np.log(Re/Nr)
        j = 0.086 * Re**p3 * Nr**p4 * (Fp/Dc)**p5 * (Fp/Dh)**p6 * (Fp/Pt)**(-0.93)

    return max(0.002, j)


def _j_wavy_staggered(Re, spec, geo=None):
    """Wavy (herringbone) fin — Wang et al. (1999) IJHMT 42:1945, Table 2
    Direct 상관식. Pd = wavy_height (amplitude)."""
    Dc = geo['Dc'] if geo and 'Dc' in geo else spec.tube_do + 2*spec.fin_thickness
    Fp = spec.fin_pitch; Pt = spec.tube_pitch_t; Pl = spec.tube_pitch_l
    Nr = spec.tube_rows
    Pd = getattr(spec, 'wavy_height', 1.5e-3)  # amplitude
    Re = np.clip(Re, 300, 10000)

    if Nr == 1:
        ln_Re = np.log(Re)
        j1 = -0.383 - 0.159*np.log(max(Nr,1)*(Fp/Dc)**0.14)
        j = 1.201 / (ln_Re**2.921) * (Fp/Dc)**j1 * (Pd/Pl)**(-0.27)
    else:
        j1 = -0.329 - 0.147*np.log(Nr*(Fp/Dc)**0.44*(Pd/Dc)**0.09*(Pt/Pl)**(-0.04))
        j = 0.394 * Re**j1 * (Fp/Dc)**(-0.031) * (Pd/Dc)**0.346 * Nr**(-0.107)
    return max(0.003, j)


def _j_louvered_staggered(Re, spec, geo=None):
    """Louvered fin — Wang, Lee & Chang (1999) IJHMT 42(1):1, Table 1
    Direct 상관식. Re 입력 = Re_Dc → 내부에서 Re_Lp로 변환.
    θ는 radian 사용 (원논문 검증: degree 사용 시 j 폭발)."""
    Dc = geo['Dc'] if geo and 'Dc' in geo else spec.tube_do + 2*spec.fin_thickness
    Lp = getattr(spec, 'ft_louver_pitch', 1.7e-3)
    Fp = spec.fin_pitch; Pt = spec.tube_pitch_t; Pl = spec.tube_pitch_l
    Nr = spec.tube_rows
    theta_r = getattr(spec, 'ft_louver_angle', 28.0) * np.pi / 180
    Fl = Pl  # fin length ≈ tube_pitch_l
    Td = getattr(spec, 'D', 0.045)
    Re_Lp = np.clip(Re * Lp / Dc, 100, 5000)

    if Nr == 1:
        j1 = -0.49 + 0.021*theta_r
        j = 0.455 * Re_Lp**j1 * (Fp/Pl)**(-0.46) * (Lp/Fl)**1.14 * max(Nr,1)**(-0.34)
    else:
        j2 = -0.545 + 0.0538*theta_r - 0.0244*Nr
        j = 0.394 * Re_Lp**j2 * (Fp/Pl)**(-0.39) * (Td/Pt)**0.46 * (Lp/Fl)**0.33
    return max(0.003, j)


def _j_slit_staggered(Re, spec, geo=None):
    """Slit fin — Plain × E_slit (문헌 비율 기반)
    Wang (2001) direct f 상관식 불안정 → etype 비율 사용.
    Yan & Sheen (2000): j_slit/j_plain ≈ 1.3~1.4"""
    j_plain = _j_plain_staggered(Re, spec, geo)
    Ns = getattr(spec, 'slit_num', 6)
    Sh = getattr(spec, 'slit_height', 1e-3)
    Fp = spec.fin_pitch
    E_s = 1.15 + 0.12 * min(Ns, 10)**0.30 * (Sh / Fp)**0.20
    return max(0.003, j_plain * np.clip(E_s, 1.20, 1.50))


def _j_ft_staggered(Re, spec, geo=None):
    """FT Staggered j-factor — 핀 타입(4종) × 상관식 모드(2종) 분기"""
    ft = getattr(spec, 'fin_type', 'plain').lower()
    mode = getattr(spec, 'corr_mode', 'etype').lower()

    if mode == 'direct':
        return _j_ft_staggered_direct(Re, spec, ft, geo)
    else:
        if ft == 'wavy':      return _j_wavy_staggered(Re, spec, geo)
        elif ft in ('louvered', 'louver'): return _j_louvered_staggered(Re, spec, geo)
        elif ft == 'slit':    return _j_slit_staggered(Re, spec, geo)
        else:                 return _j_plain_staggered(Re, spec, geo)


def _j_ft_staggered_direct(Re, spec, ft, geo=None):
    """직접 상관식 모드 — 핀별 독립 상관식"""
    import warnings

    if ft == 'plain':
        return _j_plain_staggered(Re, spec, geo)
    elif ft == 'wavy':
        warnings.warn("Wavy direct 미구현 → etype fallback", stacklevel=3)
        return _j_wavy_staggered(Re, spec, geo)
    elif ft in ('louvered', 'louver'):
        warnings.warn("Louvered direct 미구현 → etype fallback", stacklevel=3)
        return _j_louvered_staggered(Re, spec, geo)
    elif ft == 'slit':
        warnings.warn("Slit direct 미구현 → etype fallback", stacklevel=3)
        return _j_slit_staggered(Re, spec, geo)
    else:
        return _j_plain_staggered(Re, spec, geo)


def _j_ft_inline(Re, spec):
    """FT Inline j-factor — wake shielding × 핀 강화 (독립 작용)"""
    Nr = spec.tube_rows
    j_stag = _j_ft_staggered(Re, spec)
    if Nr < 2:   return j_stag
    elif Nr == 2: return j_stag * 0.82
    else:         return j_stag * 0.75


def air_htc_ft(spec: FinTubeSpec, geo: dict, V_face: float) -> float:
    """FT 공기측 HTC — Wang(2000) 정석, Dc 기반 Re"""
    rho_air=1.18; mu_air=1.85e-5; cp_air=1006; Pr_air=0.71
    G   = rho_air * V_face / geo['sigma']
    Dc  = geo.get('Dc', spec.tube_do + 2*spec.fin_thickness)
    Re  = G * Dc / mu_air
    layout = getattr(spec, 'tube_layout', 'staggered').lower()
    if layout in ('inline', 'parallel'):
        j = _j_ft_inline(Re, spec)
    else:
        j = _j_ft_staggered(Re, spec, geo)
    return j * G * cp_air / Pr_air**(2/3)


def air_htc_mchx(spec: MCHXSpec, geo: dict, V_face: float) -> float:
    """MCHX 루버 핀 공기측 HTC — Chang & Wang (1997)"""
    rho_air=1.18; mu_air=1.85e-5; cp_air=1006; Pr_air=0.71
    G      = rho_air * V_face / geo['sigma']
    Re_lp  = G * spec.louver_pitch * 1e-3 / mu_air
    alpha  = spec.louver_angle * np.pi / 180
    j = 0.5 * Re_lp**(-0.49) * (alpha/90)**0.27 * \
        (spec.fin_pitch/spec.louver_pitch)**(-0.14)
    j = max(j, 0.004)
    return j * G * cp_air / Pr_air**(2/3)


def refrigerant_htc_evap(spec, geo: dict, ref: RefrigerantState, Q_evap: float) -> float:
    """증발기 냉매측 HTC — Shah (1982) 유동비등"""
    if hasattr(spec, 'tube_di'):
        D_ch = spec.tube_di
        A_cs = np.pi * D_ch**2 / 4 * geo['N_tubes']
    else:
        D_ch = geo['Dh']
        A_cs = spec.ch_width * spec.ch_height * geo['N_ch'] * spec.n_slabs

    x_mean = 0.5
    G_ref  = max(Q_evap / (ref.h_fg_evap * A_cs), 50)
    Re_l   = max(G_ref * (1 - x_mean) * D_ch / ref.mu_l_evap, 2300)
    Nu_l   = 0.023 * Re_l**0.8 * ref.Pr_l_evap**0.4
    h_l    = Nu_l * ref.k_l_evap / D_ch
    Co = ((1-x_mean)/x_mean)**0.8 * (ref.rho_v_evap/ref.rho_l_evap)**0.5
    Bo = max(Q_evap / (geo['A_i'] * G_ref * ref.h_fg_evap + 1e-3), 1e-5)
    if Co > 0.65:
        psi = 1.8 / Co**0.8
    else:
        psi = max(1.8/Co**0.8, 0.6683*Co**(-0.2) + 1058*Bo**0.7)
    return h_l * psi


def refrigerant_htc_cond(spec, geo: dict, ref: RefrigerantState, Q_cond: float) -> float:
    """응축기 냉매측 HTC — Shah (1979) 응축"""
    if hasattr(spec, 'tube_di'):
        D_ch = spec.tube_di
        A_cs = np.pi * D_ch**2 / 4 * geo['N_tubes']
    else:
        D_ch = geo['Dh']
        A_cs = spec.ch_width * spec.ch_height * geo['N_ch'] * spec.n_slabs

    x_mean = 0.5
    Tc   = ref.T_sat_cond + 273.15
    mu_l = CP.PropsSI('V','T',Tc,'Q',0, ref.refrigerant)
    k_l  = CP.PropsSI('L','T',Tc,'Q',0, ref.refrigerant)
    Pr_l = CP.PropsSI('Prandtl','T',Tc,'Q',0, ref.refrigerant)
    G_ref = max(Q_cond / ((ref.h_cond_v - ref.h_cond_l) * A_cs), 50)
    Re_l  = G_ref * D_ch / mu_l
    h_l   = 0.023 * Re_l**0.8 * Pr_l**0.4 * k_l / D_ch
    Z     = (1/x_mean - 1)**0.8 * (ref.P_cond / CP.PropsSI('Pcrit', ref.refrigerant))**0.4
    return h_l * (1 + 3.8 / Z**0.95)


# ─── 단상 HTC: Gnielinski (1976) ────────────────────────────────

def gnielinski_htc(Re, Pr, k, D):
    """Gnielinski (1976) 단상 관내 HTC.
    유효범위: 2300 < Re < 5e6, 0.5 < Pr < 2000
    과열 증기 / 과냉 액체 공통."""
    Re = max(Re, 2300)
    f = (0.790 * np.log(Re) - 1.64)**(-2)
    Nu = (f/8) * (Re - 1000) * Pr / (1 + 12.7 * np.sqrt(f/8) * (Pr**(2/3) - 1))
    Nu = max(Nu, 3.66)  # 층류 하한
    return Nu * k / D


def refrigerant_htc_auto(spec, geo, ref, side, x, m_ref, Q_seg=None):
    """건도(x) 기반 냉매 HTC 자동 분기

    x < 0:    과냉 액체 → Gnielinski(liquid)
    0 ≤ x ≤ 1: 이상 영역 → Shah(1982) evap / Shah(1979) cond
    x > 1:    과열 증기 → Gnielinski(vapor)
    전이구간: 블렌딩

    Shah(1982) — 정석 구현:
      h_tp = ψ × h_lo
      h_lo = Dittus-Boelter with Re_lo = G×D/μ_l (total mass flux)
      ψ = f(Co, Bo, Fr_l) — 3-regime chart correlation
    """
    rf = ref.refrigerant
    if hasattr(spec, 'tube_di'):
        D = spec.tube_di
        A_cs_one = np.pi * D**2 / 4
    else:
        D = geo.get('Dh', 2e-3)
        A_cs_one = spec.ch_width * spec.ch_height * getattr(spec, 'n_ports', 10)

    G_ref = max(m_ref / A_cs_one, 50.0)

    if side == 'evap':
        T_sat_K = ref.T_sat_evap + 273.15
        P_sat = ref.P_evap; h_fg = ref.h_fg_evap
    else:
        T_sat_K = ref.T_sat_cond + 273.15
        P_sat = ref.P_cond; h_fg = ref.h_cond_v - ref.h_cond_l

    # ── 과냉 (x < -0.05) ──
    if x < -0.05:
        mu_l = CP.PropsSI('V','T',T_sat_K,'Q',0,rf)
        k_l  = CP.PropsSI('L','T',T_sat_K,'Q',0,rf)
        Pr_l = CP.PropsSI('Prandtl','T',T_sat_K,'Q',0,rf)
        return gnielinski_htc(G_ref*D/mu_l, Pr_l, k_l, D), 'subcooled'

    # ── 과열 (x > 1.05) ──
    if x > 1.05:
        mu_v = CP.PropsSI('V','T',T_sat_K,'Q',1,rf)
        k_v  = CP.PropsSI('L','T',T_sat_K,'Q',1,rf)
        Pr_v = CP.PropsSI('Prandtl','T',T_sat_K,'Q',1,rf)
        return gnielinski_htc(G_ref*D/mu_v, Pr_v, k_v, D), 'superheated'

    # ── 이상 영역 (0 ≤ x ≤ 1) ──
    x_clip = np.clip(x, 0.05, 0.95)

    # 공통 물성
    mu_l = ref.mu_l_evap if side == 'evap' else CP.PropsSI('V','T',T_sat_K,'Q',0,rf)
    k_l  = ref.k_l_evap  if side == 'evap' else CP.PropsSI('L','T',T_sat_K,'Q',0,rf)
    Pr_l = ref.Pr_l_evap if side == 'evap' else CP.PropsSI('Prandtl','T',T_sat_K,'Q',0,rf)
    rho_l = ref.rho_l_evap if side == 'evap' else CP.PropsSI('D','T',T_sat_K,'Q',0,rf)
    rho_v = ref.rho_v_evap if side == 'evap' else CP.PropsSI('D','T',T_sat_K,'Q',1,rf)

    # h_lo: 전체 유량이 액체일 때 (Shah 기준)
    Re_lo = G_ref * D / mu_l
    h_lo = 0.023 * max(Re_lo, 2300)**0.8 * Pr_l**0.4 * k_l / D

    if side == 'evap':
        # ── Shah (1982) 증발 — Chart Correlation ──
        # Convection number
        Co = ((1-x_clip)/x_clip)**0.8 * (rho_v/rho_l)**0.5

        # Boiling number — 정석: q'' / (G × h_fg)
        Nr = max(getattr(spec, 'tube_rows', 2), 1)
        Nt = max(getattr(spec, 'tube_cols', 1), 1)
        N_seg = 10  # 기본 세그먼트 수
        A_i_seg = geo['A_i'] / (Nr * Nt * N_seg)
        q_flux = max((Q_seg or 15.0) / (A_i_seg + 1e-9), 100.0)
        Bo = q_flux / (G_ref * h_fg)

        # Froude number (수평 튜브)
        Fr_l = G_ref**2 / (rho_l**2 * 9.81 * D)

        # N parameter (Shah 1982)
        if Fr_l >= 0.04:
            N = Co
        else:
            N = Co * 0.38 * Fr_l**(-0.3)

        # ψ (Chart correlation, Shah 1982 Table 2)
        # Regime: CB (convective boiling) vs NB (nucleate boiling)
        psi_cb = 1.8 / max(N, 0.01)**0.8

        if N > 1.0:
            # Nucleate boiling dominant regime
            if Bo > 0.3e-4:
                psi_nb = 230 * Bo**0.5
            else:
                psi_nb = 1.0 + 46 * Bo**0.5
            psi = max(psi_nb, psi_cb)
        else:
            # Convective boiling dominant regime (N ≤ 1.0)
            if 0.1 < N <= 1.0:
                F_s = 14.7 if Bo >= 11e-4 else 15.43
                psi_bs = F_s * Bo**0.5 * np.exp(2.74 * N**(-0.1))
            else:  # N ≤ 0.1
                F_s = 14.7 if Bo >= 11e-4 else 15.43
                psi_bs = F_s * Bo**0.5 * np.exp(2.74 * N**(-0.15))
            psi = max(psi_cb, psi_bs)

        h_2ph = h_lo * psi

    else:
        # ── Shah (1979) 응축 ──
        Re_lo = G_ref * D / mu_l
        h_lo = 0.023 * max(Re_lo, 2300)**0.8 * Pr_l**0.4 * k_l / D
        Z = (1/x_clip - 1)**0.8 * (P_sat / CP.PropsSI('Pcrit', rf))**0.4
        h_2ph = h_lo * (1 + 3.8 / max(Z, 0.01)**0.95)

    # ── 전이 블렌딩 (x = 0.90~1.05) ──
    if 0.90 < x <= 1.05:
        mu_v = CP.PropsSI('V','T',T_sat_K,'Q',1,rf)
        k_v  = CP.PropsSI('L','T',T_sat_K,'Q',1,rf)
        Pr_v = CP.PropsSI('Prandtl','T',T_sat_K,'Q',1,rf)
        h_vapor = gnielinski_htc(G_ref*D/mu_v, Pr_v, k_v, D)
        w = np.clip((x - 0.90) / 0.15, 0, 1)
        return (1 - w) * h_2ph + w * h_vapor, 'transition'

    # ── 과냉 전이 (-0.05 ~ 0) ──
    if -0.05 <= x < 0:
        h_sub = gnielinski_htc(Re_lo, Pr_l, k_l, D)
        w = np.clip((x + 0.05) / 0.05, 0, 1)
        return (1 - w) * h_sub + w * h_2ph, 'transition'

    return h_2ph, 'two_phase'


def _update_eta_o(spec, geo: dict, h_o: float) -> float:
    """
    실제 h_o 기반 핀 효율 재계산 (반복 수렴)

    기존: h_o=55 고정으로 η_fin 계산 → wavy/louvered에서 h_o=160일 때 η_o 과대
    수정: 실제 h_o로 η_fin을 재계산

    m = √(2×h_o / (k_fin × δ_fin))
    η_fin = tanh(m×L) / (m×L)
    η_o = 1 - (A_fin/A_total) × (1 - η_fin)
    """
    if spec.hx_type != "FT":
        return geo['eta_o']  # MCHX는 별도 처리

    k_fin = K_AL if spec.fin_material == "Al" else K_CU
    m_fin = np.sqrt(2 * h_o / (k_fin * spec.fin_thickness))
    # Schmidt 등가 반경법
    r_i = spec.tube_do / 2
    # Schmidt 등가 반경 (staggered plate fin)
    Xm = spec.tube_pitch_t / 2
    XL = np.sqrt((spec.tube_pitch_t / 2)**2 + spec.tube_pitch_l**2) / 2
    r_eq = max(1.27 * Xm * np.sqrt(XL / Xm - 0.3), r_i * 1.01)
    phi = (r_eq / r_i - 1) * (1 + 0.35 * np.log(r_eq / r_i))
    mL = m_fin * r_i * phi
    if mL > 0.01:
        eta_fin = np.tanh(mL) / mL
    else:
        eta_fin = 1.0
    eta_o = 1 - (geo.get('A_fin', geo.get('A_louver', geo['A_total']*0.85)) / geo['A_total']) * (1 - eta_fin)
    return max(eta_o, 0.3)


def compute_UA(spec, geo: dict, ref: RefrigerantState, V_face: float, side: str) -> dict:
    """
    전체 UA [W/K] 계산: 1/UA = R_o + R_wall + R_i
    side: 'evap' or 'cond'

    [v3 보정] h_o ↔ η_o 반복 수렴 (최대 5회)
    기존: geo['eta_o'] (h_o=55 기반 고정) → h_o가 높은 핀에서 과대
    수정: 실제 h_o로 η_o 재계산 → UA 정확도 개선
    """
    Q_ref = 3000.0
    if spec.hx_type == "FT":
        h_o = air_htc_ft(spec, geo, V_face)
        h_i = refrigerant_htc_evap(spec, geo, ref, Q_ref) if side=='evap' \
              else refrigerant_htc_cond(spec, geo, ref, Q_ref)
    else:
        h_o = air_htc_mchx(spec, geo, V_face)
        h_i = refrigerant_htc_evap(spec, geo, ref, Q_ref) if side=='evap' \
              else refrigerant_htc_cond(spec, geo, ref, Q_ref)

    # η_o 반복 수렴 (FT만, MCHX는 geo 값 유지)
    eta_o = geo['eta_o']
    for _ in range(5):
        eta_o_new = _update_eta_o(spec, geo, h_o)
        if abs(eta_o_new - eta_o) < 0.001:
            break
        eta_o = eta_o_new

    R_o    = 1.0 / (eta_o * h_o * geo['A_total'])
    R_wall = geo['R_wall']
    R_i    = 1.0 / (h_i * geo['A_i'])
    UA     = 1.0 / (R_o + R_wall + R_i)
    return dict(UA=UA, h_o=h_o, h_i=h_i, R_o=R_o, R_wall=R_wall, R_i=R_i,
                eta_o=eta_o, eta_o_init=geo['eta_o'],
                A_o=geo['A_total'], A_i=geo['A_i'])


def compute_coil_performance(spec, geo: dict, ua_result: dict,
                             T_in: float, RH_in: float, T_wall: float,
                             V_face: float) -> dict:
    """
    공통 Threlkeld 습면/건면 코일 성능 계산 — 정석 구현

    ASHRAE Fundamentals Ch.23 / Threlkeld 기반
    모든 모듈(A/B/C)이 공유하는 단일 열전달 엔진

    Parameters
    ----------
    spec    : FinTubeSpec or MCHXSpec
    geo     : compute_ft_geometry() 결과
    ua_result : compute_UA() 결과 (h_o, h_i, eta_o, R_wall, R_i, UA, A_o, A_i)
    T_in    : 입구 공기 온도 [°C]
    RH_in   : 입구 상대 습도 [-]
    T_wall  : 코일 표면 온도 [°C] (= T_sat_evap)
    V_face  : 전면 풍속 [m/s]

    Returns
    -------
    dict with: Q_total, Q_sen, Q_lat, SHR, T_out, W_in, W_out, T_dp,
               q_cond, dehumid_rate, eps, NTU, UA_eff, eta_o_eff,
               h_in, h_out, m_dot, is_wet
    """
    rho_air = 1.18; cp_air = 1006.0
    P_atm = 101325.0

    # ── STEP 1: 공통 물성 ──
    A_face = spec.W * spec.H
    m_dot = rho_air * V_face * A_face

    try:
        from CoolProp.HumidAirProp import HAPropsSI
        _CP = True
    except ImportError:
        _CP = False

    if _CP:
        T_in_K = T_in + 273.15
        T_w_K  = T_wall + 273.15
        W_in     = HAPropsSI('W', 'T', T_in_K, 'R', RH_in, 'P', P_atm)
        h_in     = HAPropsSI('H', 'T', T_in_K, 'R', RH_in, 'P', P_atm)
        T_dp     = HAPropsSI('D', 'T', T_in_K, 'R', RH_in, 'P', P_atm) - 273.15
        W_sat_w  = HAPropsSI('W', 'T', T_w_K, 'R', 1.0, 'P', P_atm)
        h_sat_w  = HAPropsSI('H', 'T', T_w_K, 'R', 1.0, 'P', P_atm)
    else:
        def _ps(T):
            Tk=T+273.15
            return np.exp(-5.8002206e3/Tk+1.3914993-4.864e-2*Tk+4.1765e-5*Tk**2-1.4452e-8*Tk**3+6.546*np.log(Tk))
        def _W(T,RH): Pw=RH*_ps(T); return 0.62198*Pw/(P_atm-Pw)
        def _h(T,W): return 1006.0*T+W*(2501000.0+1860.0*T)
        def _Td(T,RH):
            Pw=RH*_ps(T)/1000.0
            if Pw<=0: return -40
            lp=np.log(Pw); return max(-40,6.54+14.526*lp+0.7389*lp**2+0.09486*lp**3+0.4569*Pw**0.1984)
        W_in    = _W(T_in, RH_in)
        h_in    = _h(T_in, W_in)
        T_dp    = _Td(T_in, RH_in)
        W_sat_w = _W(T_wall, 1.0)
        h_sat_w = _h(T_wall, W_sat_w)

    # ── STEP 2: 건면 UA, NTU, ε (compute_UA 결과 재사용) ──
    h_o    = ua_result['h_o']
    h_i    = ua_result['h_i']
    eta_o  = ua_result['eta_o']
    UA_dry = ua_result['UA']
    R_wall = ua_result['R_wall']
    R_i    = ua_result['R_i']
    A_total = ua_result['A_o']
    A_i     = ua_result['A_i']

    C_air   = m_dot * cp_air
    NTU_dry = UA_dry / C_air
    eps_dry = 1.0 - np.exp(-NTU_dry)

    # ── STEP 3: 건면/습면 분기 ──
    if T_wall >= T_dp:
        # ────── 건면 ──────
        eps = eps_dry
        NTU = NTU_dry
        UA_eff = UA_dry
        eta_o_eff = eta_o

        Q_total = eps * C_air * (T_in - T_wall)
        Q_total = max(Q_total, 0.0)
        T_out = T_in - eps * (T_in - T_wall)
        W_out = W_in
        Q_lat = 0.0
        Q_sen = Q_total
        SHR   = 1.0
        q_cond = 0.0
        is_wet = False

        if _CP:
            h_out = HAPropsSI('H', 'T', T_out+273.15, 'W', W_out, 'P', P_atm)
        else:
            h_out = _h(T_out, W_out)

    else:
        # ────── 습면: Threlkeld 정석 ──────
        h_fg = 2501000.0 - 2381.0 * T_wall

        # ξ: 습면 보정 계수 (T_mean에서 평가)
        T_m = (T_in + T_wall) / 2.0
        if _CP:
            dT = 0.5
            W_sp = HAPropsSI('W', 'T', (T_m+dT)+273.15, 'R', 1.0, 'P', P_atm)
            W_sm = HAPropsSI('W', 'T', (T_m-dT)+273.15, 'R', 1.0, 'P', P_atm)
            dWs_dT_m = (W_sp - W_sm) / (2*dT)
            W_sp2 = HAPropsSI('W', 'T', (T_wall+dT)+273.15, 'R', 1.0, 'P', P_atm)
            W_sm2 = HAPropsSI('W', 'T', (T_wall-dT)+273.15, 'R', 1.0, 'P', P_atm)
            dWs_dT_w = (W_sp2 - W_sm2) / (2*dT)
        else:
            dT = 0.5
            dWs_dT_m = (_W(T_m+dT,1.0) - _W(T_m-dT,1.0)) / (2*dT)
            dWs_dT_w = (_W(T_wall+dT,1.0) - _W(T_wall-dT,1.0)) / (2*dT)

        xi = h_fg * dWs_dT_m / cp_air

        # 습면 핀효율
        k_fin = K_AL if getattr(spec, 'fin_material', 'Al') == 'Al' else K_CU
        df = spec.fin_thickness if hasattr(spec, 'fin_thickness') else 0.1e-3
        m_dry = np.sqrt(2 * h_o / (k_fin * df))
        m_wet = m_dry * np.sqrt(1 + max(xi, 0))

        # Schmidt 등가반경 (compute_ft_geometry와 동일)
        if hasattr(spec, 'tube_do'):
            r_i = spec.tube_do / 2
            Pt = spec.tube_pitch_t if hasattr(spec, 'tube_pitch_t') else 25.4e-3
            Pl = spec.tube_pitch_l if hasattr(spec, 'tube_pitch_l') else 22.0e-3
            Xm = Pt / 2
            XL = np.sqrt((Pt/2)**2 + Pl**2) / 2
            r_eq = max(1.27 * Xm * np.sqrt(XL / Xm - 0.3), r_i * 1.01)
            phi = (r_eq / r_i - 1) * (1 + 0.35 * np.log(r_eq / r_i))
            arg_wet = m_wet * r_i * phi
        else:
            # MCHX fallback
            L_fin = getattr(spec, 'slab_pitch', 8e-3) / 2
            arg_wet = m_wet * L_fin

        eta_fin_wet = np.tanh(arg_wet) / (arg_wet + 1e-9)
        A_fin = geo.get('A_fin', geo.get('A_louver', A_total * 0.85))
        eta_o_wet = 1.0 - (A_fin / A_total) * (1.0 - eta_fin_wet)
        eta_o_wet = max(eta_o_wet, 0.1)

        # 습면 유효비열
        cp_s = cp_air + h_fg * dWs_dT_w
        b = cp_s / cp_air

        # 습면 UA (R_i, R_wall 포함)
        UA_o_wet = b * eta_o_wet * h_o * A_total
        UA_wet = 1.0 / (1.0 / UA_o_wet + R_wall + R_i)

        # 습면 ε-NTU — ★ NTU = UA / (m_dot × cp_a) ★
        NTU_wet = UA_wet / C_air
        eps_wet = 1.0 - np.exp(-NTU_wet)

        eps = eps_wet
        NTU = NTU_wet
        UA_eff = UA_wet
        eta_o_eff = eta_o_wet

        # Q (엔탈피 기반)
        Q_total = max(eps * m_dot * (h_in - h_sat_w), 0.0)

        # W_out (Lewis 유사)
        W_out = W_in - eps * (W_in - W_sat_w)
        W_out = max(W_out, W_sat_w)

        # T_out
        T_out = T_in - eps * (T_in - T_wall)

        # Q 분해
        Q_lat = m_dot * (W_in - W_out) * h_fg
        Q_sen = max(Q_total - Q_lat, 0.0)
        SHR = Q_sen / (Q_total + 1e-9)

        # 응결수 [g/s]
        q_cond = m_dot * (W_in - W_out) * 1000.0
        is_wet = True

        if _CP:
            h_out = HAPropsSI('H', 'T', T_out+273.15, 'W', min(W_out, W_in), 'P', P_atm)
        else:
            h_out = _h(T_out, W_out)

    # 제습량 [g/h]
    dehumid_rate = q_cond * 3600.0

    return dict(
        Q_total=Q_total, Q_sen=Q_sen, Q_lat=Q_lat, SHR=SHR,
        T_out=T_out, W_in=W_in, W_out=W_out, T_dp=T_dp,
        q_cond=q_cond, dehumid_rate=dehumid_rate,
        eps=eps, NTU=NTU, UA_eff=UA_eff, eta_o_eff=eta_o_eff,
        h_in=h_in, h_out=h_out, h_sat_w=h_sat_w,
        m_dot=m_dot, is_wet=is_wet,
        # 건면 참조값
        eps_dry=eps_dry, NTU_dry=NTU_dry, UA_dry=UA_dry,
    )


def _split_geo_per_row(geo: dict, Nr: int) -> dict:
    """기하 데이터를 Row당 분할 — 직렬 연결. FT/MCHX 공용."""
    A_fin = geo.get('A_fin', geo.get('A_louver', geo['A_total'] * 0.85))
    A_tube = geo.get('A_tube', geo['A_total'] - A_fin)
    N_tubes = geo.get('N_tubes', geo.get('N_ch', 1))
    return dict(
        Ac=geo['Ac'],                     # 유로 단면적: Row에 무관
        Afr=geo['Afr'],                   # 전면 면적: 동일
        A_fin=A_fin / Nr,                 # 핀 면적: Row당 1/Nr
        A_tube=A_tube / Nr,               # 튜브 노출 면적: 1/Nr
        A_total=geo['A_total'] / Nr,      # 총 공기측 면적: 1/Nr
        A_i=geo['A_i'] / Nr,              # 냉매측 면적: 1/Nr
        eta_o=geo['eta_o'],               # 핀효율: Row 무관
        eta_fin=geo['eta_fin'],
        R_wall=geo['R_wall'] * Nr,        # 벽면 열저항: 직렬이므로 ×Nr
        sigma=geo['sigma'],               # 유로비: 동일
        N_tubes=N_tubes // Nr,            # Row당 튜브 수
        depth=geo.get('depth', 0.045) / Nr,
    )


def _split_ua_per_row(ua: dict, Nr: int) -> dict:
    """UA 데이터를 Row당 분할"""
    return dict(
        UA=ua['UA'] / Nr,                 # 직렬 → UA ∝ A → 1/Nr
        h_o=ua['h_o'],                    # h_o: Re에만 의존 → 동일
        h_i=ua['h_i'],                    # h_i: 냉매측 → 동일
        R_o=ua['R_o'] * Nr,              # R_o = 1/(η_o×h_o×A) → ×Nr
        R_wall=ua['R_wall'] * Nr,        # 직렬 → ×Nr
        R_i=ua['R_i'] * Nr,              # R_i = 1/(h_i×A_i) → ×Nr
        eta_o=ua['eta_o'],               # 핀효율: 동일
        eta_o_init=ua.get('eta_o_init', ua['eta_o']),
        A_o=ua['A_o'] / Nr,             # 공기측 면적: 1/Nr
        A_i=ua['A_i'] / Nr,             # 냉매측 면적: 1/Nr
    )


def compute_coil_performance_segmented(spec, geo: dict, ua_result: dict,
                                        T_in: float, RH_in: float, T_wall: float,
                                        V_face: float) -> dict:
    """
    Tube-Row Segmented 코일 성능 계산

    증발기를 Nr개 Row로 분할하여 순차 계산.
    각 Row 출구 상태가 다음 Row 입구가 됨.

    장점: 부분 습면 전환 (앞열 습면 → 뒷열 건면) 포착
    비용: compute_coil_performance × Nr 회 (Nr=2~3, 무시 가능)

    인터페이스는 compute_coil_performance()과 동일 → drop-in 교체 가능
    """
    Nr = getattr(spec, 'tube_rows', None)
    if Nr is None:
        Nr = getattr(spec, 'n_slabs', 1)  # MCHX: n_slabs 사용
    if Nr <= 1:
        return compute_coil_performance(spec, geo, ua_result,
                                        T_in, RH_in, T_wall, V_face)

    geo_row = _split_geo_per_row(geo, Nr)
    ua_row  = _split_ua_per_row(ua_result, Nr)

    P_atm = 101325.0
    try:
        from CoolProp.HumidAirProp import HAPropsSI
        _CP = True
    except ImportError:
        _CP = False

    # Row별 순차 계산
    T_air = T_in
    RH_air = RH_in
    rows = []

    for i in range(Nr):
        row = compute_coil_performance(spec, geo_row, ua_row,
                                       T_air, RH_air, T_wall, V_face)
        rows.append(row)

        # 다음 Row 입구 = 현재 Row 출구
        T_air = row['T_out']
        W_air = row['W_out']

        # RH 역산 (다음 Row 입구)
        if _CP:
            try:
                RH_air = HAPropsSI('R', 'T', T_air + 273.15, 'W', W_air, 'P', P_atm)
                RH_air = min(max(RH_air, 0.01), 1.0)
            except Exception:
                RH_air = min(W_air / (0.62198 * 611.2 * np.exp(17.67 * T_air / (T_air + 243.5)) / (P_atm - 611.2 * np.exp(17.67 * T_air / (T_air + 243.5))) + 1e-9), 1.0)
        else:
            def _ps(T):
                Tk = T + 273.15
                return np.exp(-5.8002206e3/Tk + 1.3914993 - 4.864e-2*Tk + 4.1765e-5*Tk**2 - 1.4452e-8*Tk**3 + 6.546*np.log(Tk))
            Pw = W_air * P_atm / (0.62198 + W_air)
            RH_air = min(max(Pw / _ps(T_air), 0.01), 1.0)

    # ── 전체 결과 합산 ──
    Q_total = sum(r['Q_total'] for r in rows)
    Q_sen   = sum(r['Q_sen'] for r in rows)
    Q_lat   = sum(r['Q_lat'] for r in rows)
    q_cond  = sum(r['q_cond'] for r in rows)
    SHR     = Q_sen / (Q_total + 1e-9)
    dehumid = q_cond * 3600.0

    # 출구 상태: 마지막 Row 출구
    last = rows[-1]
    first = rows[0]

    # 전체 등가 ε: Q_total / Q_max
    rho_air = 1.18; cp_air = 1006.0
    m_dot = rho_air * V_face * spec.W * spec.H
    if _CP:
        h_in = HAPropsSI('H', 'T', T_in + 273.15, 'R', RH_in, 'P', P_atm)
        W_sat_w = HAPropsSI('W', 'T', T_wall + 273.15, 'R', 1.0, 'P', P_atm)
        h_sat_w = HAPropsSI('H', 'T', T_wall + 273.15, 'R', 1.0, 'P', P_atm)
    else:
        def _W(T, RH):
            Pw = RH * _ps(T); return 0.62198 * Pw / (P_atm - Pw)
        def _h(T, W): return 1006.0*T + W*(2501000.0 + 1860.0*T)
        h_in = first['h_in']
        W_sat_w = _W(T_wall, 1.0)
        h_sat_w = _h(T_wall, W_sat_w)

    Q_max = m_dot * (h_in - h_sat_w) if (h_in - h_sat_w) > 0 else 1.0
    eps_eff = Q_total / Q_max if Q_max > 0 else 0.0
    eps_eff = min(eps_eff, 0.999)
    NTU_eff = -np.log(1 - eps_eff) if eps_eff < 0.999 else 7.0

    # Row별 습면/건면 정보
    wet_rows = sum(1 for r in rows if r['is_wet'])

    return dict(
        Q_total=Q_total, Q_sen=Q_sen, Q_lat=Q_lat, SHR=SHR,
        T_out=last['T_out'], W_in=first['W_in'], W_out=last['W_out'],
        T_dp=first['T_dp'],
        q_cond=q_cond, dehumid_rate=dehumid,
        eps=eps_eff, NTU=NTU_eff,
        UA_eff=ua_result['UA'], eta_o_eff=ua_result['eta_o'],
        h_in=first['h_in'], h_out=last['h_out'],
        h_sat_w=h_sat_w if '_CP' in dir() else first.get('h_sat_w', 0),
        m_dot=m_dot, is_wet=(wet_rows > 0),
        # 건면 참조값 (Lumped 기준)
        eps_dry=first['eps_dry'], NTU_dry=first['NTU_dry'] * Nr,
        UA_dry=ua_result['UA'],
        # Segment 상세
        rows=rows, wet_rows=wet_rows, Nr=Nr,
    )


# ═══════════════════════════════════════════════════════════════════
#  Level 1 코일 모델 — Row × Phase-zone 세그먼트 (건도 추적)
# ═══════════════════════════════════════════════════════════════════

def compute_coil_v2(spec, geo, ref, T_air_in, RH_in, V_face,
                    m_ref, x_in, side='evap'):
    """
    Level 1 Row-by-Row + Phase-zone 세그먼트 계산 (v2.1)

    [v2.0 대비 개선]
    - 공기밀도 ρ(T) 보정 (1.18 고정 → 이상기체)
    - b-factor 상한 3.5 (습면 UA 과대 방지)
    - Q_row 상한 = dx_max × m_ref × h_fg (Row당 건도 변화 제한)
      → Row간 열량 분배 정상화, 과열 진입 Row 특정

    Parameters
    ----------
    spec, geo, ref : 기하/냉매 스펙
    T_air_in, RH_in, V_face : 공기 조건
    m_ref : 냉매 질량유량 [kg/s]
    x_in  : 입구 건도 (0~1: 이상, >1: 과열, <0: 과냉)
    side  : 'evap' or 'cond'
    """
    Nr = getattr(spec, 'tube_rows', getattr(spec, 'n_slabs', 1))
    Nr = max(Nr, 1)
    Nt = getattr(spec, 'tube_cols', 1)

    cp_air = 1006.0; P_atm = 101325.0
    h_fg_lat = 2501000.0
    B_MAX = 3.5          # b-factor 상한
    DX_MAX_PER_TUBE = 0.12  # 튜브 1개당 최대 건도 변화

    try:
        from CoolProp.HumidAirProp import HAPropsSI
        _CP = True
    except ImportError:
        _CP = False

    # ── 냉매 물성 ──
    rf = ref.refrigerant
    if side == 'evap':
        T_sat = ref.T_sat_evap; T_sat_K = T_sat + 273.15
        h_fg = ref.h_fg_evap
    else:
        T_sat = ref.T_sat_cond; T_sat_K = T_sat + 273.15
        h_fg = ref.h_cond_v - ref.h_cond_l

    cp_v = CP.PropsSI('C', 'T', T_sat_K, 'Q', 1, rf)
    cp_l = CP.PropsSI('C', 'T', T_sat_K, 'Q', 0, rf)

    # Row당 최대 Q (건도 변화 제한 기반)
    # 한 Row = Nt개 튜브 직렬 → dx_max = DX_MAX_PER_TUBE × Nt
    dx_max_row = DX_MAX_PER_TUBE * Nt
    Q_max_2ph = dx_max_row * m_ref * h_fg  # 이상 영역 Row당 Q 상한

    # ── 공기측 h_o ──
    h_o = air_htc_ft(spec, geo, V_face) if spec.hx_type == 'FT' \
          else air_htc_mchx(spec, geo, V_face)

    eta_o = geo['eta_o']
    for _ in range(5):
        eta_o_new = _update_eta_o(spec, geo, h_o)
        if abs(eta_o_new - eta_o) < 0.001: break
        eta_o = eta_o_new

    # ── Row당 면적 ──
    A_o_row = geo['A_total'] / Nr
    A_i_row = geo['A_i'] / Nr
    R_wall_row = geo['R_wall'] * Nr

    # ── 습공기 함수 ──
    if _CP:
        def get_W(T, RH): return HAPropsSI('W', 'T', T+273.15, 'R', max(RH,0.01), 'P', P_atm)
        def get_h(T, RH): return HAPropsSI('H', 'T', T+273.15, 'R', max(RH,0.01), 'P', P_atm)
        def get_RH(T, W):
            try: return min(max(HAPropsSI('R', 'T', T+273.15, 'W', W, 'P', P_atm), 0.01), 1.0)
            except: return 0.5
        def get_Tdp(T, RH):
            try: return HAPropsSI('D', 'T', T+273.15, 'R', max(RH,0.01), 'P', P_atm) - 273.15
            except: return T - 10
        def get_Wsat(T): return HAPropsSI('W', 'T', T+273.15, 'R', 1.0, 'P', P_atm)
        def get_hsat(T): return HAPropsSI('H', 'T', T+273.15, 'R', 1.0, 'P', P_atm)
    else:
        def _ps(T):
            Tk=T+273.15
            return np.exp(-5.8002206e3/Tk+1.3914993-4.864e-2*Tk+4.1765e-5*Tk**2-1.4452e-8*Tk**3+6.546*np.log(Tk))
        def get_W(T, RH): Pw=max(RH,0.01)*_ps(T); return 0.62198*Pw/(P_atm-Pw)
        def get_h(T, RH): W=get_W(T,RH); return 1006*T+W*(2501000+1860*T)
        def get_RH(T, W): Pw=W*P_atm/(0.62198+W); return min(max(Pw/_ps(T),0.01),1.0)
        def get_Tdp(T, RH):
            Pw=max(RH,0.01)*_ps(T)
            return 243.5*np.log(Pw/611.2)/(17.67-np.log(Pw/611.2)) if Pw>0 else T-20
        def get_Wsat(T): return get_W(T, 1.0)
        def get_hsat(T): W=get_Wsat(T); return 1006*T+W*(2501000+1860*T)

    # ── 초기 상태 ──
    T_air = T_air_in
    RH_air = RH_in
    W_air = get_W(T_air, RH_air)
    x_ref = x_in
    T_ref = T_sat
    if x_ref > 1.0:
        T_ref = T_sat + (x_ref - 1.0) * h_fg / cp_v
    elif x_ref < 0:
        T_ref = T_sat + x_ref * h_fg / cp_l

    rows = []
    Q_total = 0; Q_sen_total = 0; Q_lat_total = 0

    for i_row in range(Nr):
        # ── ρ(T) 보정 ──
        rho_air = P_atm / (287.05 * (T_air + 273.15))
        m_air = rho_air * V_face * spec.W * spec.H

        # ── Q 추정 (Bo 안정화) ──
        T_wall_est = T_sat if x_ref <= 1.0 else T_ref
        dT_est = max(T_air - T_wall_est, 0.1) if side == 'evap' else max(T_wall_est - T_air, 0.1)
        Q_est = min(m_air * cp_air * dT_est * 0.3, m_ref * h_fg * 0.3)

        # ── h_i (건도 기반) ──
        h_i, phase = refrigerant_htc_auto(spec, geo, ref, side, x_ref, m_ref, Q_est)

        # ── T_wall ──
        if phase in ('two_phase', 'transition') and x_ref <= 1.0:
            T_wall = T_sat
        else:
            T_wall = T_ref

        # ── UA ──
        R_o = 1.0 / (eta_o * h_o * A_o_row + 1e-9)
        R_i = 1.0 / (h_i * A_i_row + 1e-9)
        UA_row = 1.0 / (R_o + R_wall_row + R_i)
        NTU = UA_row / (m_air * cp_air)
        eps = 1.0 - np.exp(-NTU)

        # ── 열량 계산 ──
        T_dp_air = get_Tdp(T_air, RH_air)
        is_wet = (T_dp_air > T_wall) and (side == 'evap')

        if is_wet:
            h_air_in = get_h(T_air, RH_air)
            h_sat_w = get_hsat(T_wall)
            W_sat_w = get_Wsat(T_wall)

            # b-factor (상한 적용)
            T_mid = (T_air + T_wall) / 2
            W_s1 = get_Wsat(T_mid)
            W_s2 = get_Wsat(T_mid + 0.5)
            dWs_dT = (W_s2 - W_s1) / 0.5
            b = 1.0 + h_fg_lat * dWs_dT / cp_air
            b = np.clip(b, 1.0, B_MAX)  # ★ 상한 제한

            h_o_wet = h_o * b
            A_fin_ratio = geo['A_fin'] / (geo['A_total'] + 1e-9)
            eta_o_wet = 1.0 - A_fin_ratio * (1 - eta_o)
            eta_o_wet = max(min(eta_o_wet, 1.0), 0.3)
            UA_o_wet = eta_o_wet * h_o_wet * A_o_row
            UA_wet = 1.0 / (1.0/UA_o_wet + R_wall_row + R_i)
            NTU_wet = UA_wet / (m_air * cp_air)
            eps_wet = 1.0 - np.exp(-NTU_wet)

            Q_row = max(eps_wet * m_air * (h_air_in - h_sat_w), 0.0)
            W_out = W_air - eps_wet * (W_air - W_sat_w)
            W_out = max(W_out, W_sat_w)
            T_out = T_air - eps_wet * (T_air - T_wall)
            Q_lat_row = m_air * (W_air - W_out) * h_fg_lat
            Q_sen_row = max(Q_row - Q_lat_row, 0.0)
        else:
            if side == 'evap':
                dT = max(T_air - T_wall, 0.0)
            else:
                dT = max(T_wall - T_air, 0.0)
            Q_row = eps * m_air * cp_air * dT
            if Q_row > 0:
                T_out = T_air - Q_row/(m_air*cp_air) if side=='evap' else \
                        T_air + Q_row/(m_air*cp_air)
            else:
                T_out = T_air
            T_out = np.clip(T_out, -40, 120)
            W_out = W_air
            W_sat_w = get_Wsat(T_wall) if T_wall > -40 else 0
            Q_sen_row = Q_row
            Q_lat_row = 0.0

        # ── ★ Q 상한 적용 (이상 영역: dx 제한) ──
        if phase in ('two_phase', 'transition') and x_ref <= 1.0:
            if Q_row > Q_max_2ph:
                ratio = Q_max_2ph / (Q_row + 1e-9)
                Q_row *= ratio
                Q_sen_row *= ratio
                Q_lat_row *= ratio
                # 공기 출구 보정
                if is_wet:
                    T_out = T_air - ratio * eps_wet * (T_air - T_wall)
                    W_out = W_air - ratio * eps_wet * (W_air - W_sat_w)
                    W_out = max(W_out, W_sat_w)
                else:
                    T_out = T_air - ratio * eps * (T_air - T_wall) if side=='evap' else \
                            T_air + ratio * eps * (T_wall - T_air)
                    T_out = np.clip(T_out, -40, 120)

        # ── 냉매 상태 업데이트 ──
        if phase in ('two_phase', 'transition') and x_ref <= 1.0:
            dx = Q_row / (m_ref * h_fg + 1e-9)
            x_new = x_ref + dx if side == 'evap' else x_ref - dx

            if side == 'evap' and x_new > 1.0 and x_ref < 1.0:
                frac_2ph = np.clip((1.0 - x_ref) / (dx + 1e-9), 0, 1)
                Q_sh = Q_row * (1 - frac_2ph)
                T_ref = T_sat + Q_sh / (m_ref * cp_v + 1e-9)
                x_ref = 1.0 + (T_ref - T_sat) * cp_v / h_fg
            elif side == 'cond' and x_new < 0 and x_ref > 0:
                frac_2ph = np.clip(x_ref / (dx + 1e-9), 0, 1)
                Q_sc = Q_row * (1 - frac_2ph)
                T_ref = T_sat - Q_sc / (m_ref * cp_l + 1e-9)
                x_ref = -(T_sat - T_ref) * cp_l / h_fg
            else:
                x_ref = x_new
                T_ref = T_sat

        elif phase == 'superheated':
            dT_ref = Q_row / (m_ref * cp_v + 1e-9)
            T_ref = T_ref + dT_ref if side == 'evap' else T_ref - dT_ref
            x_ref = 1.0 + (T_ref - T_sat) * cp_v / h_fg

        elif phase == 'subcooled':
            dT_ref = Q_row / (m_ref * cp_l + 1e-9)
            T_ref = T_ref + dT_ref if side == 'evap' else T_ref - dT_ref
            x_ref = -(T_sat - T_ref) * cp_l / h_fg

        # ── 공기 업데이트 ──
        T_air = T_out
        W_air = W_out
        RH_air = get_RH(T_air, W_air)

        Q_total += Q_row
        Q_sen_total += Q_sen_row
        Q_lat_total += Q_lat_row

        rows.append(dict(
            row=i_row, phase=phase, Q=Q_row, Q_sen=Q_sen_row, Q_lat=Q_lat_row,
            h_i=h_i, h_o=h_o, eta_o=eta_o, UA=UA_row,
            T_air_in=T_air_in if i_row==0 else rows[-1]['T_air_out'],
            T_air_out=T_out, T_wall=T_wall,
            x_in=x_in if i_row==0 else rows[i_row-1]['x_out'],
            x_out=x_ref, T_ref=T_ref, is_wet=is_wet,
        ))

    # ── 결과 ──
    SHR = Q_sen_total / (Q_total + 1e-9)
    n_2ph = sum(1 for r in rows if r['phase'] == 'two_phase')
    n_sh  = sum(1 for r in rows if r['phase'] == 'superheated')
    n_sc  = sum(1 for r in rows if r['phase'] == 'subcooled')
    n_tr  = sum(1 for r in rows if r['phase'] == 'transition')

    rho_in = P_atm / (287.05 * (T_air_in + 273.15))
    m_air_in = rho_in * V_face * spec.W * spec.H

    return dict(
        Q_total=Q_total, Q_sen=Q_sen_total, Q_lat=Q_lat_total, SHR=SHR,
        T_out=T_air, W_out=W_air, T_ref_out=T_ref, x_out=x_ref,
        m_dot=m_air_in, m_ref=m_ref,
        rows=rows, Nr=Nr,
        phase_summary=dict(two_phase=n_2ph, superheated=n_sh,
                          subcooled=n_sc, transition=n_tr),
        eps=Q_total/(m_air_in*cp_air*abs(T_air_in-T_sat)+1e-9) if abs(T_air_in-T_sat)>0.1 else 0,
        NTU=0, UA_eff=0, eta_o_eff=eta_o,
        T_dp=get_Tdp(T_air_in, RH_in),
        h_in=get_h(T_air_in, RH_in), h_out=get_h(T_air, RH_air),
        h_sat_w=get_hsat(T_sat),
        W_in=get_W(T_air_in, RH_in),
        is_wet=any(r['is_wet'] for r in rows),
        q_cond=m_air_in*(get_W(T_air_in,RH_in)-W_air)*1000 if Q_lat_total > 0 else 0,
        dehumid_rate=m_air_in*(get_W(T_air_in,RH_in)-W_air)*3600*1000 if Q_lat_total > 0 else 0,
        eps_dry=0, NTU_dry=0, UA_dry=0,
    )


# ═══════════════════════════════════════════════════════════════════
#  Level 2 코일 모델 — Tube-by-Tube × Segment (CoilDesigner 수준)
# ═══════════════════════════════════════════════════════════════════




def compute_coil_v3(spec, geo, ref, T_air_in, RH_in, V_face,
                    m_ref, x_in, side='evap', N_seg=10, flow='counter'):
    """Level 2 Tube-Segment coil model. flow='counter' or 'parallel'."""
    Nr = getattr(spec, 'tube_rows', getattr(spec, 'n_slabs', 1))
    Nt = getattr(spec, 'tube_cols', 1)
    Nr = max(Nr, 1); Nt = max(Nt, 1)
    N_tubes = Nr * Nt
    cp_air = 1006.0; P_atm = 101325.0; h_fg_lat = 2501000.0
    try:
        from CoolProp.HumidAirProp import HAPropsSI
        _CP = True
    except ImportError:
        _CP = False
    rf = ref.refrigerant
    if side == 'evap':
        T_sat = ref.T_sat_evap; T_sat_K = T_sat + 273.15; h_fg = ref.h_fg_evap
    else:
        T_sat = ref.T_sat_cond; T_sat_K = T_sat + 273.15; h_fg = ref.h_cond_v - ref.h_cond_l
    cp_v = CP.PropsSI('C', 'T', T_sat_K, 'Q', 1, rf)
    cp_l = CP.PropsSI('C', 'T', T_sat_K, 'Q', 0, rf)
    h_o = air_htc_ft(spec, geo, V_face) if spec.hx_type == 'FT' else air_htc_mchx(spec, geo, V_face)
    eta_o = geo['eta_o']
    for _ in range(5):
        eta_o_new = _update_eta_o(spec, geo, h_o)
        if abs(eta_o_new - eta_o) < 0.001: break
        eta_o = eta_o_new
    A_o_seg = geo['A_total'] / (N_tubes * N_seg)
    A_i_seg = geo['A_i'] / (N_tubes * N_seg)
    R_wall_seg = geo['R_wall'] * N_tubes * N_seg
    A_fin_ratio = geo['A_fin'] / (geo['A_total'] + 1e-9)
    if _CP:
        def get_W(T,RH): return HAPropsSI('W','T',T+273.15,'R',max(RH,0.01),'P',P_atm)
        def get_h(T,RH): return HAPropsSI('H','T',T+273.15,'R',max(RH,0.01),'P',P_atm)
        def get_RH(T,W):
            try: return min(max(HAPropsSI('R','T',T+273.15,'W',W,'P',P_atm),0.01),1.0)
            except: return 0.5
        def get_Tdp(T,RH):
            try: return HAPropsSI('D','T',T+273.15,'R',max(RH,0.01),'P',P_atm)-273.15
            except: return T-10
        def get_Wsat(T): return HAPropsSI('W','T',T+273.15,'R',1.0,'P',P_atm)
        def get_hsat(T): return HAPropsSI('H','T',T+273.15,'R',1.0,'P',P_atm)
    else:
        def _ps(T):
            Tk=T+273.15
            return np.exp(-5.8002206e3/Tk+1.3914993-4.864e-2*Tk+4.1765e-5*Tk**2-1.4452e-8*Tk**3+6.546*np.log(Tk))
        def get_W(T,RH): Pw=max(RH,0.01)*_ps(T); return 0.62198*Pw/(P_atm-Pw)
        def get_h(T,RH): W=get_W(T,RH); return 1006*T+W*(2501000+1860*T)
        def get_RH(T,W): Pw=W*P_atm/(0.62198+W); return min(max(Pw/_ps(T),0.01),1.0)
        def get_Tdp(T,RH):
            Pw=max(RH,0.01)*_ps(T)
            return 243.5*np.log(Pw/611.2)/(17.67-np.log(Pw/611.2)) if Pw>0 else T-20
        def get_Wsat(T): return get_W(T,1.0)
        def get_hsat(T): W=get_Wsat(T); return 1006*T+W*(2501000+1860*T)

    def calc_seg(T_air, W_air, RH_air, x_ref, T_ref, m_air):
        """단일 세그먼트 — T_wall 반복 수렴 (Newton-Raphson식)

        에너지 보존: Q_air(T_wall) = Q_ref = (T_wall - T_ref_base) × h_i × A_i
        → T_wall을 찾아 Q_air = Q_ref 동시 만족
        """
        Q_est = min(m_air*cp_air*max(abs(T_air-T_sat),0.1)*0.1, m_ref*h_fg*0.05)
        h_i, phase = refrigerant_htc_auto(spec, geo, ref, side, x_ref, m_ref, max(Q_est,1.0))

        # 냉매 기준 온도 (T_wall의 하한)
        if phase in ('two_phase','transition') and x_ref <= 1.0:
            T_ref_base = T_sat
        else:
            T_ref_base = T_ref

        R_i = 1.0 / (h_i * A_i_seg + 1e-9)

        # ── 공기측 Q 계산 함수 (T_wall 의존) ──
        def Q_air_func(Tw):
            """주어진 T_wall에서 공기측 Q 계산"""
            R_o_loc = 1.0/(eta_o*h_o*A_o_seg+1e-9)
            UA_loc = 1.0/(R_o_loc + R_wall_seg + R_i)
            NTU_loc = UA_loc/(m_air*cp_air+1e-9)
            eps_loc = 1.0-np.exp(-min(NTU_loc,7.0))

            T_dp_loc = get_Tdp(T_air, RH_air)
            wet = (T_dp_loc > Tw) and (side == 'evap')

            if wet:
                h_a = get_h(T_air, RH_air)
                h_sw = get_hsat(Tw); W_sw = get_Wsat(Tw)
                T_mid = (T_air + Tw) / 2
                dWs = (get_Wsat(T_mid+0.5) - get_Wsat(T_mid)) / 0.5
                b_loc = max(1.0 + h_fg_lat*dWs/cp_air, 1.0)
                k_f = 230.0 if getattr(spec,'fin_material','Al')=='Al' else 386.0
                df = getattr(spec,'fin_thickness',0.1e-3)
                m_w = np.sqrt(2*h_o*b_loc/(k_f*df+1e-9))
                if hasattr(spec,'tube_do'):
                    ri=spec.tube_do/2; Xm=spec.tube_pitch_t/2
                    XL=np.sqrt((spec.tube_pitch_t/2)**2+spec.tube_pitch_l**2)/2
                    r_eq=max(1.27*Xm*np.sqrt(max(XL/Xm-0.3,0.01)),ri*1.01)
                    phi=(r_eq/ri-1)*(1+0.35*np.log(r_eq/ri)); mL=m_w*ri*phi
                else: mL=m_w*getattr(spec,'slab_pitch',8e-3)/2
                eta_fw = np.tanh(mL)/(mL+1e-9)
                eta_ow = max(1.0-A_fin_ratio*(1-eta_fw),0.2)
                UA_w = 1.0/(1.0/(eta_ow*h_o*A_o_seg)+R_wall_seg+R_i)
                eps_w = 1.0-np.exp(-min(UA_w/(m_air*cp_air+1e-9),7.0))
                Q_a = max(eps_w*m_air*(h_a - h_sw), 0.0)
                W_o = max(W_air - eps_w*(W_air - W_sw), W_sw)
                T_o = T_air - eps_w*(T_air - Tw)
            else:
                if side == 'evap':
                    dT = max(T_air - Tw, 0.0)
                else:
                    dT = max(Tw - T_air, 0.0)
                Q_a = eps_loc*m_air*cp_air*dT
                T_o = T_air-Q_a/(m_air*cp_air) if side=='evap' and Q_a>0 else \
                      T_air+Q_a/(m_air*cp_air) if side=='cond' and Q_a>0 else T_air
                T_o = np.clip(T_o, -40, 120)
                W_o = W_air; W_sw = get_Wsat(Tw) if Tw > -40 else 0
            return Q_a, T_o, W_o, wet

        # ── T_wall 반복 수렴 (successive substitution + relaxation) ──
        # 원리: Q_air(T_wall) = Q_ref = (T_wall - T_ref_base) / R_i
        #   → T_wall = T_ref_base + Q_air(T_wall) × R_i
        T_wall = (T_ref_base + T_air) / 2  # 초기값: 중간점
        Q_converged = 0.0
        T_o_conv = T_air; W_o_conv = W_air; is_wet_conv = False
        alpha = 0.3  # 완화 계수 (0.3 = 보수적 수렴)

        for _iter in range(20):
            Q_a, T_o_a, W_o_a, wet_a = Q_air_func(T_wall)

            # 목표 T_wall: T_ref_base + Q_air × R_i
            T_wall_target = T_ref_base + Q_a * R_i
            # T_wall 범위 제한
            T_wall_target = np.clip(T_wall_target, T_ref_base,
                                    T_air if side == 'evap' else T_ref_base + 80)

            # 완화 업데이트
            T_wall_new = alpha * T_wall_target + (1 - alpha) * T_wall

            if abs(T_wall_new - T_wall) < 0.05:  # 0.05°C 이내 수렴
                Q_converged = Q_a
                T_o_conv = T_o_a; W_o_conv = W_o_a; is_wet_conv = wet_a
                T_wall = T_wall_new
                break

            T_wall = T_wall_new
            Q_converged = Q_a
            T_o_conv = T_o_a; W_o_conv = W_o_a; is_wet_conv = wet_a

        Q = Q_converged

        # ── 냉매 상태 업데이트 ──
        x_n = x_ref; T_n = T_ref
        if phase in ('two_phase','transition') and x_ref <= 1.0:
            dx = Q/(m_ref*h_fg+1e-9)
            if side == 'evap':
                x_n = x_ref + dx
                if x_n > 1.0 and x_ref < 1.0:
                    frac = np.clip((1.0-x_ref)/(dx+1e-9), 0, 1)
                    T_n = T_sat + Q*(1-frac)/(m_ref*cp_v+1e-9)
                    x_n = 1.0 + (T_n-T_sat)*cp_v/h_fg
                else: T_n = T_sat
            else:
                x_n = x_ref - dx
                if x_n < 0 and x_ref > 0:
                    frac = np.clip(x_ref/(dx+1e-9), 0, 1)
                    T_n = T_sat - Q*(1-frac)/(m_ref*cp_l+1e-9)
                    x_n = -(T_sat-T_n)*cp_l/h_fg
                else: T_n = T_sat
        elif phase == 'superheated':
            dTr = Q/(m_ref*cp_v+1e-9)
            T_n = T_ref+dTr if side=='evap' else T_ref-dTr
            x_n = 1.0+(T_n-T_sat)*cp_v/h_fg
        elif phase == 'subcooled':
            dTr = Q/(m_ref*cp_l+1e-9)
            T_n = T_ref+dTr if side=='evap' else T_ref-dTr
            x_n = -(T_sat-T_n)*cp_l/h_fg

        return Q, T_o_conv, W_o_conv, get_RH(T_o_conv, W_o_conv), x_n, T_n, phase, h_i, T_wall, is_wet_conv

    def build_tube_order(row_seq):
        order = []
        for i, row in enumerate(row_seq):
            tubes = list(range(row*Nt, (row+1)*Nt))
            if i % 2 == 1: tubes.reverse()
            order.extend(tubes)
        return order

    if flow == 'counter':
        ref_row_seq = list(range(Nr-1, -1, -1))
    else:
        ref_row_seq = list(range(Nr))
    tube_order = build_tube_order(ref_row_seq)

    T_air_rows = [T_air_in] * Nr
    W_air_rows = [get_W(T_air_in, RH_in)] * Nr
    RH_air_rows = [RH_in] * Nr
    N_iter = 5  # 공기↔냉매 수렴 반복 (counter/parallel 모두)
    final = None

    for iteration in range(N_iter):
        x_ref = x_in; T_ref = T_sat
        if x_ref > 1.0: T_ref = T_sat + (x_ref-1.0)*h_fg/cp_v
        elif x_ref < 0: T_ref = T_sat + x_ref*h_fg/cp_l
        all_segs = []; tube_results = []
        Q_per_row = [0.0]*Nr; Q_lat_per_row = [0.0]*Nr

        for tube_id in tube_order:
            row = tube_id // Nt
            rho_a = P_atm/(287.05*(T_air_rows[row]+273.15))
            m_air = rho_a*V_face*spec.W*spec.H
            m_air_seg = m_air / (Nt * N_seg)
            Q_t=0; Q_st=0; Q_lt=0; segs=[]
            for seg in range(N_seg):
                Q_s,T_o,W_o,RH_o,x_n,T_n,phase,h_i,T_w,wet = \
                    calc_seg(T_air_rows[row],W_air_rows[row],RH_air_rows[row],x_ref,T_ref,m_air_seg)
                Q_lat_s = m_air_seg*max(W_air_rows[row]-W_o,0)*h_fg_lat if wet else 0.0
                Q_lat_s = min(Q_lat_s, Q_s); Q_sen_s = max(Q_s-Q_lat_s, 0.0)
                segs.append(dict(tube=tube_id,seg=seg,row=row,x_in=x_ref,x_out=x_n,
                    T_ref=T_n,T_wall=T_w,Q=Q_s,Q_sen=Q_sen_s,Q_lat=Q_lat_s,
                    h_i=h_i,phase=phase,is_wet=wet))
                Q_t+=Q_s; Q_st+=Q_sen_s; Q_lt+=Q_lat_s; x_ref=x_n; T_ref=T_n
            all_segs.extend(segs)
            tube_results.append(dict(tube=tube_id,row=row,Q=Q_t,Q_sen=Q_st,Q_lat=Q_lt,
                x_in=segs[0]['x_in'],x_out=segs[-1]['x_out'],T_ref_out=segs[-1]['T_ref'],
                h_i_avg=np.mean([s['h_i'] for s in segs]),
                phase_dominant=max(set(s['phase'] for s in segs),
                    key=lambda p: sum(1 for s in segs if s['phase']==p)),
                segments=segs))
            Q_per_row[row] += Q_t; Q_lat_per_row[row] += Q_lt

        T_a = T_air_in; W_a = get_W(T_air_in, RH_in)
        for row in range(Nr):
            T_air_rows[row] = T_a; W_air_rows[row] = W_a
            RH_air_rows[row] = get_RH(T_a, W_a)
            if Q_per_row[row] > 0:
                rho_a = P_atm/(287.05*(T_a+273.15))
                m_a = rho_a*V_face*spec.W*spec.H
                T_a -= Q_per_row[row]/(m_a*cp_air+1e-9)
                T_a = np.clip(T_a, -40, 120)
                if Q_lat_per_row[row] > 0:
                    dW = Q_lat_per_row[row]/(m_a*h_fg_lat+1e-9)
                    W_a = max(W_a-dW, get_Wsat(T_sat) if side=='evap' else 0)
                # ★ 포화 보정: W > Wsat(T)이면 과포화 → 응결
                W_sat_T = get_Wsat(T_a)
                if W_a > W_sat_T:
                    W_a = W_sat_T
        final = (all_segs, tube_results, Q_per_row, Q_lat_per_row, x_ref, T_ref, T_a, W_a)

    all_segs, tube_results, Q_per_row, _, x_ref, T_ref, T_air_out, W_air_out = final
    Q_total = sum(Q_per_row)
    Q_sen_total = sum(s['Q_sen'] for s in all_segs)
    Q_lat_total = sum(s['Q_lat'] for s in all_segs)
    SHR = Q_sen_total/(Q_total+1e-9)
    phases = [s['phase'] for s in all_segs]
    rho_in = P_atm/(287.05*(T_air_in+273.15))
    m_air_in = rho_in*V_face*spec.W*spec.H
    return dict(
        Q_total=Q_total, Q_sen=Q_sen_total, Q_lat=Q_lat_total, SHR=SHR,
        T_out=T_air_out, W_out=W_air_out, T_ref_out=T_ref, x_out=x_ref,
        m_dot=m_air_in, m_ref=m_ref, flow=flow,
        tubes=tube_results, Nr=Nr, Nt=Nt, N_seg=N_seg,
        N_total_seg=len(all_segs), Q_per_row=Q_per_row,
        phase_summary=dict(two_phase=phases.count('two_phase'),
            superheated=phases.count('superheated'),
            subcooled=phases.count('subcooled'),
            transition=phases.count('transition')),
        eps=Q_total/(m_air_in*cp_air*abs(T_air_in-T_sat)+1e-9) if abs(T_air_in-T_sat)>0.1 else 0,
        T_dp=get_Tdp(T_air_in,RH_in), h_in=get_h(T_air_in,RH_in),
        W_in=get_W(T_air_in,RH_in), is_wet=any(s['is_wet'] for s in all_segs),
        rows=[], eps_dry=0, NTU=0, NTU_dry=0, UA_eff=0, UA_dry=0, eta_o_eff=eta_o,
        h_sat_w=get_hsat(T_sat), h_out=get_h(T_air_out,get_RH(T_air_out,W_air_out)),
        q_cond=m_air_in*(get_W(T_air_in,RH_in)-W_air_out)*1000 if Q_lat_total>0 else 0,
        dehumid_rate=m_air_in*(get_W(T_air_in,RH_in)-W_air_out)*3600*1000 if Q_lat_total>0 else 0,
    )


# ─── 시각화 스타일 ─────────────────────────────────────
import os, matplotlib.pyplot as plt
from matplotlib import font_manager as fm

DARK = dict(bg="#ffffff", panel="#f5f5f8", border="#c8ccd4",
            text="#1a1a2e", dim="#6a7080", grid="#e0e2e8")
C    = ["#0077cc", "#d63031", "#e67e22", "#00875a", "#8e44ad", "#e17055"]
RISK_CLR = {"안전": "#2a9d5c", "주의": "#ffc940", "위험": "#ff7c43", "매우위험": "#e03030"}
FLUX_CAUTION = 1.0; FLUX_DANGER = 4.0; FLUX_SEVERE = 12.0

_NOTO = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
if os.path.exists(_NOTO):
    plt.rcParams["font.family"] = "Noto Sans CJK JP"
plt.rcParams["axes.unicode_minus"] = False

def apply_style():
    plt.rcParams.update({
        "figure.facecolor": DARK["bg"], "axes.facecolor": DARK["panel"],
        "axes.edgecolor": DARK["border"], "axes.labelcolor": DARK["text"],
        "text.color": DARK["text"], "xtick.color": DARK["dim"],
        "ytick.color": DARK["dim"], "grid.color": DARK["grid"],
        "grid.linestyle": "--", "grid.linewidth": 0.5,
        "legend.facecolor": DARK["panel"], "legend.edgecolor": DARK["border"],
        "axes.titlepad": 8,
    })


# ═══════════════════════════════════════════════════════════════════
#  Level 1 코일 모델 — Row × Phase-zone 세그먼트
# ═══════════════════════════════════════════════════════════════════

def _ref_vapor_props(ref, T_vapor_K):
    """과열 증기 물성 (CoolProp)"""
    rf = ref.refrigerant
    P = ref.P_evap
    try:
        mu = CP.PropsSI('V', 'T', T_vapor_K, 'P', P, rf)
        k  = CP.PropsSI('L', 'T', T_vapor_K, 'P', P, rf)
        Pr = CP.PropsSI('Prandtl', 'T', T_vapor_K, 'P', P, rf)
        cp = CP.PropsSI('C', 'T', T_vapor_K, 'P', P, rf)
        rho = CP.PropsSI('D', 'T', T_vapor_K, 'P', P, rf)
    except Exception:
        # fallback: 포화 증기 물성
        Ts = ref.T_sat_evap + 273.15
        mu = CP.PropsSI('V', 'T', Ts, 'Q', 1, rf)
        k  = CP.PropsSI('L', 'T', Ts, 'Q', 1, rf)
        Pr = CP.PropsSI('Prandtl', 'T', Ts, 'Q', 1, rf)
        cp = CP.PropsSI('C', 'T', Ts, 'Q', 1, rf)
        rho = CP.PropsSI('D', 'T', Ts, 'Q', 1, rf)
    return dict(mu=mu, k=k, Pr=Pr, cp=cp, rho=rho)


def _ref_liquid_props(ref, T_liquid_K):
    """과냉 액체 물성 (CoolProp)"""
    rf = ref.refrigerant
    P = ref.P_cond
    try:
        mu = CP.PropsSI('V', 'T', T_liquid_K, 'P', P, rf)
        k  = CP.PropsSI('L', 'T', T_liquid_K, 'P', P, rf)
        Pr = CP.PropsSI('Prandtl', 'T', T_liquid_K, 'P', P, rf)
        cp = CP.PropsSI('C', 'T', T_liquid_K, 'P', P, rf)
        rho = CP.PropsSI('D', 'T', T_liquid_K, 'P', P, rf)
    except Exception:
        Ts = ref.T_sat_cond + 273.15
        mu = CP.PropsSI('V', 'T', Ts, 'Q', 0, rf)
        k  = CP.PropsSI('L', 'T', Ts, 'Q', 0, rf)
        Pr = CP.PropsSI('Prandtl', 'T', Ts, 'Q', 0, rf)
        cp = CP.PropsSI('C', 'T', Ts, 'Q', 0, rf)
        rho = CP.PropsSI('D', 'T', Ts, 'Q', 0, rf)
    return dict(mu=mu, k=k, Pr=Pr, cp=cp, rho=rho)


def refrigerant_htc_v2(x, D_ch, G_ref, ref, side='evap', T_ref_K=None):
    """
    건도 기반 냉매 HTC — 상 영역 자동 분기 + 전이 블렌딩

    x < 0:          과냉 → Gnielinski (액체)
    0 ≤ x ≤ 0.95:   이상 → Shah (1982/1979)
    0.95 < x < 1.05: 전이 → 블렌딩 (이상↔단상)
    x ≥ 1.05:       과열 → Gnielinski (증기)

    Returns: h_i [W/m²K], phase ('subcooled'/'two_phase'/'transition'/'superheated')
    """
    X_BLEND_LO = 0.95
    X_BLEND_HI = 1.05

    if side == 'evap':
        T_sat = ref.T_sat_evap
        P = ref.P_evap
        h_fg = ref.h_fg_evap
    else:
        T_sat = ref.T_sat_cond
        P = ref.P_cond
        h_fg = ref.h_cond_v - ref.h_cond_l

    T_sat_K = T_sat + 273.15

    # ── 과냉 (x < 0) ──
    if x < 0:
        T_liq_K = T_ref_K if T_ref_K else T_sat_K - 5.0
        props = _ref_liquid_props(ref, min(T_liq_K, T_sat_K - 0.5))
        Re = G_ref * D_ch / props['mu']
        h_i = gnielinski_htc(Re, props['Pr'], props['k'], D_ch)
        return h_i, 'subcooled'

    # ── 완전 과열 (x ≥ 1.05) ──
    if x >= X_BLEND_HI:
        if T_ref_K is None:
            T_ref_K = T_sat_K + 10.0
        T_ref_K = max(T_ref_K, T_sat_K + 1.0)
        props = _ref_vapor_props(ref, T_ref_K)
        Re = G_ref * D_ch / props['mu']
        h_i = gnielinski_htc(Re, props['Pr'], props['k'], D_ch)
        return h_i, 'superheated'

    # ── 이상 (0 ≤ x ≤ 0.95) ──
    x_2ph = min(x, 0.95)
    x_2ph = max(x_2ph, 0.05)  # Shah 안정 범위

    if side == 'evap':
        # Shah (1982)
        Re_l = max(G_ref * (1 - x_2ph) * D_ch / ref.mu_l_evap, 2300)
        Nu_l = 0.023 * Re_l**0.8 * ref.Pr_l_evap**0.4
        h_l = Nu_l * ref.k_l_evap / D_ch
        Co = ((1 - x_2ph) / x_2ph)**0.8 * (ref.rho_v_evap / ref.rho_l_evap)**0.5
        # Bo 추정: q''/(G×h_fg), q'' ≈ h_o_typ × ΔT_typ ≈ 100×20 = 2000 W/m²
        Bo = max(2000.0 / (G_ref * h_fg + 1e-3), 1e-6)
        if Co > 0.65:
            psi = 1.8 / Co**0.8
        else:
            psi = max(1.8 / Co**0.8, 0.6683 * Co**(-0.2) + 1058 * Bo**0.7)
        h_2ph = h_l * psi
    else:
        # Shah (1979)
        Tc = T_sat_K
        mu_l = CP.PropsSI('V', 'T', Tc, 'Q', 0, ref.refrigerant)
        k_l = CP.PropsSI('L', 'T', Tc, 'Q', 0, ref.refrigerant)
        Pr_l = CP.PropsSI('Prandtl', 'T', Tc, 'Q', 0, ref.refrigerant)
        Re_l = G_ref * D_ch / mu_l
        h_l = 0.023 * Re_l**0.8 * Pr_l**0.4 * k_l / D_ch
        P_crit = CP.PropsSI('Pcrit', ref.refrigerant)
        Z = (1 / x_2ph - 1)**0.8 * (P / P_crit)**0.4
        h_2ph = h_l * (1 + 3.8 / (Z**0.95 + 1e-9))

    if x <= X_BLEND_LO:
        return h_2ph, 'two_phase'

    # ── 전이 구간 (0.95 < x < 1.05) → 블렌딩 ──
    T_sh_K = T_sat_K + 2.0  # 전이 구간 약간의 과열
    props = _ref_vapor_props(ref, T_sh_K)
    Re_v = G_ref * D_ch / props['mu']
    h_vapor = gnielinski_htc(Re_v, props['Pr'], props['k'], D_ch)

    w = (x - X_BLEND_LO) / (X_BLEND_HI - X_BLEND_LO)
    w = max(0, min(1, w))
    h_i = (1 - w) * h_2ph + w * h_vapor
    return h_i, 'transition'

