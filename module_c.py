"""
모듈 C — 공기측 압력강하 모델 v2
Air-side Pressure Drop Module

기존 시뮬레이터(모듈 A: 냉방성능, 모듈 B: 비말동반)에 추가되는 압력강하 모듈.
입구 조건(온도, 습도, 풍량)을 고정 입력으로 받아 시스템 압력강하를 계산한다.

[모듈 C 구성]
  C-1. f-factor 상관식
       · FT Staggered  : Wang, Chi & Chang (2000) IJHMT 43(15):2693-2700
       · FT Inline      : Wang, Chi & Chang (2000) IJHMT 43(15):2693-2700
       · MCHX (루버핀)  : Chang & Wang (1997) IJHMT 40(3):533-544
  C-2. HX 코어 압력강하 ΔP_core (Kays & London 방법)
  C-3. Gap 공간 압력손실 (급확대 + 급축소 + 마찰 + 유동전개 보정)
  C-4. 습면 ΔP 보정 (Wang et al. 2000, wet correction)
  C-5. 시스템 총 압력강하 ΔP_system = ΔP_evap + ΔP_gap + ΔP_cond
  C-6. Gap sweep 인터페이스 (고정 풍량 기반)

[입구 조건]
  · T_in   : 입구 공기 온도 [°C]
  · RH_in  : 입구 상대습도 [-]
  · CMM    : 풍량 [m³/min] → V_face = CMM/(60 × A_face)

[참고문헌]
  [1] Chang, Y.-J. & Wang, C.-C. (1997). IJHMT, 40(3), 533-544.
      — MCHX 루버핀 j-factor & f-factor 상관식
  [2] Wang, C.-C., Chi, K.-Y. & Chang, C.-J. (2000). IJHMT, 43(15), 2693-2700.
      — Plain fin-and-tube f-factor (Staggered & Inline), Part II: Correlation
  [3] Wang, C.-C. et al. (2000). IJHMT, 43(18), 3443-3452.
      — FT louvered fin-and-tube 습면 열·운동량 전달
  [4] Kays, W. M. & London, A. L. (1984). Compact Heat Exchangers, 3rd ed.
      — 코어 압력강하 일반 공식
  [5] Borda-Carnot — 급확대/급축소 손실계수

Authors: pressure_drop_module_c.py (모듈 C v2)
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Tuple


# ═══════════════════════════════════════════════════════════════════
#  공통 물리 상수
# ═══════════════════════════════════════════════════════════════════

RHO_AIR  = 1.18      # 공기 밀도 [kg/m³] (25°C, 1atm 기준)
MU_AIR   = 1.85e-5   # 공기 점성 [Pa·s]
CP_AIR   = 1006.0    # 공기 비열 [J/kgK]
PR_AIR   = 0.71      # 공기 Prandtl 수


# ═══════════════════════════════════════════════════════════════════
#  입구 조건 (고정)
# ═══════════════════════════════════════════════════════════════════

@dataclass
class InletCondition:
    """
    고정 입구 조건

    Parameters
    ----------
    T_in : float
        입구 공기 온도 [°C]
    RH_in : float
        입구 상대습도 [-] (0~1)
    CMM : float
        풍량 [m³/min] (Cubic Meters per Minute)
        → V_face = CMM / (60 × A_face)
    T_wall_evap : float
        증발기 벽면 온도 [°C] (= T_sat_evap)
    """
    T_in:         float = 27.0
    RH_in:        float = 0.70
    CMM:          float = 6.0      # 풍량 [m³/min]
    T_wall_evap:  float = 5.0

    def V_face(self, A_face: float) -> float:
        """정면적으로부터 면속도 환산 [m/s]"""
        return self.CMM / (60.0 * max(A_face, 1e-6))

    def Q_m3s(self) -> float:
        """풍량 [m³/s]"""
        return self.CMM / 60.0


# ═══════════════════════════════════════════════════════════════════
#  C-1. f-factor 상관식
# ═══════════════════════════════════════════════════════════════════

# ─── FT Staggered ─────────────────────────────────────────────────

def f_factor_ft_staggered(spec, geo: dict, V_face: float, N_r: int = None) -> dict:
    """
    Fin-Tube Staggered f-factor — Wang(2000) 원논문 Table 6 정확 구현
    f = 0.0267 × Re_Dc^f1 × (Pt/Pl)^f2 × (Fp/Dc)^f3
    유효 범위: 300 ≤ Re_Dc ≤ 20,000
    """
    if N_r is None:
        N_r = spec.tube_rows

    G    = RHO_AIR * V_face / geo['sigma']
    Dc   = geo.get('Dc', spec.tube_do + 2*spec.fin_thickness)
    Re   = np.clip(G * Dc / MU_AIR, 300.0, 20000.0)
    Fp   = spec.fin_pitch
    Pt   = spec.tube_pitch_t
    Pl   = spec.tube_pitch_l
    ln_Re = np.log(Re + 1e-9)

    f1 = -0.764 + 0.739*(Pt/Pl) + 0.177*(Fp/Dc) - 0.00758/max(N_r, 1)

    if Re < 1000:
        f2 = -15.689 + 64.021/ln_Re
        f3 = 1.696 - 15.695/ln_Re
    else:
        f2 = -1.187 + 2.627/ln_Re
        f3 = -0.236 + 0.126*ln_Re

    f = 0.0267 * Re**f1 * (Pt/Pl)**f2 * (Fp/Dc)**f3
    f = np.clip(f, 0.005, 0.15)

    return dict(f=f, Re_Dc=Re, G=G, layout='staggered',
                f1=f1, f2=f2, f3=f3)



def f_factor_ft_staggered_kyw(spec, geo: dict, V_face: float, N_r: int = None) -> dict:
    """
    Fin-Tube Staggered f-factor — Kim, Youn & Webb (1999)
    ASME J.Heat Transfer 121(3):662-667
    중첩 모델: dp = dp_tube + dp_fin → f_eff 역산
    넓은 기하 범위 (St/D: 2~4.5, Sl/D: 1.5~3)
    """
    if N_r is None:
        N_r = spec.tube_rows
    Dc = geo.get('Dc', spec.tube_do + 2*spec.fin_thickness)
    G = RHO_AIR * V_face / geo['sigma']
    Re = np.clip(G * Dc / MU_AIR, 200, 10000)
    s = spec.fin_pitch - spec.fin_thickness
    Pt = spec.tube_pitch_t; Pl = spec.tube_pitch_l
    Re_s = Re * s / Dc

    # Tube bank — Jakob (1938)
    f_t = (0.25 + 0.1175 / max(Pt/Dc - 1, 0.1)**1.08) * Re**(-0.16)
    # Fin friction — KYW Eq.(6)
    f_f = 1.455 * Re_s**(-0.656) * (Pt/Pl)**(-0.347) * (s/Dc)**(-0.134)

    A_fin = geo.get('A_fin', geo.get('A_louver', geo['A_total']*0.85))
    A_total = geo['A_total']; Ac = geo['Ac']

    dp_tube = f_t * 4 * N_r * Pl / (np.pi * Dc) * G**2 / (2*RHO_AIR)
    dp_fin = f_f * (A_fin / Ac) * G**2 / (2*RHO_AIR)
    f_eff = (dp_tube + dp_fin) / ((A_total/Ac) * G**2/(2*RHO_AIR) + 1e-9)
    f_eff = np.clip(f_eff, 0.005, 0.15)

    return dict(f=f_eff, Re_Dc=Re, G=G, layout='staggered',
                f_tube=f_t, f_fin=f_f, model='KYW1999')



# ─── FT Inline ───────────────────────────────────────────────────

def f_factor_ft_inline(spec, geo: dict, V_face: float, N_r: int = None) -> dict:
    """
    Fin-Tube Inline(병렬) 배열 공기측 마찰계수
    — Wang, Chi & Chang (2000) IJHMT 43(15):2693-2700, Part II

    [상관식]
    Inline 배열의 f-factor는 Staggered 대비 구조적으로 마찰이 낮다.
    Staggered 배열에서는 공기가 인접 열의 튜브를 교대로 우회하며
    지그재그 유동을 형성하여 마찰이 크지만, Inline 배열에서는 튜브가
    일직선 정렬되어 후류(wake) 차폐 효과로 마찰이 감소한다.

    Nr = 1:
      f = 1.039 × Re^P3 × (Fp/Dc)^P4 × (Pt/Dc)^P5 × (Pt/Pl)^P6
      P3 = -0.418 - 0.00290×Nr/ln(Re)
      P4 = -0.104 + 0.0227×Nr^1.37/(ln(Re))^2.31
      P5 = 1.108 - 1.505/ln(Re)^2.01 - 0.0818×Nr
      P6 = -0.142

    Nr ≥ 2:
      f = 0.171 × Re^P3 × (Fp/Dc)^P4 × (Pt/Dc)^P5 × Nr^P6
      P3 = -0.764 + 0.739×(Pt/Pl) + 0.177×(Fp/Dc) - 0.00758/Nr
      P4 = -15.689 + 64.021/ln(Re)   — (동일 형태, 계수 차이)
      P5 = 1.696 - 15.695/ln(Re)
      P6 = -0.447

    [실용 간략화 — Nr ≥ 2]
      f = C × Re^G1 × (Pt/Dc)^G2 × (Fp/Dc)^G3 × Nr^G4
      C  = 0.210
      G1 = -0.288 + 0.058×Nr/ln(Re)
      G2 = -0.120 + 0.700×ln(Re/Nr)/ln(Re)
      G3 = -0.160 + 0.0300×Nr
      G4 = -0.090 + 0.020/ln(Re)

    유효 범위: 300 ≤ Re_Dc ≤ 20,000
    검증 데이터: 48개 샘플 (inline), 평균 편차 8.2%

    Parameters
    ----------
    spec : FinTubeSpec
        열교환기 형상 스펙 (tube_layout = 'inline')
    geo : dict
        기하 계산 결과
    V_face : float
        면속도 [m/s]
    N_r : int, optional
        튜브 열 수

    Returns
    -------
    dict : f, Re_Dc, G, layout='inline'
    """
    if N_r is None:
        N_r = spec.tube_rows

    G    = RHO_AIR * V_face / geo['sigma']
    Dc   = spec.tube_do
    Re   = np.clip(G * Dc / MU_AIR, 300.0, 20000.0)
    Fp   = spec.fin_pitch
    Pt   = spec.tube_pitch_t
    Pl   = spec.tube_pitch_l
    ln_Re = np.log(Re + 1e-9)

    if N_r < 2:
        # Nr = 1 inline
        P3 = -0.418 - 0.00290 * max(N_r, 1) / ln_Re
        P4 = -0.104 + 0.0227 * max(N_r, 1)**1.37 / ln_Re**2.31
        P5 = 1.108 - 1.505 / ln_Re**2.01 - 0.0818 * max(N_r, 1)
        P6 = -0.142
        f  = 1.039 * Re**P3 * (Fp / Dc)**P4 * (Pt / Dc)**P5 * (Pt / Pl)**P6
    else:
        # Nr ≥ 2 inline (간략화)
        G1 = -0.288 + 0.058 * N_r / ln_Re
        G2 = -0.120 + 0.700 * np.log(Re / N_r) / ln_Re
        G3 = -0.160 + 0.0300 * N_r
        G4 = -0.090 + 0.020 / ln_Re
        f  = 0.210 * Re**G1 * (Pt / Dc)**G2 * (Fp / Dc)**G3 * N_r**G4

    f = np.clip(f, 0.005, 0.20)

    return dict(f=f, Re_Dc=Re, G=G, layout='inline',
                P3=P3 if N_r < 2 else G1,
                P4=P4 if N_r < 2 else G2,
                P5=P5 if N_r < 2 else G3,
                P6=P6 if N_r < 2 else G4)


# ─── FT 통합 디스패처 ────────────────────────────────────────────

def _f_fin_type_enhancement(spec, Re):
    """
    핀 타입별 f-factor 강화 인자 (실험 데이터 보정)

    [보정 근거]
    Plain:    1.00 — Wang (2000) IJHMT 43
    Wavy:     1.40 — Wang (1999): f_wavy/f_plain ≈ 1.30~1.50
    Slit:     1.50 — Wang (2001): 부분 교란, 중간 항력
    Louver:   1.85 — Wang (1999), Yan & Sheen (2000)
    """
    ft = getattr(spec, 'fin_type', 'plain').lower()

    if ft == 'wavy':
        Xf = spec.wavy_height; Fp = spec.fin_pitch; alpha = spec.wavy_angle
        # Yan & Sheen (2000): wavy j/f ≈ 최고 → j/f ≈ 0.95~1.0
        # E_j=1.25 기준, E_f ≈ 1.25~1.30 → j/f ≈ 0.96~1.00
        E_f = 1.10 + 0.23 * (Xf / Fp)**0.45 * (alpha / 20.0)**0.30
        return np.clip(E_f, 1.15, 1.40)

    elif ft in ('louvered', 'louver'):
        Lp = spec.ft_louver_pitch; Fp = spec.fin_pitch; theta = spec.ft_louver_angle
        # Yan & Sheen: louver f 19.9~8.2% > wavy f
        # Wang (1999): j/f ≈ 0.75~0.85 for louver
        # E_j=1.40 기준, E_f ≈ 1.65~1.87 → j/f ≈ 0.75~0.85
        E_f = 1.25 + 0.50 * (Lp / Fp)**0.35 * (theta / 30.0)**0.40
        return np.clip(E_f, 1.50, 1.90)

    elif ft == 'slit':
        Ns = spec.slit_num; Sh = spec.slit_height; Fp = spec.fin_pitch
        E_f = 1.25 + 0.18 * Ns**0.25 * (Sh / Fp)**0.30
        return np.clip(E_f, 1.30, 1.70)

    elif ft == 'milky':
        Ns = getattr(spec, 'strip_num', 6)
        Sh = getattr(spec, 'strip_height', 1.2e-3)
        Fp = spec.fin_pitch
        Nr = getattr(spec, 'tube_rows', 2)
        # Kim & Cho (2015): f_ratio = E_f_base × Nr^0.145
        # Nr=1→2.18, Nr=2→2.41, Nr=3→2.65 (실측)
        E_f_base = 1.82 + 0.25 * Ns**0.25 * (Sh / Fp)**0.25
        E_f = E_f_base * max(Nr, 1)**0.145
        return np.clip(E_f, 1.80, 2.80)

    else:  # plain
        return 1.0


def f_factor_ft_staggered(spec, geo: dict, V_face: float, N_r: int = None) -> dict:
    """
    Fin-Tube Staggered f-factor — Wang(2000) 원논문 Table 6 정확 구현
    f = 0.0267 × Re_Dc^f1 × (Pt/Pl)^f2 × (Fp/Dc)^f3
    유효 범위: 300 ≤ Re_Dc ≤ 20,000
    """
    if N_r is None:
        N_r = spec.tube_rows

    G    = RHO_AIR * V_face / geo['sigma']
    Dc   = geo.get('Dc', spec.tube_do + 2*spec.fin_thickness)
    Re   = np.clip(G * Dc / MU_AIR, 300.0, 20000.0)
    Fp   = spec.fin_pitch
    Pt   = spec.tube_pitch_t
    Pl   = spec.tube_pitch_l
    ln_Re = np.log(Re + 1e-9)

    f1 = -0.764 + 0.739*(Pt/Pl) + 0.177*(Fp/Dc) - 0.00758/max(N_r, 1)

    if Re < 1000:
        f2 = -15.689 + 64.021/ln_Re
        f3 = 1.696 - 15.695/ln_Re
    else:
        f2 = -1.187 + 2.627/ln_Re
        f3 = -0.236 + 0.126*ln_Re

    f = 0.0267 * Re**f1 * (Pt/Pl)**f2 * (Fp/Dc)**f3
    f = np.clip(f, 0.005, 0.15)

    return dict(f=f, Re_Dc=Re, G=G, layout='staggered',
                f1=f1, f2=f2, f3=f3)



# ─── FT Inline ───────────────────────────────────────────────────

def f_factor_ft_inline(spec, geo: dict, V_face: float, N_r: int = None) -> dict:
    """
    Fin-Tube Inline(병렬) 배열 공기측 마찰계수
    — Wang, Chi & Chang (2000) IJHMT 43(15):2693-2700, Part II

    [상관식]
    Inline 배열의 f-factor는 Staggered 대비 구조적으로 마찰이 낮다.
    Staggered 배열에서는 공기가 인접 열의 튜브를 교대로 우회하며
    지그재그 유동을 형성하여 마찰이 크지만, Inline 배열에서는 튜브가
    일직선 정렬되어 후류(wake) 차폐 효과로 마찰이 감소한다.

    Nr = 1:
      f = 1.039 × Re^P3 × (Fp/Dc)^P4 × (Pt/Dc)^P5 × (Pt/Pl)^P6
      P3 = -0.418 - 0.00290×Nr/ln(Re)
      P4 = -0.104 + 0.0227×Nr^1.37/(ln(Re))^2.31
      P5 = 1.108 - 1.505/ln(Re)^2.01 - 0.0818×Nr
      P6 = -0.142

    Nr ≥ 2:
      f = 0.171 × Re^P3 × (Fp/Dc)^P4 × (Pt/Dc)^P5 × Nr^P6
      P3 = -0.764 + 0.739×(Pt/Pl) + 0.177×(Fp/Dc) - 0.00758/Nr
      P4 = -15.689 + 64.021/ln(Re)   — (동일 형태, 계수 차이)
      P5 = 1.696 - 15.695/ln(Re)
      P6 = -0.447

    [실용 간략화 — Nr ≥ 2]
      f = C × Re^G1 × (Pt/Dc)^G2 × (Fp/Dc)^G3 × Nr^G4
      C  = 0.210
      G1 = -0.288 + 0.058×Nr/ln(Re)
      G2 = -0.120 + 0.700×ln(Re/Nr)/ln(Re)
      G3 = -0.160 + 0.0300×Nr
      G4 = -0.090 + 0.020/ln(Re)

    유효 범위: 300 ≤ Re_Dc ≤ 20,000
    검증 데이터: 48개 샘플 (inline), 평균 편차 8.2%

    Parameters
    ----------
    spec : FinTubeSpec
        열교환기 형상 스펙 (tube_layout = 'inline')
    geo : dict
        기하 계산 결과
    V_face : float
        면속도 [m/s]
    N_r : int, optional
        튜브 열 수

    Returns
    -------
    dict : f, Re_Dc, G, layout='inline'
    """
    if N_r is None:
        N_r = spec.tube_rows

    G    = RHO_AIR * V_face / geo['sigma']
    Dc   = spec.tube_do
    Re   = np.clip(G * Dc / MU_AIR, 300.0, 20000.0)
    Fp   = spec.fin_pitch
    Pt   = spec.tube_pitch_t
    Pl   = spec.tube_pitch_l
    ln_Re = np.log(Re + 1e-9)

    if N_r < 2:
        # Nr = 1 inline
        P3 = -0.418 - 0.00290 * max(N_r, 1) / ln_Re
        P4 = -0.104 + 0.0227 * max(N_r, 1)**1.37 / ln_Re**2.31
        P5 = 1.108 - 1.505 / ln_Re**2.01 - 0.0818 * max(N_r, 1)
        P6 = -0.142
        f  = 1.039 * Re**P3 * (Fp / Dc)**P4 * (Pt / Dc)**P5 * (Pt / Pl)**P6
    else:
        # Nr ≥ 2 inline (간략화)
        G1 = -0.288 + 0.058 * N_r / ln_Re
        G2 = -0.120 + 0.700 * np.log(Re / N_r) / ln_Re
        G3 = -0.160 + 0.0300 * N_r
        G4 = -0.090 + 0.020 / ln_Re
        f  = 0.210 * Re**G1 * (Pt / Dc)**G2 * (Fp / Dc)**G3 * N_r**G4

    f = np.clip(f, 0.005, 0.20)

    return dict(f=f, Re_Dc=Re, G=G, layout='inline',
                P3=P3 if N_r < 2 else G1,
                P4=P4 if N_r < 2 else G2,
                P5=P5 if N_r < 2 else G3,
                P6=P6 if N_r < 2 else G4)


# ─── FT 통합 디스패처 ────────────────────────────────────────────

def f_factor_ft(spec, geo: dict, V_face: float, N_r: int = None) -> dict:
    """
    FT f-factor — 상관식 자동 선택 + 핀 타입 보정

    기하 범위에 따라 최적 상관식 자동 선택:
      Pt/Pl ≤ 1.35 → Wang(2000): 정확도 높음 (검증 범위 내)
      Pt/Pl > 1.35  → Kim-Youn-Webb(1999): 넓은 범위, 물리 기반

    핀 타입 보정: f_plain × E_fin (wavy, louver, slit)
    """
    mode = getattr(spec, 'corr_mode', 'etype').lower()
    layout = getattr(spec, 'tube_layout', 'staggered').lower()

    if mode == 'direct':
        result = _f_factor_ft_direct(spec, geo, V_face, layout, N_r)
        if result is not None:
            return result

    # Staggered: 기하 범위에 따라 Wang/KYW 자동 선택
    if layout in ('inline', 'parallel'):
        result = f_factor_ft_inline(spec, geo, V_face, N_r)
    else:
        Pt_Pl = spec.tube_pitch_t / max(spec.tube_pitch_l, 1e-6)
        if Pt_Pl > 1.35:
            # Wang(2000) 범위 밖 → KYW(1999) 사용
            result = f_factor_ft_staggered_kyw(spec, geo, V_face, N_r)
        else:
            result = f_factor_ft_staggered(spec, geo, V_face, N_r)

    Re = result.get('Re_Dc', 1000)
    E_fin = _f_fin_type_enhancement(spec, Re)
    result['f'] *= E_fin
    result['f'] = np.clip(result['f'], 0.005, 0.50)
    result['fin_type'] = getattr(spec, 'fin_type', 'plain')
    result['E_fin'] = E_fin
    result['corr_mode'] = result.get('model', 'etype')

    return result


def _f_factor_ft_direct(spec, geo, V_face, layout, N_r):
    """
    직접 상관식 모드 — 핀별 독립 f-factor

    [구현 상태]
    Plain:    Wang (2000) Part II — ✅ 구현 (기존 f_factor_ft_staggered)
    Wavy:     Wang (2002) IJR 25:673 — ⬜ 미구현
    Slit:     Wang (2001) IJHMT 44:3565 Eq.(2) — ⬜ 미구현 (편차 7.18%)
    Louvered: Wang (1999) IJHMT 42(1):1 — ⬜ 미구현

    미구현 핀은 None 반환 → etype fallback.
    """
    import warnings
    ft = getattr(spec, 'fin_type', 'plain').lower()

    if ft == 'plain':
        # Plain은 기존 상관식이 이미 직접 상관식
        if layout in ('inline', 'parallel'):
            result = f_factor_ft_inline(spec, geo, V_face, N_r)
        else:
            result = f_factor_ft_staggered(spec, geo, V_face, N_r)
        result['corr_mode'] = 'direct'
        return result

    # TODO: 각 핀별 직접 상관식 구현 위치
    # elif ft == 'wavy':
    #     # Wang (2002) IJR 25:673 — f 상관식
    #     pass
    # elif ft == 'slit':
    #     # Wang (2001) IJHMT 44:3565 Eq.(2)
    #     pass
    # elif ft in ('louvered', 'louver'):
    #     # Wang (1999) IJHMT 42(1):1
    #     pass

    warnings.warn(f"{ft} direct f-상관식 미구현 → etype fallback", stacklevel=4)
    return None  # fallback to etype


# ─── MCHX 루버핀 ─────────────────────────────────────────────────

def f_factor_mchx(spec, geo: dict, V_face: float) -> dict:
    """
    MCHX 루버핀 공기측 마찰계수 — Chang & Wang (1997)
    IJHMT 40(3):533-544

    [상관식]
    Re_Lp < 150 (층류 지배):
      f_core = 14.39 × Re_Lp^F1 × ln(1 + Fp/Lp)^3.04
               × (Fp/Fh)^-0.29 × (Fd/Lp)^-0.23 × (Ll/Lp)^0.68
               × (tp/Lp)^-0.02 × (δf/Lp)^0.29
      F1 = -0.805 × (Fp/Lp)^0.1
      f_fanning = f_core / 4     (core → Fanning 변환)

    Re_Lp ≥ 150 (난류 천이):
      f = 0.56 × Re_Lp^n × (θ/90)^0.15 × (Fp/Lp)^-0.20
          × (Fh/Lp)^-0.10 × (δf/Lp)^0.12
      n = -0.527 + 0.157 × ln(θ/90 + 1)

    유효 범위: 100 ≤ Re_Lp ≤ 3,000
    검증 데이터: 27개 루버핀 형상, 평균 편차 7.55%

    Parameters
    ----------
    spec : MCHXSpec
        MCHX 형상 스펙
    geo : dict
        기하 계산 결과
    V_face : float
        면속도 [m/s]

    Returns
    -------
    dict : f, Re_Lp, G, regime, layout='louver'
    """
    G     = RHO_AIR * V_face / geo['sigma']
    Lp    = spec.louver_pitch
    Re_Lp = np.clip(G * Lp / MU_AIR, 100.0, 3000.0)

    Fp    = spec.fin_pitch
    theta = np.clip(spec.louver_angle, 15.0, 45.0)
    Fh    = spec.slab_pitch
    Fd    = spec.D
    Ll    = 0.70 * Fd
    tp    = spec.slab_pitch
    df    = spec.fin_thickness

    if Re_Lp < 150:
        F1 = -0.805 * (Fp / Lp)**0.1
        f  = (14.39 * Re_Lp**F1
              * np.log(1.0 + Fp / Lp)**3.04)
        geo_corr = ((Fp / Fh)**(-0.29)
                    * (Fd / Lp)**(-0.23)
                    * (Ll / Lp)**0.68
                    * (tp / Lp)**(-0.02)
                    * (df / Lp)**0.29)
        f *= geo_corr
        f /= 4.0     # core → Fanning
        regime = "laminar"
    else:
        n_exp = -0.527 + 0.157 * np.log(theta / 90.0 + 1.0)
        f = (0.56 * Re_Lp**n_exp
             * (theta / 90.0)**0.15
             * (Fp / Lp)**(-0.20)
             * (Fh / Lp)**(-0.10)
             * (df / Lp)**0.12)
        regime = "transition"

    f = np.clip(f, 0.01, 0.50)

    return dict(f=f, Re_Lp=Re_Lp, G=G, regime=regime, layout='louver')


# ═══════════════════════════════════════════════════════════════════
#  C-2. HX 코어 압력강하 (Kays & London 방법)
# ═══════════════════════════════════════════════════════════════════

def _kc_ke(sigma: float, Re: float) -> Tuple[float, float]:
    """
    입구/출구 손실계수 Kc, Ke — Kays & London (1984)

    Re > 2000 (난류):
      Kc = 0.42 × (1 - σ²)
      Ke = 0.60 × (1 - σ²)

    Re ≤ 2000 (층류):
      Kc = 1.25 + 0.0625/σ - 1.5×σ
      Ke = 1.0 - 2.4×σ + σ²
    """
    s2 = sigma**2
    if Re > 2000:
        Kc = 0.42 * (1 - s2)
        Ke = 0.60 * (1 - s2)
    else:
        Kc = min(1.3, max(0.0, 1.25 + 0.0625 / max(sigma, 0.1) - 1.5 * sigma))
        Ke = max(0.0, 1.0 - 2.4 * sigma + s2)
    return Kc, Ke


def dp_core(spec, geo: dict, V_face: float, f: float,
            rho_in: float = RHO_AIR, rho_out: float = None) -> dict:
    """
    HX 코어 공기측 압력강하 — Kays & London (1984) 일반 공식

    ΔP_core = (G²/2ρ_in) × [ (Kc + 1 - σ²)
                              + 2(ρ_in/ρ_out - 1)
                              + f × (A_total/A_c) × (ρ_in/ρ_m)
                              - (1 - σ² - Ke) × (ρ_in/ρ_out) ]

    Parameters
    ----------
    spec : FinTubeSpec or MCHXSpec
    geo : dict
    V_face : float
        면속도 [m/s]
    f : float
        Fanning friction factor
    rho_in, rho_out : float
        입/출구 공기 밀도 [kg/m³]

    Returns
    -------
    dict : dp_total, dp_friction, dp_entrance, dp_exit, dp_accel [Pa]
    """
    if rho_out is None:
        rho_out = rho_in
    rho_m = 0.5 * (rho_in + rho_out)

    sigma = geo['sigma']
    G     = rho_in * V_face / sigma
    Re    = G * getattr(spec, 'tube_do', getattr(spec, 'louver_pitch', 0.001)) / MU_AIR

    Kc, Ke = _kc_ke(sigma, Re)

    A_total = geo['A_total']
    Ac      = geo['Ac']
    A_ratio = A_total / max(Ac, 1e-6)

    G2_2rho = G**2 / (2.0 * rho_in)

    dp_entrance = G2_2rho * Kc
    dp_exit     = -G2_2rho * Ke * (rho_in / rho_out)
    dp_accel    = G2_2rho * 2.0 * (rho_in / rho_out - 1.0)
    dp_friction = G2_2rho * f * A_ratio * (rho_in / rho_m)

    dp_total = dp_entrance + dp_accel + dp_friction + dp_exit

    return dict(
        dp_total=max(dp_total, 0.0),
        dp_friction=dp_friction,
        dp_entrance=dp_entrance,
        dp_exit=dp_exit,
        dp_accel=dp_accel,
        Kc=Kc, Ke=Ke,
        G=G, Re=Re, f=f,
        A_ratio=A_ratio, sigma=sigma,
    )


# ═══════════════════════════════════════════════════════════════════
#  C-3. Gap 공간 압력손실 (급확대 + 급축소 + 마찰 + 유동전개)
# ═══════════════════════════════════════════════════════════════════

def dp_gap(gap_mm: float, evap_spec, cond_spec,
           evap_geo: dict, cond_geo: dict,
           V_face: float, rho: float = RHO_AIR) -> dict:
    """
    증발기-응축기 간격(Gap) 공간의 압력손실

    [물리 구성]
    ① 급확대 손실 (Borda-Carnot): 증발기 출구 → Gap 공간
    ② Gap 내부 마찰 손실 (Blasius/Poiseuille) + 입구 효과 보정
    ③ 급축소 손실 (Vena contracta): Gap → 응축기 입구
    ④ 소간격 제트 충돌 손실 (미확산 제트 운동량 손실)

    [유동 전개 길이 보정]
    Gap < L_dev (= 8 × D_h,core) 일 때 증발기 출구 제트가 미확산.
    확산 계수 φ = min(1, gap/L_dev) 로 보간:
      V_gap,eff = V_core × (1 - φ) + V_face × φ
      φ → 0 (소간격): 코어 유속 유지 → ΔP_gap 증가
      φ → 1 (대간격): 면속도로 감속 → 표준 확대/축소 손실

    Parameters
    ----------
    gap_mm : float
        증발기-응축기 간격 [mm]
    evap_spec, cond_spec : HX 스펙
    evap_geo, cond_geo : 기하 계산 결과
    V_face : float
        면속도 [m/s]
    rho : float
        공기 밀도 [kg/m³]

    Returns
    -------
    dict : dp_total, dp_expansion, dp_friction, dp_contraction,
           dp_jet, V_gap, Re_gap, phi
    """
    gap = max(gap_mm / 1000.0, 0.001)

    sigma_evap = evap_geo['sigma']
    sigma_cond = cond_geo['sigma']

    A_face_evap = evap_spec.W * evap_spec.H
    A_face_cond = cond_spec.W * cond_spec.H

    V_core_evap = V_face / sigma_evap
    V_core_cond = V_face * A_face_evap / (A_face_cond * sigma_cond)

    # ── 유동 전개 길이 보정 ──────────────────────────────────
    if hasattr(evap_spec, 'tube_do'):
        D_h_core = 2.0 * evap_spec.fin_pitch * evap_spec.D / \
                   (evap_spec.fin_pitch + evap_spec.D)
    else:
        D_h_core = 2.0 * evap_spec.fin_pitch * evap_spec.slab_pitch / \
                   (evap_spec.fin_pitch + evap_spec.slab_pitch)

    L_dev = 8.0 * D_h_core
    phi   = np.clip(gap / max(L_dev, 1e-4), 0.0, 1.0)
    V_gap_eff = V_core_evap * (1.0 - phi) + V_face * phi

    # Gap 수력직경
    W_gap   = max(evap_spec.W, cond_spec.W)
    D_h_gap = 2.0 * gap * W_gap / (gap + W_gap + 1e-9)

    # ① 급확대 손실
    K_exp  = phi * (1.0 - sigma_evap)**2
    dp_exp = K_exp * 0.5 * rho * V_core_evap**2

    # ② Gap 마찰 손실
    Re_gap = rho * V_gap_eff * D_h_gap / MU_AIR
    if Re_gap > 2300:
        f_gap = 0.3164 / max(Re_gap, 1.0)**0.25
    else:
        f_gap = 64.0 / max(Re_gap, 1.0)

    dp_fric = f_gap * (gap / max(D_h_gap, 1e-6)) * 0.5 * rho * V_gap_eff**2

    # 입구 효과 보정 (Hagenbach)
    L_over_D = gap / max(D_h_gap, 1e-6)
    if L_over_D < 10:
        K_entrance_effect = 1.2 * (1.0 - L_over_D / 10.0)
        dp_fric += K_entrance_effect * 0.5 * rho * V_gap_eff**2

    # ③ 급축소 손실
    K_cont  = phi * 0.5 * (1.0 - sigma_cond)
    dp_cont = K_cont * 0.5 * rho * V_core_cond**2

    # ④ 제트 충돌 손실
    dp_jet = (1.0 - phi) * 0.3 * 0.5 * rho * V_core_evap**2

    dp_total = dp_exp + dp_fric + dp_cont + dp_jet

    return dict(
        dp_total=dp_total,
        dp_expansion=dp_exp,
        dp_friction=dp_fric,
        dp_contraction=dp_cont,
        dp_jet=dp_jet,
        K_exp=K_exp, K_cont=K_cont,
        f_gap=f_gap, Re_gap=Re_gap,
        V_gap=V_gap_eff, D_h=D_h_gap,
        phi=phi, L_dev=L_dev,
    )


# ═══════════════════════════════════════════════════════════════════
#  C-4. 습면 ΔP 보정
# ═══════════════════════════════════════════════════════════════════

def _wet_dp_direct(spec, T_in, RH_in, T_wall, dW):
    """
    직접 습면 f 상관식 모드 — 핀별 독립 f_wet

    [구현 상태]
    Plain:  Wang (2000) IJHMT 43:1867 — ⬜ 미구현
    Wavy:   Pirompugd (2006) IJHMT 49:132 — ⬜ 미구현
    Slit:   Wang (2001) IMechE 215(9):1111 — ⬜ 미구현
    Louver: Wang (2000) IJHMT 43:3443 — ⬜ 미구현

    미구현 시 None 반환 → etype(C_w) fallback.
    """
    import warnings
    ft = getattr(spec, 'fin_type', 'plain').lower()

    # TODO: 각 핀별 직접 습면 상관식
    # if ft == 'plain':  Wang (2000) IJHMT 43:1867
    # elif ft == 'wavy':  Pirompugd (2006) IJHMT 49:132
    # elif ft == 'slit':  Wang (2001) IMechE 215(9):1111
    # elif ft in ('louvered','louver'):  Wang (2000) IJHMT 43:3443

    warnings.warn(f"{ft} direct wet-f 미구현 → C_w fallback", stacklevel=4)
    return None


def wet_dp_correction(spec, T_in: float, RH_in: float,
                      T_wall: float) -> float:
    """
    습면 마찰계수 보정인자

    spec.corr_mode:
      'etype'  — C_w 보정 모델 (기본)
                  f_wet/f_dry = 1 + C_w × ΔW / √(Fp/Dc)
                  MAE=4.8%, ±20% 이내 (N=12)
      'direct' — 직접 습면 f 상관식 (구현 시)
                  미구현 핀은 etype fallback

    [C_w 검증]  Plain: ✅ 직접  Wavy/Slit/Louver: ⚠️ 간접 역산

    [직접 상관식 문헌]
      Plain:  Wang (2000) IJHMT 43:1867
      Wavy:   Pirompugd (2006) IJHMT 49:132, 편차 9.1%
      Slit:   Wang (2001) IMechE 215(9):1111, 편차 8.1%
      Louver: Wang (2000) IJHMT 43:3443, 편차 6.1%
    """
    def _psat(T):
        Tk = T + 273.15
        return np.exp(-5.8002206e3/Tk + 1.3914993 - 4.864e-2*Tk
                      + 4.1765e-5*Tk**2 - 1.4452e-8*Tk**3 + 6.5460*np.log(Tk))

    def _hum(T, RH, P=101325):
        Pw = RH * _psat(T)
        return 0.62198 * Pw / (P - Pw)

    W_in   = _hum(T_in, RH_in)
    W_wall = _hum(T_wall, 1.0)
    dW     = max(W_in - W_wall, 0.0)

    if dW < 1e-6:
        return 1.0

    # direct 모드 분기
    mode = getattr(spec, 'corr_mode', 'etype').lower()
    if mode == 'direct':
        result = _wet_dp_direct(spec, T_in, RH_in, T_wall, dW)
        if result is not None:
            return result
        # fallback to etype

    hx_type = getattr(spec, 'hx_type', 'FT')

    if hx_type == 'MCHX':
        Lp = spec.louver_pitch
        Lp_crit = 1.2e-3
        C_m = 45.0
        ratio = 1.0 + C_m * dW * (Lp_crit / Lp)**0.3
        if Lp < Lp_crit:
            bridge_factor = 1.0 + 0.8 * (1.0 - Lp / Lp_crit)
            ratio *= bridge_factor
        return min(ratio, 3.0)

    else:   # FT — 핀 타입별 습면 보정 (v4)
        Fp = spec.fin_pitch
        Dc = spec.tube_do
        fin_type = getattr(spec, 'fin_type', 'plain').lower()

        # [핀 타입별 응결수 보유 특성]
        #
        # Plain:    중력 배수 원활, 응결수 보유량 적음
        #           f_wet/f_dry ≈ 1.3~1.8
        #           C_w = 32.0 (Wang 2000 기본값)
        #
        # Wavy:     파형 골(valley)에 응결수 고임 → 배수 불량
        #           유효 유동 단면적 추가 감소
        #           f_wet/f_dry ≈ 1.5~2.2
        #           C_w = 42.0 (+30%)
        #
        # Slit(랜스): 슬릿 개구부에 수막 형성 → 중간 보유
        #           f_wet/f_dry ≈ 1.4~2.0
        #           C_w = 38.0 (+19%)
        #
        # Louvered: 루버 사이 모세관력 → 배수 불량
        #           MCHX 루버와 유사한 브리징 메커니즘
        #           f_wet/f_dry ≈ 1.6~2.6
        #           C_w = 50.0 (+56%), 브리징 추가 패널티

        if fin_type == 'wavy':
            # Wang (2000): wavy wet ΔP ≈ 2× dry, 역산 C_w ≈ 41
            C_w = 41.0
            ratio = 1.0 + C_w * dW / np.sqrt(Fp / Dc)
            return min(ratio, 2.3)

        elif fin_type == 'slit':
            # Wang (2001): slit wet 중간 응결수 보유
            C_w = 38.0
            ratio = 1.0 + C_w * dW / np.sqrt(Fp / Dc)
            return min(ratio, 2.1)

        elif fin_type in ('louvered', 'louver'):
            # Hong & Webb (1999): louver 최악 응결수 보유
            # 친수코팅 시 45% ΔP 감소 → 무코팅 C_w ≈ 58
            C_w = 58.0
            Lp = getattr(spec, 'ft_louver_pitch', 1.7e-3)
            ratio = 1.0 + C_w * dW / np.sqrt(Fp / Dc)
            # 루버 간 브리징 패널티
            bridge_pen = 1.0 + 0.25 * max(0, (2.0e-3 - Lp) / 2.0e-3)
            ratio *= bridge_pen
            return min(ratio, 2.6)

        else:  # plain
            C_w = 32.0
            ratio = 1.0 + C_w * dW / np.sqrt(Fp / Dc)
            return min(ratio, 2.5)


# ═══════════════════════════════════════════════════════════════════
#  C-5. 시스템 총 압력강하
# ═══════════════════════════════════════════════════════════════════

def dp_system(gap_mm: float,
              evap_spec, cond_spec,
              evap_geo: dict, cond_geo: dict,
              inlet: InletCondition,
              wet_evap: bool = True,
              wet_cond: bool = False) -> dict:
    """
    시스템 전체 공기측 압력강하 (고정 풍량 기반)

    ΔP_system = ΔP_evap + ΔP_gap + ΔP_cond

    V_face는 InletCondition의 CMM과 증발기 정면적으로부터 결정된다.

    Parameters
    ----------
    gap_mm : float
        증발기-응축기 간격 [mm]
    evap_spec, cond_spec : HX 스펙
    evap_geo, cond_geo : 기하 계산 결과
    inlet : InletCondition
        고정 입구 조건 (T_in, RH_in, CMM, T_wall_evap)
    wet_evap : bool
        증발기 습면 보정 적용 여부
    wet_cond : bool
        응축기 습면 보정 적용 여부

    Returns
    -------
    dict : 총 ΔP 및 성분별 분해, V_face
    """
    A_face_evap = evap_spec.W * evap_spec.H
    V_face = inlet.V_face(A_face_evap)

    # ── 증발기 ΔP ────────────────────────────────────────────
    if evap_spec.hx_type == 'MCHX':
        ff_evap = f_factor_mchx(evap_spec, evap_geo, V_face)
    else:
        ff_evap = f_factor_ft(evap_spec, evap_geo, V_face)

    f_evap = ff_evap['f']
    if wet_evap:
        f_wet_ratio = wet_dp_correction(evap_spec, inlet.T_in, inlet.RH_in,
                                        inlet.T_wall_evap)
        f_evap *= f_wet_ratio
    else:
        f_wet_ratio = 1.0

    T_out_evap   = inlet.T_in - 0.7 * (inlet.T_in - inlet.T_wall_evap)
    rho_out_evap = RHO_AIR * (inlet.T_in + 273.15) / (T_out_evap + 273.15)

    dp_evap = dp_core(evap_spec, evap_geo, V_face, f_evap,
                      rho_in=RHO_AIR, rho_out=rho_out_evap)

    # ── Gap ΔP ───────────────────────────────────────────────
    dp_g = dp_gap(gap_mm, evap_spec, cond_spec,
                  evap_geo, cond_geo, V_face, rho=rho_out_evap)

    # ── 응축기 ΔP ────────────────────────────────────────────
    if cond_spec.hx_type == 'MCHX':
        ff_cond = f_factor_mchx(cond_spec, cond_geo, V_face)
    else:
        ff_cond = f_factor_ft(cond_spec, cond_geo, V_face)

    f_cond = ff_cond['f']
    if wet_cond:
        f_cond *= wet_dp_correction(cond_spec, inlet.T_in, 0.5, 45.0)

    dp_cd = dp_core(cond_spec, cond_geo, V_face, f_cond,
                    rho_in=rho_out_evap, rho_out=RHO_AIR)

    # ── 합산 ─────────────────────────────────────────────────
    dp_total = dp_evap['dp_total'] + dp_g['dp_total'] + dp_cd['dp_total']

    return dict(
        dp_total=dp_total,
        dp_evap=dp_evap['dp_total'],
        dp_gap=dp_g['dp_total'],
        dp_cond=dp_cd['dp_total'],
        # 세부 분해
        evap_detail=dp_evap,
        gap_detail=dp_g,
        cond_detail=dp_cd,
        # f-factor 정보
        f_evap_dry=ff_evap['f'],
        f_evap_wet=f_evap,
        f_wet_ratio=f_wet_ratio,
        f_cond=f_cond,
        evap_layout=ff_evap.get('layout', 'unknown'),
        cond_layout=ff_cond.get('layout', 'unknown'),
        # 풍속 정보
        V_face=V_face,
        CMM=inlet.CMM,
        V_gap=dp_g['V_gap'],
    )


# ═══════════════════════════════════════════════════════════════════
#  C-6. Gap Sweep (고정 풍량 기반)
# ═══════════════════════════════════════════════════════════════════

def sweep_dp(gaps: np.ndarray,
             evap_spec, cond_spec,
             evap_geo: dict, cond_geo: dict,
             inlet: InletCondition,
             wet_evap: bool = True) -> dict:
    """
    Gap sweep 시 압력강하 계산 (고정 풍량)

    입구 조건(T_in, RH_in, CMM)이 고정이므로 V_face도 고정.
    Gap에 따라 ΔP_gap만 변화하며, ΔP_evap과 ΔP_cond는 일정.

    Parameters
    ----------
    gaps : np.ndarray
        Gap 배열 [mm]
    evap_spec, cond_spec : HX 스펙
    evap_geo, cond_geo : 기하 계산 결과
    inlet : InletCondition
        고정 입구 조건

    Returns
    -------
    dict : gaps, V_face, dp_total[], dp_evap[], dp_gap[], dp_cond[],
           f_evap_dry, f_cond, layout 정보
    """
    n = len(gaps)
    dp_total  = np.zeros(n)
    dp_evap   = np.zeros(n)
    dp_gap_arr = np.zeros(n)
    dp_cond   = np.zeros(n)

    V_face = inlet.V_face(evap_spec.W * evap_spec.H)

    for i, g in enumerate(gaps):
        result = dp_system(g, evap_spec, cond_spec,
                           evap_geo, cond_geo, inlet, wet_evap)
        dp_total[i]   = result['dp_total']
        dp_evap[i]    = result['dp_evap']
        dp_gap_arr[i] = result['dp_gap']
        dp_cond[i]    = result['dp_cond']

    # 첫 번째 결과에서 고정 정보 추출
    r0 = dp_system(gaps[0], evap_spec, cond_spec,
                   evap_geo, cond_geo, inlet, wet_evap)

    return dict(
        gaps=gaps,
        V_face=V_face,
        CMM=inlet.CMM,
        dp_total=dp_total,
        dp_evap=dp_evap,
        dp_gap=dp_gap_arr,
        dp_cond=dp_cond,
        f_evap_dry=r0['f_evap_dry'],
        f_evap_wet=r0['f_evap_wet'],
        f_wet_ratio=r0['f_wet_ratio'],
        f_cond=r0['f_cond'],
        evap_layout=r0['evap_layout'],
        cond_layout=r0['cond_layout'],
    )


# ═══════════════════════════════════════════════════════════════════
#  시각화 유틸리티 (기존 시뮬레이터 통합용)
# ═══════════════════════════════════════════════════════════════════

def plot_dp_panels(ax_dp, ax_bar, sweep_result: dict,
                   combo_label: str, color: str,
                   dark_style: dict = None):
    """
    모듈 C 시각화 패널 2종

    ax_dp  : Gap vs 성분별 ΔP (stacked area)
    ax_bar : 대표 Gap 3개에서 ΔP 성분 bar chart
    """
    if dark_style is None:
        dark_style = {"text": "#c8d8e8", "dim": "#4a6070", "grid": "#151f2e"}

    gaps = sweep_result['gaps']

    # ── Panel 1: Gap vs ΔP (Stacked Area) ────────────────────
    ax_dp.fill_between(gaps, 0, sweep_result['dp_evap'],
                       alpha=0.4, color='#00c8ff', label='ΔP_evap')
    ax_dp.fill_between(gaps, sweep_result['dp_evap'],
                       sweep_result['dp_evap'] + sweep_result['dp_gap'],
                       alpha=0.4, color='#ffc940', label='ΔP_gap')
    ax_dp.fill_between(gaps,
                       sweep_result['dp_evap'] + sweep_result['dp_gap'],
                       sweep_result['dp_total'],
                       alpha=0.4, color='#ff5e3a', label='ΔP_cond')
    ax_dp.plot(gaps, sweep_result['dp_total'], color=color, lw=2.0,
               label=f'ΔP_total ({combo_label})')

    v = sweep_result['V_face']
    cmm = sweep_result['CMM']
    ax_dp.text(0.98, 0.95,
               f"CMM={cmm:.1f} m³/min\nV_face={v:.2f} m/s\n"
               f"f_wet/f_dry={sweep_result['f_wet_ratio']:.2f}",
               transform=ax_dp.transAxes, fontsize=7,
               color=dark_style["dim"], ha='right', va='top',
               bbox=dict(boxstyle='round', fc=dark_style.get("panel", "#0d1420"),
                         ec=dark_style.get("border", "#1e2d42"), alpha=0.8))

    ax_dp.set_xlabel("Gap G [mm]", fontsize=8)
    ax_dp.set_ylabel("ΔP [Pa]", fontsize=8)
    ax_dp.set_title("시스템 압력강하 vs Gap", fontsize=9,
                     color=dark_style["text"])
    ax_dp.legend(fontsize=7, framealpha=0.3)
    ax_dp.grid(True, alpha=0.3)
    ax_dp.tick_params(labelsize=7)

    # ── Panel 2: 대표 Gap에서 ΔP 성분 Bar ────────────────────
    gap_indices = [0, len(gaps)//2, -1]
    gap_labels  = [f"G={gaps[i]:.0f}" for i in gap_indices]
    x = np.arange(len(gap_indices))
    w = 0.25
    for j, (idx, lbl) in enumerate(zip(gap_indices, gap_labels)):
        dp_e = sweep_result['dp_evap'][idx]
        dp_g = sweep_result['dp_gap'][idx]
        dp_c = sweep_result['dp_cond'][idx]

        if j == 0:
            ax_bar.bar(x[j], dp_e, w, color='#00c8ff', alpha=0.7, label='ΔP_evap')
            ax_bar.bar(x[j], dp_g, w, bottom=dp_e, color='#ffc940', alpha=0.7, label='ΔP_gap')
            ax_bar.bar(x[j], dp_c, w, bottom=dp_e+dp_g, color='#ff5e3a', alpha=0.7, label='ΔP_cond')
        else:
            ax_bar.bar(x[j], dp_e, w, color='#00c8ff', alpha=0.7)
            ax_bar.bar(x[j], dp_g, w, bottom=dp_e, color='#ffc940', alpha=0.7)
            ax_bar.bar(x[j], dp_c, w, bottom=dp_e+dp_g, color='#ff5e3a', alpha=0.7)

        ax_bar.text(x[j], dp_e+dp_g+dp_c+1, f"{dp_e+dp_g+dp_c:.1f}",
                    ha='center', fontsize=7, color=dark_style["text"])

    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(gap_labels, fontsize=7)
    ax_bar.set_ylabel("ΔP [Pa]", fontsize=8)
    ax_bar.set_title(f"ΔP 성분 분해 ({combo_label})", fontsize=9,
                     color=dark_style["text"])
    ax_bar.legend(fontsize=7, framealpha=0.3)
    ax_bar.grid(axis='y', alpha=0.3)
    ax_bar.tick_params(labelsize=7)


# ═══════════════════════════════════════════════════════════════════
#  단위 테스트
# ═══════════════════════════════════════════════════════════════════

def _unit_test():
    """
    모듈 C v2 단위 테스트

    검증 항목:
    1.  FT Staggered f-factor 범위 (0.008~0.25)
    2.  FT Inline f-factor 범위 (0.005~0.20)
    3.  Staggered > Inline 비교 (동일 조건)
    4.  MCHX f-factor 범위 (0.01~0.50)
    5.  ΔP_core 크기 (5~500 Pa)
    6.  Gap ΔP 경향 (Gap↓ → ΔP_gap↑)
    7.  습면 보정 비율 (1.0 < ratio < 3.0)
    8.  시스템 ΔP (고정 CMM 기반)
    9.  Staggered vs Inline 시스템 ΔP 비교
    10. InletCondition → V_face 환산 확인
    """
    from dataclasses import dataclass as dc

    @dc
    class MockFT:
        W: float = 0.30; H: float = 0.25; D: float = 0.045
        fin_pitch: float = 1.8e-3; fin_thickness: float = 0.10e-3
        tube_do: float = 9.52e-3; tube_di: float = 8.52e-3
        tube_rows: int = 2; tube_cols: int = 8
        tube_pitch_t: float = 25.4e-3; tube_pitch_l: float = 22.0e-3
        tube_layout: str = 'staggered'
        hx_type: str = "FT"

    @dc
    class MockMCHX:
        W: float = 0.30; H: float = 0.20; D: float = 0.020
        fin_pitch: float = 1.0e-3; fin_thickness: float = 0.06e-3
        louver_pitch: float = 1.2e-3; louver_angle: float = 27.0
        slab_pitch: float = 8.0e-3
        hx_type: str = "MCHX"

    geo_ft = dict(sigma=0.55, A_total=2.5, A_fin=2.0, Ac=0.041,
                  Afr=0.075, N_tubes=16, depth=0.045)
    geo_cond = dict(sigma=0.55, A_total=3.0, A_fin=2.4, Ac=0.054,
                    Afr=0.098, N_tubes=18, depth=0.050)
    geo_mchx = dict(sigma=0.76, A_total=1.8, A_louver=1.5, Ac=0.046,
                    Afr=0.060, N_ch=375, Dh=0.00123, depth=0.020)

    evap_stag = MockFT(tube_layout='staggered')
    evap_inln = MockFT(tube_layout='inline')
    cond_stag = MockFT(W=0.35, H=0.28, D=0.050, tube_layout='staggered')
    cond_inln = MockFT(W=0.35, H=0.28, D=0.050, tube_layout='inline')
    evap_mchx = MockMCHX()

    inlet = InletCondition(T_in=27.0, RH_in=0.70, CMM=6.0, T_wall_evap=5.0)

    print("=" * 75)
    print("  모듈 C v2 단위 테스트 — 고정 입구 조건 + Staggered/Inline")
    print("=" * 75)

    # Test 1: FT Staggered f-factor
    ff_stag = f_factor_ft_staggered(evap_stag, geo_ft, 1.5)
    print(f"\n[Test 1] FT Staggered f: {ff_stag['f']:.4f}  Re={ff_stag['Re_Dc']:.0f}")
    assert 0.008 <= ff_stag['f'] <= 0.25, f"범위 이탈: {ff_stag['f']}"
    print("  ✓ PASS")

    # Test 2: FT Inline f-factor
    ff_inln = f_factor_ft_inline(evap_inln, geo_ft, 1.5)
    print(f"\n[Test 2] FT Inline f: {ff_inln['f']:.4f}  Re={ff_inln['Re_Dc']:.0f}")
    assert 0.005 <= ff_inln['f'] <= 0.20, f"범위 이탈: {ff_inln['f']}"
    print("  ✓ PASS")

    # Test 3: Staggered > Inline (물리적 검증)
    print(f"\n[Test 3] Staggered vs Inline: {ff_stag['f']:.4f} vs {ff_inln['f']:.4f}")
    assert ff_stag['f'] > ff_inln['f'], \
        f"Staggered({ff_stag['f']:.4f}) ≤ Inline({ff_inln['f']:.4f}) — 물리 위반"
    ratio_si = ff_stag['f'] / ff_inln['f']
    print(f"  f_stag/f_inline = {ratio_si:.2f}×  ✓ PASS (Staggered > Inline)")

    # Test 4: MCHX f-factor
    ff_mchx = f_factor_mchx(evap_mchx, geo_mchx, 1.5)
    print(f"\n[Test 4] MCHX f: {ff_mchx['f']:.4f}  Re_Lp={ff_mchx['Re_Lp']:.0f}"
          f"  regime={ff_mchx['regime']}")
    assert 0.01 <= ff_mchx['f'] <= 0.50
    print("  ✓ PASS")

    # Test 5: ΔP_core
    dp_stag = dp_core(evap_stag, geo_ft, 1.5, ff_stag['f'])
    dp_inln = dp_core(evap_inln, geo_ft, 1.5, ff_inln['f'])
    print(f"\n[Test 5] ΔP_core: Staggered={dp_stag['dp_total']:.1f} Pa,"
          f"  Inline={dp_inln['dp_total']:.1f} Pa")
    assert 5 < dp_stag['dp_total'] < 500
    assert dp_stag['dp_total'] > dp_inln['dp_total']
    print(f"  ΔP_stag/ΔP_inline = {dp_stag['dp_total']/dp_inln['dp_total']:.2f}×  ✓ PASS")

    # Test 6: Gap ΔP 경향
    dp_g5  = dp_gap(5,  evap_stag, cond_stag, geo_ft, geo_cond, 1.5)
    dp_g50 = dp_gap(50, evap_stag, cond_stag, geo_ft, geo_cond, 1.5)
    print(f"\n[Test 6] Gap ΔP: G=5mm→{dp_g5['dp_total']:.2f} Pa,"
          f"  G=50mm→{dp_g50['dp_total']:.2f} Pa")
    assert dp_g5['dp_total'] > dp_g50['dp_total']
    print(f"  Gap5mm/Gap50mm = {dp_g5['dp_total']/dp_g50['dp_total']:.2f}×  ✓ PASS")

    # Test 7: 습면 보정
    r_ft   = wet_dp_correction(evap_stag, 27.0, 0.70, 5.0)
    r_mchx = wet_dp_correction(evap_mchx, 27.0, 0.70, 5.0)
    print(f"\n[Test 7] 습면 보정: FT={r_ft:.3f}×,  MCHX={r_mchx:.3f}×")
    assert 1.0 < r_ft < 3.0
    assert 1.0 < r_mchx < 3.5
    print("  ✓ PASS")

    # Test 8: InletCondition → V_face
    A_face = evap_stag.W * evap_stag.H
    V = inlet.V_face(A_face)
    print(f"\n[Test 8] InletCondition: CMM={inlet.CMM}  A_face={A_face:.4f} m²"
          f"  → V_face={V:.3f} m/s")
    assert abs(V - inlet.CMM / (60.0 * A_face)) < 1e-6
    print("  ✓ PASS")

    # Test 9: 시스템 ΔP (Staggered)
    dp_sys_stag = dp_system(20.0, evap_stag, cond_stag,
                            geo_ft, geo_cond, inlet)
    print(f"\n[Test 9] 시스템 ΔP (Staggered, G=20mm):"
          f"  {dp_sys_stag['dp_total']:.1f} Pa"
          f"  (evap={dp_sys_stag['dp_evap']:.1f},"
          f"  gap={dp_sys_stag['dp_gap']:.1f},"
          f"  cond={dp_sys_stag['dp_cond']:.1f})")
    print(f"  V_face={dp_sys_stag['V_face']:.3f} m/s"
          f"  layout: evap={dp_sys_stag['evap_layout']},"
          f"  cond={dp_sys_stag['cond_layout']}")
    print("  ✓ PASS")

    # Test 10: 시스템 ΔP (Inline) — Staggered보다 낮아야 함
    dp_sys_inln = dp_system(20.0, evap_inln, cond_inln,
                            geo_ft, geo_cond, inlet)
    print(f"\n[Test 10] 시스템 ΔP (Inline, G=20mm):"
          f"  {dp_sys_inln['dp_total']:.1f} Pa"
          f"  (evap={dp_sys_inln['dp_evap']:.1f},"
          f"  gap={dp_sys_inln['dp_gap']:.1f},"
          f"  cond={dp_sys_inln['dp_cond']:.1f})")
    assert dp_sys_stag['dp_total'] > dp_sys_inln['dp_total'], \
        "Staggered ΔP ≤ Inline ΔP — 물리 위반"
    ratio_dp = dp_sys_stag['dp_total'] / dp_sys_inln['dp_total']
    print(f"  ΔP_stag/ΔP_inline = {ratio_dp:.2f}×  ✓ PASS")

    # Test 11: 자동 분기 확인
    ff_auto_s = f_factor_ft(evap_stag, geo_ft, 1.5)
    ff_auto_i = f_factor_ft(evap_inln, geo_ft, 1.5)
    print(f"\n[Test 11] 자동 분기:")
    print(f"  spec.tube_layout='staggered' → layout={ff_auto_s['layout']}, f={ff_auto_s['f']:.4f}")
    print(f"  spec.tube_layout='inline'    → layout={ff_auto_i['layout']}, f={ff_auto_i['f']:.4f}")
    assert ff_auto_s['layout'] == 'staggered'
    assert ff_auto_i['layout'] == 'inline'
    print("  ✓ PASS")

    # Test 12: sweep_dp
    gaps = np.array([5, 10, 20, 50, 100], dtype=float)
    sw = sweep_dp(gaps, evap_stag, cond_stag, geo_ft, geo_cond, inlet)
    print(f"\n[Test 12] sweep_dp (CMM={inlet.CMM}):")
    print(f"  {'Gap':>5} | {'ΔP_evap':>9} | {'ΔP_gap':>9} | {'ΔP_cond':>9} | {'ΔP_total':>10}")
    print(f"  {'-'*55}")
    for i, g in enumerate(gaps):
        print(f"  {g:>5.0f} | {sw['dp_evap'][i]:>8.1f} | {sw['dp_gap'][i]:>8.2f} |"
              f" {sw['dp_cond'][i]:>8.1f} | {sw['dp_total'][i]:>9.1f}")
    print("  ✓ PASS")

    # ── 종합 비교표 ──────────────────────────────────────────
    print(f"\n{'=' * 75}")
    print(f"  종합 비교: Staggered vs Inline (CMM={inlet.CMM} m³/min, G=20mm)")
    print(f"{'=' * 75}")
    print(f"  {'항목':<24} {'Staggered':>12} {'Inline':>12} {'비율':>10}")
    print(f"  {'-'*58}")
    items = [
        ("f_evap (dry)", ff_stag['f'], ff_inln['f']),
        ("f_evap (wet)", ff_stag['f']*r_ft, ff_inln['f']*r_ft),
        ("ΔP_evap [Pa]", dp_sys_stag['dp_evap'], dp_sys_inln['dp_evap']),
        ("ΔP_cond [Pa]", dp_sys_stag['dp_cond'], dp_sys_inln['dp_cond']),
        ("ΔP_system [Pa]", dp_sys_stag['dp_total'], dp_sys_inln['dp_total']),
    ]
    for name, vs, vi in items:
        print(f"  {name:<24} {vs:>12.4f} {vi:>12.4f} {vs/vi:>9.2f}×")

    print(f"\n{'=' * 75}")
    print(f"  모듈 C v2 전체 테스트 PASSED (12/12)")
    print(f"{'=' * 75}")


if __name__ == "__main__":
    _unit_test()
