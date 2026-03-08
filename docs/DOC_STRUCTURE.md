# 증발기-응축기 간격 통합 시뮬레이터 — 모델 구조

## 1. 시뮬레이터 개요

증발기와 응축기 사이 간격(Gap)이 냉방 성능, 비말동반 위험, 압력강하에 미치는 영향을 정량 분석하는 통합 시뮬레이터.

**해석 범위**: Gap = 5~100mm, 4종 핀(plain/wavy/slit/louvered), 2종 배열(staggered/inline), FT/MCHX 지원

**핵심 출력**: Gap에 따른 Cap Retention [%], Carryover Risk, ΔP_system [Pa]


## 2. 아키텍처

```
                    ┌──────────────────────────────────┐
                    │         INPUT (JSON/GUI)          │
                    │  HX Spec, Air, Refrigerant, Gap   │
                    └─────────────┬────────────────────┘
                                  │
                    ┌─────────────▼────────────────────┐
                    │     common.py — 공통 기반          │
                    │                                   │
                    │  FinTubeSpec / MCHXSpec (기하)     │
                    │  RefrigerantState (CoolProp 물성)  │
                    │  compute_ft_geometry()             │
                    │  compute_UA() — h_o, h_i, η_o     │
                    │  compute_coil_performance()        │
                    │  compute_coil_performance_segmented│
                    │  ↑ Threlkeld 정석 열전달 엔진      │
                    └───┬──────────┬──────────┬─────────┘
                        │          │          │
              ┌─────────▼──┐ ┌────▼─────┐ ┌──▼─────────┐
              │ Module A   │ │ Module B │ │ Module C   │
              │ Gap 손실   │ │ 비말동반  │ │ 압력강하    │
              │ module_a.py│ │module_b.py│ │module_c.py │
              └─────┬──────┘ └────┬─────┘ └──┬─────────┘
                    │             │           │
              ┌─────▼─────────────▼───────────▼─────────┐
              │           run_single.py                  │
              │           run_compare.py                 │
              │           visualize.py                   │
              │           gap_sim_config.py (GUI)        │
              └─────────────────────────────────────────┘
```


## 3. 파일 구조

| 파일 | 줄 수 | 역할 |
|------|------|------|
| `common.py` | 1045 | 공통 기반: 스펙, 기하, 물성, 열전달 엔진, j/f 상관식 |
| `module_a.py` | 455 | Gap 열유동: 재순환, 복사, 전도, 혼합 손실 |
| `module_b.py` | 553 | 비말동반: Weber, Monte Carlo, Flux, Risk |
| `module_c.py` | 1382 | 압력강하: f-factor, ΔP_core, ΔP_gap, 습면 보정 |
| `visualize.py` | 1143 | 가시화: 2×3 패널, 개략도, 다중 케이스 |
| `run_single.py` | 490 | 단일 케이스 실행 + CSV 출력 |
| `run_compare.py` | 676 | 다중 케이스 비교 실행 |
| `gap_sim_config.py` | — | tkinter GUI 설정 도구 |


## 4. 공통 기반 (common.py)

### 4.1 데이터 클래스

**FinTubeSpec**: Fin-Tube 열교환기 형상
- 기본 기하: W, H, D, fin_pitch(또는 fpi), fin_thickness, tube_do/di
- 튜브 배열: tube_rows, tube_cols, tube_pitch_t/l, tube_layout
- 핀 타입: fin_type (plain/wavy/slit/louvered) + 각 타입별 파라미터
- 재질: fin_material (Al/Cu), tube_material

**MCHXSpec**: Micro-Channel 열교환기 형상
- 채널: ch_width, ch_height, ch_wall, n_ports
- 루버 핀: fin_pitch, louver_pitch, louver_angle

**RefrigerantState**: 냉매 상태 (CoolProp 기반)
- T_sat_evap/cond, P_evap/cond
- h_fg, ρ_l/v, μ_l, k_l, Pr_l (증발/응축 모두)

### 4.2 기하 계산

**compute_ft_geometry(spec)**
- Plate fin 면적: `A_fin = Nf × 2 × (W×D - Nr×Nt×π×Do²/4)`
- 튜브 노출면적: `A_tube = π×Do×W×(1-δ/Fp)×Nr×Nt`
- 핀효율: Schmidt 등가반경 (`Xm=Pt/2`, `XL=√((Pt/2)²+Pl²)/2`)
- 최소 유로비: `σ = (Fp-δ)/Fp × (Pt-Do)/Pt`

### 4.3 j-factor 상관식 (4핀 × 2배열)

| 핀 타입 | 상관식 | 출처 |
|---------|--------|------|
| Plain | Wang, Chi & Chang (2000) | IJHMT 43(15):2693 |
| Wavy | Wang et al. (1997) | IJHMT 40(4):773 |
| Louvered | Wang, Lee & Chang (1999) | IJHMT 42(1):1 |
| Slit | Wang et al. (2001) / Nakayama & Xu | — |

Inline 배열: `j_inline = j_stag × F_layout` (Nr=2: 0.82, Nr≥3: 0.75)

### 4.4 UA 계산

**compute_UA(spec, geo, ref, V_face, side)**

```
1/UA = 1/(η_o × h_o × A_total) + R_wall + 1/(h_i × A_i)
       ├── 공기측 (지배적 ~61%)    ├── 벽면     ├── 냉매측
```

- h_o: `j × G_c × cp_m / Pr^(2/3)` (핀 타입별 j 자동 분기)
- h_i: Shah (1982) 유동비등 / Shah (1979) 응축
- η_o: Schmidt 등가반경 + h_o 기반 반복 수렴 (5회)

### 4.5 Threlkeld 열전달 엔진 (핵심)

**compute_coil_performance(spec, geo, ua_result, T_in, RH_in, T_wall, V_face)**

습면/건면 자동 전환. 모든 모듈(A/B)이 공유하는 단일 엔진.

```
                T_wall < T_dp ?
                  yes / no
                 /       \
           습면 (Threlkeld)    건면
           ξ = h_fg×dWs/dT/cp  NTU = UA/(m_dot×cp)
           m_wet = m_dry×√(1+ξ)  ε = 1-exp(-NTU)
           η_o,wet → UA_wet     Q = ε×C×(T_in-T_wall)
           NTU = UA_wet/(m_dot×cp)
           ε = 1-exp(-NTU)
           Q = ε×m_dot×(h_in-h_sat,w)  ← 엔탈피 기반
```

**출력**: Q_total, Q_sen, Q_lat, SHR, T_out, W_out, q_cond, ε, NTU

### 4.6 Tube-Row Segmented 모델

**compute_coil_performance_segmented(...)**

Nr개 Row를 순차 계산. 각 Row 출구가 다음 Row 입구.

```
공기 → [Row 1: coil_perf] → T_mid, W_mid, RH_mid → [Row 2: coil_perf] → T_out, W_out
```

Row별 분할: A_total/Nr, A_i/Nr, R_wall×Nr, UA/Nr (h_o, h_i, σ 동일)

효과: 부분 습면 전환 (앞열 습면, 뒷열 건면) 포착. Nr=1이면 Lumped와 동일.


## 5. Module A — Gap 열유동 (module_a.py)

### 5.1 역할

Gap 간격에 따른 **4종 열손실**을 계산하여 순냉방능력(Q_net)과 Cap Retention 산출.

### 5.2 Gap 파라미터

**GapParams**: T_amb, RH_in, CMM, gap_mode (open/semi/sealed), seal_fraction, frame_material

### 5.3 열손실 모델

**① 외부 재순환 (ΔT_recir)**
- 응축기 배출 고온 공기가 Gap을 통해 증발기 흡입구로 역류
- `ΔT_recir = f(gap_mm, gap_mode, seal_fraction, T_cond, T_amb)`
- open: 최대 재순환, sealed: 재순환 없음, semi: seal_fraction으로 보간

**② 내부 혼합 (Q_mix)**
- Gap 내부에서 증발기 냉각 공기와 주변 공기의 혼합
- `η_mix = f(gap_mm, gap_mode)`
- `Q_mix = η_mix × Q_useful`

**③ 복사 열침입 (Q_rad)**
- 고온 응축기 표면 → 저온 증발기 표면 복사
- `Q_rad = σ × F_view × (T_cond⁴ - T_evap⁴) × A`
- 형상인수 F_view: 대향 평판 근사

**④ 전도 단락 (Q_cond)**
- 프레임을 통한 고온→저온 열전도
- `Q_cond = k_frame × A_frame × ΔT / gap`
- 프레임 재질: Al (237 W/mK) / Steel (45 W/mK)

### 5.4 성능 지표

```
Q_net = Q_useful - Q_rad - Q_cond - Q_mix
Cap [%] = Q_net / Q_ref × 100
Cap* [%] = (Q_net - Q_carry) / Q_ref × 100  (비말동반 보정)
```

Q_ref: Gap→∞ (손실 없음) 기준 능력


## 6. Module B — 비말동반 (module_b.py)

### 6.1 역할

증발기 표면 응결수가 공기 흐름에 의해 이탈→응축기에 도달하는 현상 분석.

### 6.2 물리 모델

**응결수 생성**: 공통 엔진 `compute_coil_performance_segmented`의 q_cond 사용

**Weber 수 기반 이탈 판정**:
```
We_ch = ρ_air × V_ch² × Fp / σ_water
We_crit = 0.023 / √Fp
V_onset: We_ch = We_crit를 만족하는 V
V < V_onset → η_co = 0 (SAFE)
V > V_onset → η_co > 0 (ACTIVE)
```

**Monte Carlo 궤적 계산**:
- 50개 액적 발사 (bimodal 크기 분포)
- 항력 + 중력 궤적 적분
- Gap 통과 확률: P_reach = (도달 액적 수) / (전체 액적 수)

**Carryover Flux**:
```
Flux = q_cond × η_co × P_reach / A_face × 1000
```

### 6.3 Risk 등급 (4단계)

| Flux 범위 | 등급 | 의미 |
|-----------|------|------|
| < 1.0 | SAFE | 비말동반 무시 가능 |
| 1.0 ~ 4.0 | CAUTION | 주의 필요 |
| 4.0 ~ 12.0 | DANGER | 응축기 성능 저하 |
| ≥ 12.0 | SEVERE | 심각한 성능 저하 |

### 6.4 Q_carry 응축기 페널티

```
Q_carry = q_cond × η_co × P_reach × h_fg
Cap* = Cap - Q_carry / Q_ref × 100
```


## 7. Module C — 압력강하 (module_c.py)

### 7.1 역할

증발기, Gap, 응축기 각 구간의 압력강하를 계산하여 시스템 ΔP 산출.

### 7.2 f-factor 상관식

4핀 × 2배열 조합별 f-factor (Wang 계열 상관식)

습면 보정: `f_wet = f_dry × F_wet(RH, Re)` — f_wet/f_dry ≈ 1.5~2.3

### 7.3 ΔP 분해

```
ΔP_system = ΔP_evap + ΔP_gap + ΔP_cond
```

**ΔP_core** (Kays & London):
```
ΔP = G²/(2ρ_in) × [Kc + 1 - σ² + f×(A_total/Ac)×(ρ_in/ρ_m) - (1-σ²-Ke)×(ρ_in/ρ_out)]
```

**ΔP_gap**: 급확대 → 자유공간 → 급축소 (등가 길이법)


## 8. 실행 파이프라인

### 8.1 단일 케이스 (run_single.py)

```
JSON 설정 → build_spec → compute_geometry → compute_UA
    → [A] sweep(gaps): simulate_gap × N_gaps
    → [B] analyze_combined: q_cond(공통엔진) + Weber + MC
    → [A+B] compute_carry_penalty → apply_carry_penalty
    → [C] sweep_dp(gaps): f_factor + dp_core + dp_gap
    → 가시화 (module_a/b/c.png + gap_schematic.png)
    → CSV 출력 (4파일)
```

**출력**: module_a/b/c.png, gap_schematic.png, 4개 CSV

### 8.2 다중 비교 (run_compare.py)

```
JSON × N개 → analyze_case() × N → compare_fig_a/b/c()
    → compare_results.csv
```

**사용법**: `python run_compare.py case1.json case2.json ...`


## 9. 가시화 구성

### 9.1 Module A (2×3)

| 위치 | 내용 |
|------|------|
| (0,0) | Cap Retention + Cap* (90% 도달 Gap 표시) |
| (0,1) | Q 분해 (Sensible/Latent + Q_net) |
| (0,2) | 손실 분해 (Q_rad + Q_cond + Q_mix + Q_carry) |
| (1,0) | ΔT_recir |
| (1,1) | SHR (0.75 기준선) |
| (1,2) | 제습량 [g/h] |

### 9.2 Module B (2×3)

| 위치 | 내용 |
|------|------|
| (0,0) | Cap vs Flux (dual axis, Sweet-spot G*) |
| (0,1) | Carryover Flux vs Gap (기준선 표시) |
| (0,2) | P_reach (액적 도달 확률) |
| (1,0) | We/We_crit vs V_face (운전점 마커) |
| (1,1) | Risk Zone Map (4단계 색상) |
| (1,2) | Q_carry 응축기 페널티 |

### 9.3 Module C (2×3)

| 위치 | 내용 |
|------|------|
| (0,0) | ΔP_total vs Gap |
| (0,1) | ΔP stacked (evap/gap/cond) |
| (0,2) | ΔP fraction [%] |
| (1,0) | ΔP_gap 상세 + 비율 |
| (1,1) | f_dry vs f_wet bar |
| (1,2) | 수치 요약 테이블 |

### 9.4 Gap 물리 현상 개략도

증발기-응축기 공간 개략도 + 5대 물리 현상 표시 (재순환, 공기흐름, 복사, 전도, 비말동반 궤적)


## 10. 검증 결과 요약

| 항목 | 등급 | 핵심 수치 |
|------|------|----------|
| 모듈 일관성 (A=B) | A | 0W 차이 (8조건 전수) |
| 기하 계산 | A | 수기검증 일치 |
| h_o / η_o / UA | A | 문헌 범위 내, R_o 지배적 61% |
| ε-NTU (Threlkeld) | A | NTU 0.6~1.8 |
| Gap 손실 | B+ | 전도 35% > 혼합 63% > 복사 2% |
| 비말동반 | B+ | Weber 기반 4단계 점진 전환 |
| 압력강하 | A | f_wet/f_dry 1.97 (문헌 1.7~2.3) |
| Tube-Row Segment | B+ | 부분 습면 포착, Nr=1 fallback |
| 에너지 밸런스 | A | Q_sen+Q_lat = Q_total (0.00%) |

**종합: A-**


## 11. 잔존 한계

| 항목 | 현재 상태 | 필요 사항 |
|------|----------|----------|
| 냉매 사이클 | T_sat 고정 | 압축기 모델 + 사이클 솔버 |
| 팬 연동 | CMM 고정 | 팬 P-Q 곡선 + 반복 수렴 |
| Row별 T_wall | 동일 T_sat | 냉매 회로 + tube-by-tube |
| 공기 분배 | 균일 가정 | CFD 연동 |
| Re > 5000 | 미검증 | 상관식 외삽 주의 |
| 직접 상관식 | 미구현 (etype fallback) | corr_mode='direct' 개발 |


## 12. JSON 설정 파일 구조

```json
{
  "label": "Case Name",
  "evap": {
    "W": 0.30, "H": 0.25, "D": 0.045,
    "fpi": 14, "fin_thickness": 0.10,
    "tube_do": 9.52, "tube_di": 8.52,
    "tube_rows": 2, "tube_cols": 8,
    "tube_pitch_t": 25.4, "tube_pitch_l": 22.0,
    "fin_type": "wavy", "tube_layout": "staggered"
  },
  "cond": { "..." },
  "ref": {
    "refrigerant": "R410A",
    "T_sat_evap": 5.0, "T_sat_cond": 45.0
  },
  "gap": {
    "T_amb": 27.0, "RH_in": 0.70, "CMM": 13.5,
    "gap_mode": "semi", "seal_fraction": 0.7
  },
  "sim": {
    "gap_min": 5, "gap_max": 100, "gap_points": 30
  }
}
```


## 13. 사용법

```bash
# 단일 케이스
python run_single.py gap_config.json

# 다중 비교
python run_compare.py case1.json case2.json case3.json

# GUI 설정 도구
python gap_sim_config.py
```
