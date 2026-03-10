# 열전달 상관식 정리

## 목차
1. 공기측 j-factor (FT)
2. 공기측 j-factor (MCHX)
3. 공기측 f-factor
4. 냉매측 HTC — 증발 (FT: Chen 1966)
5. 냉매측 HTC — 증발 (MCHX: Kim & Mudawar 2013)
6. 냉매측 HTC — 응축
7. 냉매측 HTC — 단상
8. 핀효율 (건면/습면)
9. 코일 성능 모델 (Level 2)
10. 상관식 자동 선택 엔진
11. 검증 결과 요약

---

## 1. 공기측 j-factor (FT)

### 1-1. Plain — Wang, Chi & Chang (2000)
**출처**: IJHMT 43(15):2693-2700, Table 3
**자동 선택**: FT staggered, Pt/Pl ≤ 1.35

### 1-2. Wavy — Wang et al. (1999)
**출처**: IJHMT 42:1945-1956
**자동 선택**: FT staggered, fin_type='wavy'

### 1-3. Louvered — Wang, Lee & Chang (1999)
**출처**: IJHMT 42(1):1-17
**자동 선택**: FT staggered, fin_type='louvered'

### 1-4. Slit — E-type 또는 Wang (2001)
**Direct**: Wang et al. (2001) IJHMT 44:3565 — Dc≥10mm 전용
**E-type**: j_plain × E_slit — 범용
**자동 선택**: Dc≥10mm & Nr≤4 → Direct, 그 외 → E-type

---

## 2. 공기측 j-factor (MCHX)

### Chang & Wang (1997)
**출처**: IJHMT 40(3):533-544
**자동 선택**: MCHX 전용

---

## 3. 공기측 f-factor

### Wang (2000) — Pt/Pl ≤ 1.35
**출처**: IJHMT 43(15):2693, Table 6

### KYW (1999) — Pt/Pl > 1.35
**출처**: Kim, Youn & Webb, ASME JHT 121(3):662

### Enhanced fin 보정
wavy/louver/slit은 f_plain × E_fin(type, Nr)

---

## 4. 냉매측 HTC — 증발 (FT: Chen 1966) ★ 기본값

**출처**: J. Heat Transfer 88(2):189-196

```
h_tp = F × h_l × C_nb

F = 2.35 × (0.213 + 1/Xtt)^0.736    (대류 강화, q'' 무관)
C_nb = 1 + (1-x) × P_r^0.4           (핵비등 보정, q'' 무관)
```

**h_l**: Dittus-Boelter (액상 유량 기준)
**Xtt**: Martinelli parameter

**핵심 설계 선택**:
- 원논문의 S×h_nb(Forster-Zuber) 대신 P_r 기반 C_nb 사용
- 이유: Forster-Zuber의 q''^0.75 항이 자기강화 루프 유발
  (높은 q'' → 높은 h_nb → 낮은 R_i → 낮은 T_wall → 높은 Δh → 높은 Q → 높은 q'' → ...)
- C_nb = 1+(1-x)×P_r^0.4: q'' 무관, CD 검증 h_i ±10% 정합

**검증** (R410A, D=4.42mm, G=298 kg/m²s):
```
x=0.3: Chen=4,170 CD=4,618 (0.90×)
x=0.5: Chen=4,827 CD=5,038 (0.96×)
x=0.7: Chen=5,376 CD=5,398 (1.00×)
```

**자동 선택**: D ≥ 3mm (evap_corr='auto' → 'chen')

### 대안: Gungor-Winterton (1986), Shah (1982)
evap_corr='gungor_winterton' 또는 'shah'로 수동 선택 가능.
단, 핵비등 q'' 의존으로 h_i 과대 가능 (이 조건에서 CD의 2×).

---

## 5. 냉매측 HTC — 증발 (MCHX: Kim & Mudawar 2013)

**출처**: IJHMT 58:718-734 (2013), 10,805 data points

```
h_tp = √(h_nb² + h_cb²)

h_nb = 2345 × (Bo×P_H/P_F)^0.7 × P_r^0.38 × (1-x)^(-0.51) × h_f
h_cb = [5.2×(Bo×P_H/P_F)^0.08 × We_fo^(-0.54) + 3.5/Xtt^0.94 × (ρ_v/ρ_l)^0.25] × h_f
```

P_H/P_F: 직사각형 채널 3면 가열 비율
유효: D_h = 0.19~6.5mm

**자동 선택**: MCHX (D < 3mm)

---

## 6. 냉매측 HTC — 응축

### FT (D ≥ 3mm): Shah (1979)
**출처**: IJHMT 22:547
```
h_cond = h_lo × (1 + 3.8/Z^0.95)
```

### MCHX (D < 3mm): Kim & Mudawar (2012)
**출처**: IJHMT 55:3246-3261
Annular vs Slug/Bubbly regime (We* 기준 자동 분기)

---

## 7. 냉매측 HTC — 단상

### Gnielinski (1976)
과열 증기 및 과냉 액체 공통:
```
Nu = (f/8)(Re-1000)Pr / [1+12.7√(f/8)(Pr^(2/3)-1)]
```

### 전이 블렌딩
- x = 0.90~1.05: 이상↔과열 선형 블렌딩
- x = -0.05~0: 과냉↔이상 선형 블렌딩

---

## 8. 핀효율 (건면/습면)

### 건면: Schmidt 등가원형핀 (FT)
```
m_dry = √(2h_o/(k_fin×δ))
η_fin = tanh(m×r_i×φ)/(m×r_i×φ)
```

### 습면: T_fin_avg 반복 수렴 ★

**핵심**: b-factor를 핀 평균 온도에서 평가 (물리적 정확)

```
① η_dry → T_fin_avg = T_air - η_dry × (T_air - T_wall)
② b(T_fin_avg) = 1 + h_fg × (dWs/dT)|T_fin / cp_air
③ m_wet = √(2×h_o×b/(k×δ)) → η_fin_wet
④ T_fin_avg = T_air - η_fin_wet × (T_air - T_wall)
⑤ 3회 반복 → 수렴
```

**b 역할 분리**:
- b → m_wet → η_fin_wet (핀효율에만 반영)
- h_o 자체에는 b 미적용
- UA에는 b 미포함 (이중 계산 방지)

**대표값** (T_air=45°C, RH=70%, T_wall=31°C):
```
T_fin_avg = 36.2°C, b = 6.79
η_fin_dry = 0.913, η_fin_wet = 0.630
η_o_dry = 0.920, η_o_wet = 0.660
```

---

## 9. 코일 성능 모델 (Level 2)

### compute_coil_v3

**구조**: Nr×Nt×N_seg 세그먼트 개별 계산

**T_wall 반복 수렴**:
```
T_wall = T_ref + Q × R_i  (α=0.7, 12회 이내)
h_i 매 반복마다 재계산 (Q ↔ Bo ↔ h_i 동시 수렴)
```

**공기 상태 Row간 업데이트** (엔탈피 기반):
```
h_out = h_in - Q_row/m_air
T_out = (h_out - W×2501000)/(1006 + W×1860)
W_out = W_in - Q_lat/(m×h_fg)
```

**냉매 상태 업데이트**:
- 이상(x≤1): dx = Q/(m_ref×h_fg)
- 과열(x>1): dT = Q/(m_ref×cp_v)
- transition(x>1): 과열로 처리 (★ 이전 버그 수정)

---

## 10. 상관식 자동 선택 엔진

`select_correlations(spec, geo, ref, side)` → dict

| 조건 | j | f | h_i (evap) | h_i (cond) |
|------|---|---|-----------|-----------|
| FT plain Pt/Pl≤1.35 | Wang(2000) | Wang(2000) | Chen(1966) | Shah(1979) |
| FT plain Pt/Pl>1.35 | Wang(2000)+⚠ | KYW(1999) | Chen(1966) | Shah(1979) |
| FT wavy | Wang(1999) | — | Chen(1966) | Shah(1979) |
| FT louver | Wang(1999) | — | Chen(1966) | Shah(1979) |
| FT slit Dc≥10 | Wang(2001) | — | Chen(1966) | Shah(1979) |
| MCHX (D<3mm) | C&W(1997) | C&W(1997) | K&M(2013) | K&M(2012) |

---

## 11. 검증 결과 요약

### CoilDesigner 비교 (R410A, Issue#1 스펙)

| 항목 | 시뮬레이터 | CD | 오차 | 등급 |
|------|-----------|-----|------|------|
| A_total | 0.776 m² | 0.788 m² | -1.5% | A |
| h_o | 102 W/m²K | 96.2 | +6% | A |
| dp | 58.9 Pa | 63 Pa | -7% | A |
| SHR | 0.330 | 0.331 | +0.2% | A+ |
| h_i (이상, 평균) | 1.02× CD | 기준 | ±10% | A |
| T_wall | ~31°C | ~33°C | -2°C | A |
| Q_total | 790W | 1372W | -42% | C+ |

**Q_total -42% 원인 분석**:
습면 핀효율 (η_o=0.66) vs CD 건면 핀효율 (η_o=0.93).
CD는 습면 핀효율 감소를 반영하지 않는 간소화 방식 사용.
물리적으로는 우리 모델이 더 정확 (습면 η 감소는 실제 현상).
