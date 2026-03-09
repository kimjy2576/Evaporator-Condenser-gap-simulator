# 열전달 상관식 정리

## 목차
1. 공기측 j-factor (FT)
2. 공기측 j-factor (MCHX)
3. 냉매측 HTC — 증발
4. 냉매측 HTC — 응축
5. 핀효율
6. 코일 성능 모델 (ε-NTU)
7. 검증 결과 요약
8. 개선 필요 사항

---

## 1. 공기측 j-factor (FT)

### 1-1. Plain — Wang, Chi & Chang (2000)

**출처**: IJHMT 43(15):2693-2700, Table 3

**방식**: Direct

**수식 (Nr≥2, Re≥1000)**:
```
j = 0.086 × Re_Dc^p3 × Nr^p4 × (Fp/Dc)^p5 × (Fp/Dh)^p6 × (Fp/Pt)^(-0.93)

p3 = -0.361 - 0.042×Nr/ln(Re) + 0.158×ln(Nr×(Fp/Dc)^0.41)
p4 = -1.224 - 0.076×(Pl/Dh)^1.42 / ln(Re)
p5 = -0.083 + 0.058×Nr/ln(Re)
p6 = -5.735 + 1.21×ln(Re/Nr)
```

**Re 정의**: Re_Dc = G×Dc/μ (Dc = Do + 2δf, collar diameter)

**유효 범위**:
| 변수 | 범위 | 비고 |
|------|------|------|
| Re_Dc | 300 ~ 20,000 | |
| Nr | 1 ~ 6 | Nr=1은 별도 수식 |
| Dc | 6.35 ~ 12.7 mm | |
| Fp | 1.19 ~ 8.7 mm | |
| Pt | 17.7 ~ 31.75 mm | |
| Pl | 12.4 ~ 27.5 mm | |

**검증 정확도**: 평균 편차 7.5% (74개 샘플)

**현재 h_o @ V=1.5, 기본 스펙**: 61 W/m²K (j=0.016)

---

### 1-2. Wavy — Wang et al. (1999)

**출처**: IJHMT 42:1945-1956, Table 2

**방식**: Direct

**수식 (Nr≥2)**:
```
j = 0.394 × Re_Dc^j1 × (Fp/Dc)^(-0.031) × (Pd/Dc)^0.346 × Nr^(-0.107)

j1 = -0.329 - 0.147×ln(Nr×(Fp/Dc)^0.44 × (Pd/Dc)^0.09 × (Pt/Pl)^(-0.04))
```

**핵심 변수**: Pd = wavy height (amplitude, 반파고)

**유효 범위**: Re_Dc 300~10,000, Pd 1.18~1.58 mm

**현재 h_o**: 87 W/m²K (j/j_plain = 1.40)

**문헌 비율**: j_wavy/j_plain = 1.1~1.4 → ✅

---

### 1-3. Louvered — Wang, Lee & Chang (1999)

**출처**: IJHMT 42(1):1-8, Table 1

**방식**: Direct (Re_Lp 기반, θ radian)

**수식 (Nr≥2)**:
```
j = 0.394 × Re_Lp^j2 × (Fp/Pl)^(-0.39) × (Td/Pt)^0.46 × (Lp/Fl)^0.33

j2 = -0.545 + 0.0538×θ_rad - 0.0244×Nr

Re_Lp = G × Lp / μ  (Re_Dc → Re_Lp 변환: Re_Lp = Re_Dc × Lp/Dc)
θ_rad = θ_deg × π/180
```

**주의**: 원논문에서 θ 단위가 모호 — degree 사용 시 j 폭발. radian 사용 시 검증됨.

**유효 범위**: Re_Lp 100~5,000, θ 15~35°

**현재 h_o**: 92 W/m²K (j/j_plain = 1.48)

**문헌 비율**: j_louver/j_plain = 1.3~1.8 → ✅

---

### 1-4. Slit — E-type (문헌 비율)

**방식**: j_plain × E_slit

**수식**:
```
j_slit = j_plain × E_s

E_s = 1.15 + 0.12 × min(Ns,10)^0.30 × (Sh/Fp)^0.20
E_s ∈ [1.20, 1.50]
```

**근거**: Yan & Sheen (2000), Wang (2001) — direct 상관식은 Fp/Dc 범위 불안정으로 사용 불가

**현재 h_o**: 81 W/m²K (j/j_plain = 1.33)

**문헌 비율**: j_slit/j_plain = 1.2~1.6 → ✅

---

## 2. 공기측 j-factor (MCHX)

### Chang & Wang (1997)

**출처**: 간략화 상관식

**수식**:
```
j = 0.5 × Re_Lp^(-0.49) × (α/90)^0.27 × (Fp/Lp)^(-0.14)

Re_Lp = G × Lp / μ
α = louver angle [rad]
```

**유효 범위**: Re_Lp 100~1,000

**현재 h_o**: 921 W/m²K (j=0.31)

**⚠️ 문제점**:
- 간략화 상관식 — 정확도 낮음
- 더 정확한 대안: Chang & Wang (2006) 일반화 상관식
- 또는 Kim & Bullard (2002)

---

## 3. 냉매측 HTC — 증발

### Shah (1982) 유동비등

**출처**: Shah, M.M. (1982) ASHRAE Trans. 88(1):185-196

**수식**:
```
h_evap = h_l × ψ

h_l = 0.023 × Re_l^0.8 × Pr_l^0.4 × k_l / D  (Dittus-Boelter)
Re_l = G_ref × (1-x) × D / μ_l

Co = ((1-x)/x)^0.8 × (ρ_v/ρ_l)^0.5     (Convection number)
Bo = Q / (A_i × G_ref × h_fg)            (Boiling number)

Co > 0.65:  ψ = 1.8 / Co^0.8
Co ≤ 0.65:  ψ = max(1.8/Co^0.8, 0.6683×Co^(-0.2) + 1058×Bo^0.7)
```

**현재 설정**: x_mean = 0.5 (고정), G_ref = Q/(h_fg × A_cs)

**⚠️ 문제점**:
- x = 0.5 고정 → 실제는 입구~출구 건도 변화
- 과열 영역 미반영 (CoilDesigner: 59% 과열)
- G_ref 추정이 Q에 의존 → 순환 참조 가능성

**현재 h_i**: 2,888 W/m²K (FT), 16,215 W/m²K (MCHX)

---

## 4. 냉매측 HTC — 응축

### Shah (1979) 응축

**출처**: Shah, M.M. (1979) IJHMT 22:547-556

**수식**:
```
h_cond = h_l × (1 + 3.8 / Z^0.95)

h_l = 0.023 × Re_l^0.8 × Pr_l^0.4 × k_l / D  (Dittus-Boelter)
Z = (1/x - 1)^0.8 × (P/P_crit)^0.4
```

**현재 설정**: x_mean = 0.5 (고정)

**⚠️ 문제점**: 증발기와 동일 (x 고정, 과냉 미반영)

**현재 h_i**: 1,446 W/m²K (FT), 2,806 W/m²K (MCHX)

---

## 5. 핀효율

### Schmidt 등가 반경법 (FT)

**수식**:
```
r_eq = 1.27 × Xm × √(XL/Xm - 0.3)
Xm = Pt/2
XL = √((Pt/2)² + Pl²) / 2

φ = (r_eq/r_i - 1) × (1 + 0.35×ln(r_eq/r_i))
m = √(2h_o / (k_fin × δf))

η_fin = tanh(m×r_i×φ) / (m×r_i×φ)
η_o = 1 - (A_fin/A_total) × (1 - η_fin)
```

**반복 수렴**: h_o → η_o → UA → h_o (최대 5회)

### 직선 핀 근사 (MCHX)

```
L_fin = slab_pitch / 2
η_fin = tanh(m×L_fin) / (m×L_fin)
```

---

## 6. 코일 성능 모델 (Threlkeld ε-NTU)

### 구조

```
Tube-Row Segmented: Nr개 Row 순차 계산
각 Row: Threlkeld 습면/건면 ε-NTU

1) T_dp > T_wall → 습면
   UA_wet = f(b, η_o_wet, h_o, A, R_wall, R_i)
   b = cp_s / cp_a = (cp_a + h_fg × dWs/dT) / cp_a
   ε_wet = 1 - exp(-NTU_wet)
   Q = ε × m_dot × (h_in - h_sat_w)

2) T_dp ≤ T_wall → 건면
   ε_dry = 1 - exp(-NTU_dry)
   Q = ε × m_dot × cp_a × (T_in - T_wall)
```

### ⚠️ 핵심 가정

| 가정 | 영향 | 개선 방법 |
|------|------|-----------|
| **T_wall = T_sat_evap (고정)** | Q 과대 (3×) | 냉매 순환 모델 |
| **x_mean = 0.5 (고정)** | h_i 과대 | 세그먼트별 x 추적 |
| **ρ_air = 1.18 (20°C 고정)** | m_dot ±6% | T 기반 밀도 |
| **Row 단위 분할** | 부분 습면 포착 | Tube-Segment 분할 |

---

## 7. 검증 결과 요약

### CoilDesigner 비교 (Issue#1 스펙)

| 항목 | 시뮬레이터 | CoilDesigner | 오차 |
|------|-----------|-------------|------|
| A_total (증발기) | 0.776 m² | 0.788 m² | **-1.5%** ✅ |
| A_total (응축기) | 1.272 m² | 1.293 m² | **-1.6%** ✅ |
| h_o (증발기) | 102 | 96 W/m²K | **+6%** ✅ |
| h_o (응축기) | 128 | 128 W/m²K | **-0.2%** ✅ |
| dp (증발기) | 58.9 Pa | 63 Pa | **-6.5%** ✅ |
| dp (응축기) | 211 Pa | 223 Pa | **-5.1%** ✅ |
| Q_total | 4,396 W | 1,372 W | **3.2×** ❌ |

### 문헌 j/f 비율 검증

| 핀 | j/j_plain | 문헌 | f/f_plain | 문헌 |
|----|-----------|------|-----------|------|
| Plain | 1.00 | (기준) | 1.00 | (기준) |
| Wavy | 1.40 | 1.1~1.4 ✅ | 1.30 | 1.2~1.5 ✅ |
| Louver | 1.48 | 1.3~1.8 ✅ | 1.73 | 1.5~2.5 ✅ |
| Slit | 1.33 | 1.2~1.6 ✅ | 1.49 | 1.3~1.7 ✅ |

---

## 8. 개선 필요 사항

### 우선순위 높음

| # | 항목 | 현재 | 목표 | 영향 |
|---|------|------|------|------|
| 1 | **공기 밀도** | ρ=1.18 고정 | T 기반 계산 | m_dot ±6%, dp ±12% |
| 2 | **MCHX j-factor** | Chang&Wang(1997) 간략화 | Chang&Wang(2006) 일반화 또는 Kim&Bullard(2002) | h_o 정확도 |
| 3 | **h_i x 의존성** | x=0.5 고정 | Row별 x 추적 | h_i 정확도 |
| 4 | **습면 f-factor** | 건면 f 사용 | 습면 보정 f_wet/f_dry | dp 정확도 |

### 우선순위 중간

| # | 항목 | 현재 | 목표 |
|---|------|------|------|
| 5 | FT inline j-factor | 간략화 | Wang(2000) 원논문 |
| 6 | 핀효율 습면 보정 | Threlkeld 내부 | 독립 함수화 |
| 7 | 접촉 열저항 | 미반영 | R_contact 추가 |

### 우선순위 낮음 (구조 변경 필요)

| # | 항목 | 현재 | 목표 |
|---|------|------|------|
| 8 | T_wall 고정 해소 | T_sat_evap | 냉매 순환 연동 |
| 9 | Tube-Segment 분할 | Row 단위 | Tube×10 세그먼트 |
| 10 | 과열/과냉 영역 | 미반영 | 상 분기 h_i |
