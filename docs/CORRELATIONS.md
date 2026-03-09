# 열전달 상관식 정리

## 목차
1. 공기측 j-factor (FT)
2. 공기측 j-factor (MCHX)
3. 공기측 f-factor
4. 냉매측 HTC — 증발
5. 냉매측 HTC — 응축
6. 냉매측 HTC — 단상 (과열/과냉)
7. 핀효율
8. 코일 성능 모델 (Level 0/1/2)
9. 상관식 자동 선택 엔진
10. 검증 결과 요약

---

## 1. 공기측 j-factor (FT)

### 1-1. Plain — Wang, Chi & Chang (2000)

**출처**: IJHMT 43(15):2693-2700, Table 3

**수식 (Nr≥2, Re≥1000)**:
```
j = 0.086 × Re_Dc^p3 × Nr^p4 × (Fp/Dc)^p5 × (Fp/Dh)^p6 × (Fp/Pt)^(-0.93)
```

**Re 정의**: Re_Dc = G×Dc/μ (Dc = Do + 2δf)

**유효 범위**: Re 300~20000, Nr 1~6, Dc 6.35~12.7mm

**자동 선택 조건**: FT staggered, Pt/Pl ≤ 1.35

---

### 1-2. Wavy — Wang et al. (1999)

**출처**: IJHMT 42:1945-1956, Table 2

**자동 선택 조건**: FT staggered, fin_type='wavy'

---

### 1-3. Louvered — Wang, Lee & Chang (1999)

**출처**: IJHMT 42(1):1-17, Table 1

**주의**: θ는 radian 사용 (degree 사용 시 j 발산)

**자동 선택 조건**: FT staggered, fin_type='louvered'

---

### 1-4. Slit — E-type 또는 Wang (2001)

**E-type**: j_slit = j_plain × E_slit (E_slit ≈ 1.33)
- 범용, Dc 제한 없음

**Direct — Wang et al. (2001)**: IJHMT 44:3565, Eq.(1)
- Do ≈ 10mm 전용. Dc<10mm에서 j/j_plain=2.0 (문헌 1.2~1.6 초과)

**자동 선택**: Dc ≥ 10mm & Nr ≤ 4 → Direct, 그 외 → E-type

---

## 2. 공기측 j-factor (MCHX)

### Chang & Wang (1997)

**출처**: IJHMT 40(3):533-544

```
j = Re_Lp^J1 × (θ/90)^0.27 × (Fp/Lp)^(-0.14) × ...
J1 = -0.49 × (θ/90)^0.27
```

**자동 선택**: MCHX 전용

---

## 3. 공기측 f-factor

### 3-1. Wang (2000) — Pt/Pl ≤ 1.35

**출처**: IJHMT 43(15):2693, Table 6

```
f = 0.0267 × Re^f1 × (Pt/Pl)^f2 × (Fp/Dc)^f3
```

**자동 선택 조건**: FT staggered, Pt/Pl ≤ 1.35

### 3-2. KYW (1999) — Pt/Pl > 1.35

**출처**: Kim, Youn & Webb (1999) ASME JHT 121(3):662

```
f = f_tube(Jakob 1938) + f_fin(KYW Eq.6)    # 중첩 모델
```

**자동 선택 조건**: FT staggered, Pt/Pl > 1.35

### 3-3. Enhanced fin 보정

Enhanced 핀(wavy/louver/slit)은 plain f에 E_fin 보정:
```
f = f_plain × E_fin(type, Nr)
wavy:    E = (1.10 + ...) × (Nr/2)^0.10
louver:  E = (1.25 + ...) × (Nr/2)^0.15
slit:    E = (1.25 + ...) × (Nr/2)^0.30
```

---

## 4. 냉매측 HTC — 증발

### 4-1. Gungor & Winterton (1986) — 기본값

**출처**: IJHMT 29(3):351-358

```
h_tp = E × h_l + S × h_pool

E = 1 + 24000 × Bo^1.16 + 1.37 × (1/Xtt)^0.86    # Enhancement
S = 1 / (1 + 1.15e-6 × E² × Re_l^1.17)            # Suppression
h_pool = Cooper(1984) pool boiling correlation
```

**Martinelli parameter**: Xtt = ((1-x)/x)^0.9 × (ρ_v/ρ_l)^0.5 × (μ_l/μ_v)^0.1

**수평 튜브 보정 (1987 update)**: Fr < 0.05 시 E, S 보정

**장점**: x > 0.8에서 안정적 (Shah 대비), 물리 기반 (핵비등+대류비등 중첩)

**자동 선택**: D ≥ 6mm → 기본, D 3~6mm → 경고, D < 3mm → 미니채널 경고

### 4-2. Shah (1982) — 선택 가능

**출처**: ASHRAE Trans 88(1):185

```
h_tp = ψ × h_lo
h_lo = 0.023 × Re_lo^0.8 × Pr_l^0.4 × k_l/D    # 전체 유량 기준 ★
```

**3-regime chart correlation**:
- N > 1.0 (핵비등 지배): ψ = max(ψ_nb, ψ_cb)
- N ≤ 1.0 (대류비등 지배): ψ = max(ψ_bs, ψ_cb)
- ψ_cb = 1.8 / N^0.8

**주의**: h_lo 기준 (h_l((1-x)G) 아님). Bo = q''/(G×h_fg) (세그먼트 열유속 기반)

**문제점**: x > 0.8에서 N → 0 → ψ → ∞ (발산)

**자동 선택**: `evap_corr='shah'`로 수동 지정 가능

### 4-3. Bo (Boiling number) 계산

```
Bo = q'' / (G × h_fg)
q'' = Q_seg / A_i_seg        # 세그먼트당 열유속
A_i_seg = A_i / (Nr × Nt × N_seg)
```

**이전 오류**: A_i/Nr 사용 → Bo 40배 과소 → 핵비등 미반영
**현재**: A_i_seg 사용 → 정확한 Bo

---

## 5. 냉매측 HTC — 응축

### Shah (1979)

**출처**: IJHMT 22:547

```
h_cond = h_lo × (1 + 3.8 / Z^0.95)
Z = (1/x - 1)^0.8 × (P/Pcrit)^0.4
```

---

## 6. 냉매측 HTC — 단상

### Gnielinski (1976)

**과열 증기 및 과냉 액체 공통**:
```
Nu = (f/8)(Re-1000)Pr / [1 + 12.7√(f/8)(Pr^(2/3)-1)]
f = (0.790 ln(Re) - 1.64)^(-2)
```

**유효**: 2300 < Re < 5×10^6, 0.5 < Pr < 2000

### 전이 블렌딩 (x = 0.90 ~ 1.05)

```
w = (x - 0.90) / 0.15    # 0→1 선형 가중
h = (1-w) × h_2phase + w × h_vapor
```

---

## 7. 핀효율

### Schmidt 등가원형핀 (FT)
```
m = √(2h_o/(k_fin×δ_f))
r_eq = 1.27 × Xm × √(XL/Xm - 0.3)
η_fin = tanh(m×r_i×φ) / (m×r_i×φ)
```

### 습면 핀효율
```
m_wet = √(2×h_o×b / (k_fin×δ_f))    # b는 여기에만!
η_fin_wet < η_fin_dry
```

**핵심**: b-factor는 핀효율 계산(m값)에만 반영. h_o 자체는 습면/건면 동일.

---

## 8. 코일 성능 모델

### Level 0 — Row-by-Row (기존)
- Nr개 Row 분할, T_wall = T_sat 고정, h_i = x=0.5 고정
- 함수: `compute_coil_performance_segmented()`

### Level 1 — Row × Phase-zone
- Row별 건도(x) 추적, 이상/과열 자동 분기
- dx 상한, b-factor 제한, ρ(T) 보정
- 함수: `compute_coil_v2()`

### Level 2 — Tube-Segment (CoilDesigner 수준)
- Nr×Nt×N_seg = 160 세그먼트 개별 계산
- T_wall 반복 수렴 (Q_air = Q_ref 동시 만족)
- h_i 매 반복마다 재계산 (Q ↔ Bo ↔ h_i ↔ R_i ↔ T_wall)
- Counter/Parallel flow 선택
- 공기 포화 보정 (과포화 방지)
- 함수: `compute_coil_v3(spec, geo, ref, T_air, RH, V, m_ref, x_in, side, N_seg, flow, evap_corr)`

---

## 9. 상관식 자동 선택 엔진

### `select_correlations(spec, geo, ref, side)`

기하·운전조건으로 최적 상관식을 자동 선택. 사용자가 직접 고를 필요 없음.

**선택 트리**:

| 조건 | 선택 |
|------|------|
| FT staggered plain, Pt/Pl ≤ 1.35 | j: Wang(2000), f: Wang(2000) |
| FT staggered plain, Pt/Pl > 1.35 | j: Wang(2000)+⚠, f: KYW(1999) |
| FT staggered wavy | j: Wang(1999) |
| FT staggered louver | j: Wang(1999) |
| FT staggered slit, Dc ≥ 10mm | j: Wang(2001) direct |
| FT staggered slit, Dc < 10mm | j: Plain × E_slit |
| FT inline | j/f: Wang(2000) inline |
| MCHX | j/f: Chang & Wang(1997) |
| 증발 D ≥ 6mm | h_i: Gungor-Winterton(1986) |
| 증발 D < 3mm | h_i: G-W + ⚠미니채널 경고 |
| 응축 | h_i: Shah(1979) |
| 과열/과냉 | Gnielinski(1976) |

**범위 밖 경고**: Pt/Pl>1.35, Dc<10mm(slit), D<3mm(미니채널) 등

---

## 10. 검증 결과 요약

### CoilDesigner 비교 (Issue#1 스펙)

**기하/물성 검증**:

| 항목 | 시뮬레이터 | CoilDesigner | 오차 |
|------|-----------|-------------|------|
| A_total (증발기) | 0.776 m² | 0.788 m² | -1.5% ✅ |
| h_o (증발기) | 102 | 96 W/m²K | +6% ✅ |
| dp (증발기) | 58.9 Pa | 63 Pa | -6.5% ✅ |

**코일 성능 검증 (Level 2, counter, G-W)**:

| 항목 | Level 0 | Level 2 | CD | Level 2/CD |
|------|---------|---------|-----|-----------|
| Q_total | 4,306W | 952W | 1,372W | 0.69 |

**개선 이력**:
```
Level 0:        Q=4,306W (3.14×)  T_wall 고정, h_i 고정
Level 2+T_wall: Q=  952W (0.69×)  160 seg, T_wall 수렴, G-W
Level 2+Shah:   Q=1,067W (0.78×)  Shah(1982), T_wall 수렴
CD:             Q=1,372W (1.00×)  목표
```
