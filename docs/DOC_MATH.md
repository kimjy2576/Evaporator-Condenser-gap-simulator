# 📐 적용 상관식 및 수학적 모델 설명서

## 1. 공통 열전달 엔진 (Threlkeld ε-NTU)

### 1.1 건면 UA 계산

3저항 직렬 모델:

$$\frac{1}{UA_{dry}} = \underbrace{\frac{1}{\eta_o \cdot h_o \cdot A_{total}}}_{R_o\;(공기측)} + \underbrace{\frac{\ln(D_o/D_i)}{2\pi k_{tube} L N_t}}_{R_{wall}} + \underbrace{\frac{1}{h_i \cdot A_i}}_{R_i\;(냉매측)}$$

### 1.2 공기측 HTC (j-factor)

$$h_o = j \cdot G_c \cdot c_{p,m} / Pr^{2/3}$$

여기서 $G_c = \rho_{air} \cdot V_{face} / \sigma$, $Re_{D_c} = G_c \cdot D_c / \mu_{air}$

### 1.3 핀 효율 (Schmidt 등가반경법)

$$m_{dry} = \sqrt{\frac{2 h_o}{k_{fin} \cdot \delta_f}}$$

Staggered 배열 등가반경:

$$X_M = P_t/2, \quad X_L = \frac{\sqrt{(P_t/2)^2 + P_l^2}}{2}$$

$$r_{eq} = 1.27 \cdot X_M \cdot \sqrt{X_L/X_M - 0.3}$$

$$\phi = (r_{eq}/r_i - 1)(1 + 0.35 \ln(r_{eq}/r_i))$$

$$\eta_{fin} = \frac{\tanh(m_{dry} \cdot r_i \cdot \phi)}{m_{dry} \cdot r_i \cdot \phi}, \quad \eta_o = 1 - \frac{A_{fin}}{A_{total}}(1 - \eta_{fin})$$

### 1.4 습면 ε-NTU (Threlkeld 정석)

습면 판정: $T_{wall} < T_{dp}$

습면 보정 계수:

$$\xi = \frac{h_{fg} \cdot (dW_s/dT)|_{T_m}}{c_{p,a}}, \quad T_m = \frac{T_{in} + T_{wall}}{2}$$

습면 핀효율:

$$m_{wet} = m_{dry} \cdot \sqrt{1 + \xi}$$

$$\eta_{o,wet} = 1 - \frac{A_{fin}}{A_{total}}(1 - \eta_{fin,wet})$$

습면 유효비열 및 UA:

$$c_{p,s} = c_{p,a} + h_{fg} \cdot (dW_s/dT)|_{T_{wall}}, \quad b = c_{p,s}/c_{p,a}$$

$$UA_{o,wet} = b \cdot \eta_{o,wet} \cdot h_o \cdot A_{total}$$

$$\frac{1}{UA_{wet}} = \frac{1}{UA_{o,wet}} + R_{wall} + R_i$$

ε-NTU:

$$NTU_{wet} = \frac{UA_{wet}}{\dot{m} \cdot c_{p,a}}, \quad \varepsilon_{wet} = 1 - e^{-NTU_{wet}}$$

Q 계산 (엔탈피 기반):

$$Q_{total} = \varepsilon_{wet} \cdot \dot{m} \cdot (h_{in} - h_{sat,wall})$$

$$W_{out} = W_{in} - \varepsilon_{wet} \cdot (W_{in} - W_{sat,wall})$$

$$Q_{lat} = \dot{m} \cdot (W_{in} - W_{out}) \cdot h_{fg}, \quad Q_{sen} = Q_{total} - Q_{lat}$$

### 1.5 Tube-Row Segment 모델

Nr개 Row 순차 계산. Row당 분할:

$$A_{row} = A_{total}/N_r, \quad UA_{row} = UA/N_r, \quad R_{row} = R \cdot N_r$$

각 Row 출구 → 다음 Row 입구:

$$T_{air}^{(i+1)} = T_{out}^{(i)}, \quad W_{air}^{(i+1)} = W_{out}^{(i)}$$

RH 역산: $RH^{(i+1)} = f(T_{air}^{(i+1)}, W_{air}^{(i+1)})$ → 습면/건면 재판정

---

## 2. j-factor 상관식 (4핀 × 2배열)

### 2.1 Plain — Wang, Chi & Chang (2000)

출처: IJHMT 43(15):2693-2700

$$j = 0.108 \cdot Re_{D_c}^{-0.29} \cdot \left(\frac{P_t}{P_l}\right)^{c_1} \cdot \left(\frac{F_p}{D_c}\right)^{c_2} \cdot N_r^{c_3}$$

적용 범위: 300 ≤ Re ≤ 5000, 1 ≤ Nr ≤ 6

### 2.2 Wavy — Wang et al. (1997)

출처: IJHMT 40(4):773-784

$$j = j_{plain} \cdot E_{wavy}$$

$$E_{wavy} = f(Re, P_d/D_c, X_f/D_c, \alpha_{wavy})$$

P_d: 파형 높이, X_f: 파형 피치, α: 파형 각도

적용 범위: 300 ≤ Re ≤ 7500

### 2.3 Louvered — Wang, Lee & Chang (1999)

출처: IJHMT 42(1):1-17

$$j = j_{plain} \cdot E_{louver}$$

$$E_{louver} = f(Re, L_p/F_p, \theta_l, N_r)$$

L_p: 루버 피치, θ_l: 루버 각도

적용 범위: 100 ≤ Re ≤ 3000 (Re_Lp 기준)

### 2.4 Slit — Wang et al. (2001)

$$j = j_{plain} \cdot E_{slit}$$

$$E_{slit} = f(Re, S_n, S_h/F_p)$$

S_n: 슬릿 수, S_h: 슬릿 높이

### 2.5 Inline 배열 보정

$$j_{inline} = j_{staggered} \times F_{layout}$$

| Nr | F_layout |
|----|----------|
| 1 | 1.0 |
| 2 | 0.82 |
| ≥3 | 0.75 |

물리: Wake shielding → 후열 h_o 감소

---

## 3. 냉매측 HTC

### 3.1 증발 — Chen (1966) [FT 기본값] ★

$$h_{tp} = F \cdot h_l \cdot C_{nb}$$

**F factor** (대류 강화, q'' 무관):

$$F = \begin{cases} 1.0 & 1/X_{tt} \leq 0.1 \\ 2.35 \cdot (0.213 + 1/X_{tt})^{0.736} & 1/X_{tt} > 0.1 \end{cases}$$

**C_nb** (핵비등 보정, q'' 무관): $C_{nb} = 1 + (1-x) \cdot P_r^{0.4}$

**Martinelli**: $X_{tt} = \left(\frac{1-x}{x}\right)^{0.9} \left(\frac{\rho_v}{\rho_l}\right)^{0.5} \left(\frac{\mu_l}{\mu_v}\right)^{0.1}$

원논문의 S×h_nb(Forster-Zuber) 대신 C_nb 사용. Forster-Zuber의 $q''^{0.75}$ 항이 습면에서 자기강화 루프를 유발하기 때문.

### 3.2 증발 — Kim & Mudawar (2013) [MCHX 기본값] ★

$$h_{tp} = \sqrt{h_{nb}^2 + h_{cb}^2}$$

$h_{nb} = 2345 (Bo \cdot P_H/P_F)^{0.7} P_r^{0.38} (1-x)^{-0.51} h_f$, $h_{cb} = [5.2(Bo \cdot P_H/P_F)^{0.08} We_{fo}^{-0.54} + 3.5/X_{tt}^{0.94} (\rho_v/\rho_l)^{0.25}] h_f$

유효: D_h = 0.19~6.5mm (미니채널)

### 3.3 증발 — Gungor-Winterton (1986), Shah (1982) [선택 가능]

수동 선택: evap_corr='gungor_winterton' 또는 'shah'. 핵비등 q'' 의존.

### 3.4 응축 — Shah (1979) / Kim & Mudawar (2012)

FT (D≥3mm): $h_{cond} = h_{lo} \cdot (1 + 3.8/Z^{0.95})$

MCHX (D<3mm): Kim & Mudawar (2012) Annular/Slug regime 자동 분기

### 3.5 단상 — Gnielinski (1976)

$$Nu = \frac{(f/8)(Re-1000)Pr}{1 + 12.7\sqrt{f/8}(Pr^{2/3}-1)}, \quad f = (0.790 \ln Re - 1.64)^{-2}$$

### 3.6 전이 블렌딩

- x = 0.90~1.05: $h = (1-w) h_{2\phi} + w \cdot h_{vapor}$
- transition(x>1): 과열로 처리 (T_ref 정상 업데이트)

### 3.7 T_wall 반복 수렴 (Level 2)

$$T_w^{(n+1)} = \alpha \cdot (T_{ref} + Q_{air}^{(n)} \cdot R_i) + (1-\alpha) \cdot T_w^{(n)}$$

$\alpha = 0.7$, 12회 이내 수렴 (8회 일반적).

### 3.8 습면 핀효율 — b@T_fin_avg 반복 수렴 ★

$$b(T) = 1 + \frac{h_{fg} \cdot (dW_s/dT)|_T}{c_{p,a}}$$

**반복**: $T_{fin}^{(0)} = T_{air} - \eta_{dry} (T_{air}-T_w)$ → $b(T_{fin})$ → $\eta_{wet}$ → $T_{fin} = T_{air} - \eta_{wet}(T_{air}-T_w)$ → 3회 반복

**b 역할 분리**:
- b → m_wet → η_fin_wet (핀효율에만)
- h_o에는 b 미적용, UA에도 b 미포함

### 3.9 Row간 공기 상태 업데이트 (엔탈피 기반) ★

$$h_{out} = h_{in} - Q_{row}/\dot{m}_a, \quad T_{out} = \frac{h_{out} - W \cdot 2501000}{1006 + W \cdot 1860}$$

---

## 4. f-factor 및 압력강하

### 4.1 Kays & London 코어 ΔP

$$\Delta P = \frac{G^2}{2\rho_{in}} \left[ K_c + 1 - \sigma^2 + f \cdot \frac{A_{total}}{A_c} \cdot \frac{\rho_{in}}{\rho_m} - (1-\sigma^2-K_e)\frac{\rho_{in}}{\rho_{out}} \right]$$

### 4.2 습면 f-factor 보정

$$f_{wet} = f_{dry} \cdot F_{wet}(Re, RH)$$

문헌 비율: Plain 1.5~2.0, Wavy 1.7~2.3, Louvered 1.8~2.5

### 4.3 Gap 압력강하

급확대 + 자유공간 + 급축소:

$$\Delta P_{gap} = \frac{\rho V_{gap}^2}{2} \left[ K_{exp} + f_{gap}\frac{L_{gap}}{D_{h,gap}} + K_{cont} \right]$$

---

## 5. Gap 열손실 모델 (Module A)

### 5.1 외부 재순환

$$\Delta T_{recir} = f(G, \text{mode}, \text{sf}) \cdot (T_{cond} - T_{amb})$$

open: 최대 재순환 (sf=0), sealed: 재순환 없음 (sf=1), semi: 보간

### 5.2 복사 열침입

$$Q_{rad} = \sigma_{SB} \cdot F_{view} \cdot A \cdot (T_{cond}^4 - T_{evap}^4)$$

$$F_{view} = \text{대향 평판 형상인수}(W, H, G)$$

### 5.3 프레임 전도

$$Q_{cond} = k_{frame} \cdot \frac{A_{frame}}{G/1000} \cdot (T_{cond} - T_{evap})$$

### 5.4 내부 혼합

$$Q_{mix} = \eta_{mix}(G, \text{mode}) \cdot Q_{useful}$$

### 5.5 순 냉방능력

$$Q_{net} = Q_{useful} - Q_{rad} - Q_{cond} - Q_{mix}$$

$$\text{Cap} [\%] = \frac{Q_{net}}{Q_{ref}} \times 100$$

---

## 6. 비말동반 모델 (Module B)

### 6.1 Weber 수 기반 이탈

$$We_{ch} = \frac{\rho_{air} \cdot V_{ch}^2 \cdot F_p}{\sigma_w}$$

$$We_{crit} = \frac{0.023}{\sqrt{F_p}}$$

$$V_{onset}: We_{ch}(V_{onset}) = We_{crit}$$

### 6.2 이탈 효율

$$\eta_{co} = C_{eta} \cdot (V_{face} - V_{onset})^2 \quad (V > V_{onset})$$

### 6.3 Monte Carlo 궤적

N개 액적 (bimodal 크기 분포) 발사:

$$m\ddot{x} = F_{drag,x}, \quad m\ddot{y} = F_{drag,y} - mg$$

$$F_{drag} = \frac{1}{2} C_D \rho_{air} A_p |V_{rel}| V_{rel}$$

$$P_{reach} = \frac{N_{도달}}{N_{total}}$$

### 6.4 Carryover Flux 및 Risk

$$\text{Flux} = q_{cond} \cdot \eta_{co} \cdot P_{reach} / A_{face} \times 1000$$

| Flux | 등급 |
|------|------|
| < 1.0 | SAFE |
| 1.0~4.0 | CAUTION |
| 4.0~12.0 | DANGER |
| ≥ 12.0 | SEVERE |

### 6.5 응축기 Q_carry 페널티

$$Q_{carry} = q_{cond} \cdot \eta_{co} \cdot P_{reach} \cdot h_{fg}$$

---

## 7. 습공기 물성

### 7.1 CoolProp 기반 (1차)

`HAPropsSI` — IAPWS-IF97 기반 정밀 계산

$$W = \text{HAPropsSI}('W', 'T', T, 'R', RH, 'P', P_{atm})$$
$$h = \text{HAPropsSI}('H', 'T', T, 'R', RH, 'P', P_{atm})$$
$$T_{dp} = \text{HAPropsSI}('D', 'T', T, 'R', RH, 'P', P_{atm})$$

### 7.2 Fallback 근사 (CoolProp 미설치 시)

포화 수증기압 (ASHRAE):

$$\ln P_{sat} = \frac{-5800.2206}{T_K} + 1.3915 - 0.04864 T_K + 4.176 \times 10^{-5} T_K^2 - 1.445 \times 10^{-8} T_K^3 + 6.546 \ln T_K$$

절대습도: $W = 0.62198 \cdot P_w / (P_{atm} - P_w)$

엔탈피: $h = 1006T + W(2{,}501{,}000 + 1860T)$ [J/kg]

---

## 8. 참고문헌

| # | 저자 | 제목 | 출처 |
|---|------|------|------|
| 1 | **Chen, J.C.** | **Correlation for boiling heat transfer to saturated fluids in convective flow** | **J. Heat Transfer 88(2), 1966** |
| 2 | Threlkeld, J.L. | Thermal Environmental Engineering | Prentice-Hall, 1970 |
| 3 | Wang, C.C. et al. | Heat transfer and friction correlations for plain fin-and-tube heat exchangers | IJHMT 43(15), 2000 |
| 4 | Wang, C.C. et al. | Airside performance of herringbone wavy fin-and-tube heat exchangers | IJHMT 42, 1999 |
| 5 | Wang, C.C. et al. | Heat transfer and friction correlation for compact louvered fin-and-tube heat exchangers | IJHMT 42(1), 1999 |
| 6 | Wang, C.C. et al. | An experimental study of convective heat transfer from a slit fin surface | IJHMT 44, 2001 |
| 7 | Shah, M.M. | A general correlation for heat transfer during film condensation inside pipes | IJHMT 22, 1979 |
| 8 | Shah, M.M. | Chart correlation for saturated boiling heat transfer | ASHRAE Trans 88, 1982 |
| 9 | Gungor, K.E. & Winterton, R.H.S. | A general correlation for flow boiling in tubes and annuli | IJHMT 29(3), 1986 |
| 10 | Gnielinski, V. | New equations for heat and mass transfer in turbulent pipe and channel flow | Int. Chem. Eng. 16, 1976 |
| 11 | **Kim, S.M. & Mudawar, I.** | **Universal approach to predicting saturated flow boiling heat transfer in mini/micro-channels** | **IJHMT 58, 2013** |
| 12 | **Kim, S.M. & Mudawar, I.** | **Universal approach to predicting heat transfer in condensing mini/micro-channels** | **IJHMT 55, 2012** |
| 13 | Cooper, M.G. | Heat flow rates in saturated nucleate pool boiling | Advances in Heat Transfer 16, 1984 |
| 14 | Kim, N.H., Youn, B. & Webb, R.L. | Air-side heat transfer and friction correlations for plain fin-and-tube HX | ASME JHT 121(3), 1999 |
| 15 | Chang, Y.J. & Wang, C.C. | A generalized heat transfer correlation for louver fin geometry | IJHMT 40(3), 1997 |
| 16 | Kays, W.M. & London, A.L. | Compact Heat Exchangers | McGraw-Hill, 1984 |
| 17 | Schmidt, T.E. | Heat transfer calculations for extended surfaces | Refrigerating Eng. 57, 1949 |
| 18 | ASHRAE | ASHRAE Handbook — Fundamentals, Ch.23 | ASHRAE, 2021 |
