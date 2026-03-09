# ❄️ 증발기-응축기 간격 통합 시뮬레이터

Gap between Evaporator and Condenser - Integrated Simulator

## Features
- **Module A** — Gap 열손실 → Cap Retention
- **Module B** — 비말동반 → Risk 등급
- **Module C** — 압력강하 → 시스템 ΔP
- **🔬 Level 2 코일 모델** — Tube-Segment 상세 해석 (160 element)
- **상관식 자동 선택** — 기하·운전조건 기반, 사용자 입력 불필요
- **FT (Fin-Tube) / MCHX (Micro-Channel)** 모두 지원
- Threlkeld 정석 습면/건면 ε-NTU + T_wall 반복 수렴

## Correlation Auto-Selection
`select_correlations(spec, geo, ref, side)` — 입력 스펙만으로 최적 상관식 자동 매핑

| 카테고리 | 상관식 |
|---------|--------|
| 공기 j (FT) | Wang(2000) plain / Wang(1999) wavy,louver / E-type slit |
| 공기 j (MCHX) | Chang & Wang (1997) |
| 공기 f (FT) | Wang(2000) Pt/Pl≤1.35 / KYW(1999) Pt/Pl>1.35 |
| 냉매 증발 | Gungor-Winterton(1986) / Shah(1982) 선택 가능 |
| 냉매 응축 | Shah(1979) |
| 단상 (과열/과냉) | Gnielinski(1976) |
| 핀효율 | Schmidt 등가원형핀 (FT) / 직선핀 (MCHX) |

## Level 2 Coil Model
`compute_coil_v3(spec, geo, ref, T_air, RH, V, m_ref, x_in, ...)`

- Nr×Nt×N_seg = 160 세그먼트 개별 계산
- 냉매 건도(x) 세그먼트별 추적: 과냉→이상→과열 자동 분기
- T_wall 반복 수렴 (Q_air = Q_ref, successive substitution)
- Counter/Parallel flow 선택
- 습면 Threlkeld 정석 (b-factor는 η_fin에만 반영, h_o 불변)

## Run locally
```bash
pip install -r requirements.txt
streamlit run gap_simulator_app.py
```

## Deploy
Deployed on [Streamlit Cloud](https://streamlit.io/cloud)
