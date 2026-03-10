# ❄️ 증발기-응축기 간격 통합 시뮬레이터

Gap between Evaporator and Condenser - Integrated Simulator

## Features
- **Module A** — Gap 열손실 → Cap Retention (Level 2 Tube-Segment)
- **Module B** — 비말동반 → Risk 등급
- **Module C** — 압력강하 → 시스템 ΔP
- **상관식 자동 선택** — 기하·운전조건 기반

## Correlation Engine

| 카테고리 | D ≥ 3mm (FT) | D < 3mm (MCHX) |
|---------|--------------|-----------------|
| 공기 j | Wang(2000) plain / Wang(1999) wavy,louver | Chang & Wang (1997) |
| 공기 f | Wang(2000) / KYW(1999) | Chang & Wang (1997) |
| 냉매 증발 | **Chen(1966)** F×h_l×C_nb | **Kim & Mudawar(2013)** |
| 냉매 응축 | Shah(1979) | **Kim & Mudawar(2012)** |
| 단상 | Gnielinski(1976) | Gnielinski(1976) |
| 핀효율 | Schmidt 등가원형핀 (습면: T_fin_avg 반복) | 직선핀 |

## Level 2 Coil Model
- Nr×Nt×N_seg 세그먼트 개별 계산
- T_wall 반복 수렴 (α=0.7, Q_air=Q_ref)
- 엔탈피 기반 Row간 공기 상태 업데이트
- 습면 핀효율: b@T_fin_avg 반복 수렴 (물리적 정확)
- Counter/Parallel flow 증발기·응축기 독립 설정

## Run locally
```bash
pip install -r requirements.txt
streamlit run gap_simulator_app.py
```
