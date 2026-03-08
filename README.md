# ❄️ 증발기-응축기 간격 통합 시뮬레이터

Gap between Evaporator and Condenser - Integrated Simulator

## Features
- **Module A** — Gap 열손실 → Cap Retention
- **Module B** — 비말동반 → Risk 등급
- **Module C** — 압력강하 → 시스템 ΔP
- **FT (Fin-Tube) / MCHX (Micro-Channel)** 모두 지원
- Threlkeld 정석 습면/건면 ε-NTU + Tube-Row Segmented 모델

## Run locally
```bash
pip install -r requirements.txt
streamlit run gap_simulator_app.py
```

## Deploy
Deployed on [Streamlit Cloud](https://streamlit.io/cloud)
