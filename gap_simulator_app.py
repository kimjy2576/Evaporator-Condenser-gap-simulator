"""
gap_simulator_app.py — 증발기-응축기 간격 통합 시뮬레이터 (Streamlit)

실행: streamlit run gap_simulator_app.py
"""
import streamlit as st
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings, os, sys, io, json
warnings.filterwarnings("ignore")

for _p in ["/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc",
           "C:/Windows/Fonts/malgun.ttf"]:
    if os.path.exists(_p):
        from matplotlib import font_manager as fm
        fm.fontManager.addfont(_p)
        matplotlib.rcParams['font.family'] = fm.FontProperties(fname=_p).get_name()
        break
matplotlib.rcParams['axes.unicode_minus'] = False
matplotlib.rcParams['mathtext.fontset'] = 'cm'
matplotlib.rcParams['figure.dpi'] = 200
matplotlib.rcParams['savefig.dpi'] = 200

from common import (
    FinTubeSpec, MCHXSpec, RefrigerantState, DARK, C,
    compute_ft_geometry, compute_mchx_geometry, compute_UA,
    compute_coil_performance_segmented,
    FLUX_CAUTION, FLUX_DANGER, FLUX_SEVERE,
)
from module_a import GapParams, simulate_gap, sweep
from module_b import (analyze_combined, compute_carry_penalty, apply_carry_penalty,
                      CarryoverSpec, MCHXCarryoverSpec,
                      we_ch, we_crit, v_onset, eta_co, monte_carlo)
from module_c import InletCondition, sweep_dp
from run_single import build_spec
from visualize import (make_single_fig_a, make_single_fig_b, make_single_fig_c,
                       make_single_schematic)

st.set_page_config(page_title="Gap Simulator", page_icon="❄️", layout="wide")

# ═══════════════════════════════════════════════════════════
#  비밀번호 보호
# ═══════════════════════════════════════════════════════════
import hmac

def check_password():
    def password_entered():
        try:
            correct = st.secrets["password"]
        except Exception:
            correct = "gap1234"  # 로컬 실행 시 기본 비밀번호
        if hmac.compare_digest(st.session_state["password"], correct):
            st.session_state["password_correct"] = True
            del st.session_state["password"]
        else:
            st.session_state["password_correct"] = False

    if st.session_state.get("password_correct", False):
        return True

    st.title("❄️ 증발기-응축기 간격 통합 시뮬레이터")
    st.text_input("🔒 비밀번호를 입력하세요", type="password",
                   on_change=password_entered, key="password")
    if "password_correct" in st.session_state:
        st.error("비밀번호가 틀렸습니다")
    return False

if not check_password():
    st.stop()

st.title("❄️ 증발기-응축기 간격 통합 시뮬레이터")

# ═══════════════════════════════════════════════════════════
#  JSON → session_state
# ═══════════════════════════════════════════════════════════
def apply_json_to_state(cfg):
    e = cfg.get('evap',{}); c = cfg.get('cond',{}); r = cfg.get('ref',{}); g = cfg.get('gap',{}); s = cfg.get('sim',{})
    if e.get('hx_type'): st.session_state['evap_hx_type'] = e['hx_type']
    if c.get('hx_type'): st.session_state['cond_hx_type'] = c['hx_type']
    _fin_list = ["plain","wavy","slit","louvered"]
    _lay_list = ["staggered","inline"]
    if e.get('fin_type') and e['fin_type'] in _fin_list: st.session_state['efin'] = _fin_list.index(e['fin_type'])
    if e.get('tube_layout') and e['tube_layout'] in _lay_list: st.session_state['elay'] = _lay_list.index(e['tube_layout'])
    if e.get('fpi'): st.session_state['efpi'] = int(e['fpi'])
    if e.get('tube_rows'): st.session_state['eNr'] = int(e['tube_rows'])
    if e.get('tube_cols'): st.session_state['eNt'] = int(e['tube_cols'])
    for k,sk in [('W','eW'),('H','eH'),('D','eD')]:
        if k in e: v=e[k]; st.session_state[sk]=int(v*1000) if v<1 else int(v)
    if e.get('ch_width'): st.session_state['e_ch_w']=e['ch_width']
    if e.get('ch_height'): st.session_state['e_ch_h']=e['ch_height']
    if e.get('louver_pitch'): st.session_state['e_lp']=e['louver_pitch']
    if e.get('louver_angle'): st.session_state['e_la']=e['louver_angle']
    if e.get('n_ports'): st.session_state['e_nports']=int(e['n_ports'])
    if e.get('n_slabs'): st.session_state['e_nslabs']=int(e['n_slabs'])
    if c.get('fin_type') and c['fin_type'] in _fin_list: st.session_state['cfin'] = _fin_list.index(c['fin_type'])
    if c.get('tube_layout') and c['tube_layout'] in _lay_list: st.session_state['clay'] = _lay_list.index(c['tube_layout'])
    if c.get('fpi'): st.session_state['cfpi'] = int(c['fpi'])
    if c.get('tube_rows'): st.session_state['cNr'] = int(c['tube_rows'])
    if c.get('tube_cols'): st.session_state['cNt'] = int(c['tube_cols'])
    for k,sk in [('W','cW'),('H','cH'),('D','cD')]:
        if k in c: v=c[k]; st.session_state[sk]=int(v*1000) if v<1 else int(v)
    if g.get('T_amb'): st.session_state['T_amb']=g['T_amb']
    if g.get('RH_in'): st.session_state['RH_in']=g['RH_in']
    if g.get('CMM'): st.session_state['CMM']=g['CMM']
    if g.get('gap_mode'): st.session_state['gap_mode']=["open","semi","sealed"].index(g['gap_mode'])
    if g.get('seal_fraction'): st.session_state['seal_frac']=g['seal_fraction']
    _rlist = ["R410A","R32","R134a","R290","R1234yf"]
    if r.get('refrigerant') and r['refrigerant'] in _rlist: st.session_state['refr']=_rlist.index(r['refrigerant'])
    if r.get('T_sat_evap'): st.session_state['T_evap']=r['T_sat_evap']
    if r.get('T_sat_cond'): st.session_state['T_cond']=r['T_sat_cond']
    if s.get('gap_min'): st.session_state['gap_min']=s['gap_min']
    if s.get('gap_max'): st.session_state['gap_max']=s['gap_max']
    if s.get('gap_points'): st.session_state['gap_pts']=s['gap_points']

# ═══════════════════════════════════════════════════════════
#  Sidebar
# ═══════════════════════════════════════════════════════════
with st.sidebar:
    st.header("🔧 설정")

    with st.expander("📁 JSON 파일", expanded=False):
        uploaded = st.file_uploader("JSON 불러오기", type=['json'], key='json_upload')
        if uploaded is not None:
            try:
                loaded = json.loads(uploaded.read().decode('utf-8'))
                apply_json_to_state(loaded)
                st.success(f"✅ {uploaded.name}")
                st.rerun()
            except Exception as ex:
                st.error(f"❌ {ex}")
        # JSON 저장 — build_config는 아래에서 정의되므로 placeholder
        save_json_placeholder = st.empty()

    with st.expander("증발기", expanded=True):
        evap_hx_type = st.radio("HX 타입", ["FT","MCHX"], horizontal=True, key='evap_hx_type')
        # 공통: 외형 치수
        c1,c2,c3 = st.columns(3)
        with c1: evap_W = st.number_input("W[mm]",0,1000,300,key="eW")
        with c2: evap_H = st.number_input("H[mm]",0,800,250,key="eH")
        with c3: evap_D = st.number_input("D[mm]",0,200,45,key="eD")

        if evap_hx_type == "FT":
            c1,c2 = st.columns(2)
            with c1:
                evap_fin = st.selectbox("핀 타입",["plain","wavy","slit","louvered"],index=1,key='efin')
                evap_layout = st.selectbox("배열",["staggered","inline"],key='elay')
            with c2:
                evap_Nr = st.number_input("Row",1,6,2,key="eNr")
                evap_Nt = st.number_input("Col",2,30,8,key="eNt")
            # 핀/튜브 기하
            c1,c2 = st.columns(2)
            with c1:
                evap_fpi = st.slider("FPI",0,30,14,key='efpi')
                e_fin_t = st.number_input("핀 두께 [mm]",0.0,0.30,0.10,0.01,key='e_fin_t')
                e_Do = st.number_input("튜브 Do [mm]",0.0,15.0,9.52,0.01,key='e_Do')
            with c2:
                e_Di = st.number_input("튜브 Di [mm]",0.0,14.0,8.52,0.01,key='e_Di')
                e_drain = st.number_input("드레인 높이 [mm]",0.0,50.0,25.0,1.0,key='e_drain')
            c1,c2 = st.columns(2)
            with c1: e_Pt = st.number_input("Pitch_t [mm]",0.0,40.0,25.4,0.1,key='e_Pt')
            with c2: e_Pl = st.number_input("Pitch_l [mm]",0.0,40.0,22.0,0.1,key='e_Pl')
            # 핀 타입별 상세
            if evap_fin in ['wavy']:
                c1,c2 = st.columns(2)
                with c1: e_wavy_h = st.number_input("파형 높이 [mm]",0.0,5.0,1.5,0.1,key='e_wavy_h')
                with c2: e_wavy_a = st.number_input("파형 각도 [°]",0.0,30.0,17.5,0.5,key='e_wavy_a')
            if evap_fin in ['louvered','louver']:
                c1,c2 = st.columns(2)
                with c1: e_lv_p = st.number_input("루버 피치 [mm]",0.0,4.0,1.7,0.1,key='e_lv_p')
                with c2: e_lv_a = st.number_input("루버 각도 [°]",0.0,40.0,28.0,1.0,key='e_lv_a')
            if evap_fin == 'slit':
                c1,c2 = st.columns(2)
                with c1: e_slit_n = st.number_input("슬릿 수",2,12,6,key='e_slit_n')
                with c2: e_slit_h = st.number_input("슬릿 높이 [mm]",0.0,3.0,1.0,0.1,key='e_slit_h')
            # 재질
            c1,c2 = st.columns(2)
            with c1: e_fin_mat = st.selectbox("핀 재질",["Al","Cu"],key='e_fin_mat')
            with c2: e_tube_mat = st.selectbox("튜브 재질",["Cu","Al"],key='e_tube_mat')
        else:
            # MCHX
            c1,c2 = st.columns(2)
            with c1:
                e_ch_w = st.number_input("Ch Width [mm]",0.0,3.0,0.8,0.1,key='e_ch_w')
                e_ch_h = st.number_input("Ch Height [mm]",0.0,5.0,1.5,0.1,key='e_ch_h')
                e_ch_wall = st.number_input("Ch Wall [mm]",0.0,1.0,0.4,0.1,key='e_ch_wall')
                e_lp = st.number_input("Louver Pitch [mm]",0.0,3.0,1.2,0.1,key='e_lp')
            with c2:
                e_la = st.number_input("Louver Angle [°]",0.0,40.0,27.0,1.0,key='e_la')
                e_nports = st.number_input("N_ports",4,40,12,key='e_nports')
                e_nslabs = st.number_input("N_slabs",1,6,1,key='e_nslabs')
                e_slab_p = st.number_input("Slab Pitch [mm]",0.0,15.0,8.0,0.5,key='e_slab_p')
            c1,c2 = st.columns(2)
            with c1:
                evap_fpi = st.slider("FPI",0,35,20,key='efpi_mchx')
                e_mchx_ft = st.number_input("핀 두께 [mm] ",0.0,0.20,0.06,0.01,key='e_mchx_ft')
            with c2:
                e_fin_mat = st.selectbox("핀 재질 ",["Al","Cu"],key='e_mchx_fmat')
                e_tube_mat = st.selectbox("튜브 재질 ",["Al","Cu"],key='e_mchx_tmat')

    with st.expander("응축기", expanded=False):
        cond_hx_type = st.radio("HX 타입 ",["FT","MCHX"],horizontal=True,key='cond_hx_type')
        c1,c2,c3 = st.columns(3)
        with c1: cond_W = st.number_input("W[mm] ",0,1000,350,key="cW")
        with c2: cond_H = st.number_input("H[mm] ",0,800,280,key="cH")
        with c3: cond_D = st.number_input("D[mm] ",0,200,50,key="cD")

        if cond_hx_type == "FT":
            c1,c2 = st.columns(2)
            with c1:
                cond_fin = st.selectbox("핀 타입 ",["plain","wavy","slit","louvered"],key="cfin")
                cond_layout = st.selectbox("배열 ",["staggered","inline"],key="clay")
            with c2:
                cond_Nr = st.number_input("Row ",1,6,2,key="cNr")
                cond_Nt = st.number_input("Col ",2,30,9,key="cNt")
            c1,c2 = st.columns(2)
            with c1:
                cond_fpi = st.slider("FPI ",0,30,14,key="cfpi")
                c_fin_t = st.number_input("핀 두께[mm]",0.0,0.30,0.10,0.01,key='c_fin_t')
                c_Do = st.number_input("튜브 Do[mm]",0.0,15.0,9.52,0.01,key='c_Do')
            with c2:
                c_Di = st.number_input("튜브 Di[mm]",0.0,14.0,8.52,0.01,key='c_Di')
                c_drain = st.number_input("드레인[mm]",0.0,50.0,25.0,1.0,key='c_drain')
            c1,c2 = st.columns(2)
            with c1: c_Pt = st.number_input("Pitch_t[mm]",0.0,40.0,25.4,0.1,key='c_Pt')
            with c2: c_Pl = st.number_input("Pitch_l[mm]",0.0,40.0,22.0,0.1,key='c_Pl')
            if cond_fin in ['wavy']:
                c1,c2 = st.columns(2)
                with c1: c_wavy_h = st.number_input("파형높이[mm]",0.0,5.0,1.5,0.1,key='c_wavy_h')
                with c2: c_wavy_a = st.number_input("파형각도[°]",0.0,30.0,17.5,0.5,key='c_wavy_a')
            if cond_fin in ['louvered','louver']:
                c1,c2 = st.columns(2)
                with c1: c_lv_p = st.number_input("루버피치[mm]",0.0,4.0,1.7,0.1,key='c_lv_p')
                with c2: c_lv_a = st.number_input("루버각도[°]",0.0,40.0,28.0,1.0,key='c_lv_a')
            if cond_fin == 'slit':
                c1,c2 = st.columns(2)
                with c1: c_slit_n = st.number_input("슬릿수 ",2,12,6,key='c_slit_n')
                with c2: c_slit_h = st.number_input("슬릿높이[mm]",0.0,3.0,1.0,0.1,key='c_slit_h')
            c1,c2 = st.columns(2)
            with c1: c_fin_mat = st.selectbox("핀재질",["Al","Cu"],key='c_fin_mat')
            with c2: c_tube_mat = st.selectbox("튜브재질",["Cu","Al"],key='c_tube_mat')
        else:
            c1,c2 = st.columns(2)
            with c1:
                c_ch_w = st.number_input("Ch Width[mm] ",0.0,3.0,0.8,0.1,key='c_ch_w')
                c_ch_h = st.number_input("Ch Height[mm] ",0.0,5.0,1.5,0.1,key='c_ch_h')
                c_ch_wall = st.number_input("Ch Wall[mm] ",0.0,1.0,0.4,0.1,key='c_ch_wall')
                c_lp = st.number_input("Louver Pitch[mm] ",0.0,3.0,1.2,0.1,key='c_lp')
            with c2:
                c_la = st.number_input("Louver Angle[°] ",0.0,40.0,27.0,1.0,key='c_la')
                c_nports = st.number_input("N_ports ",4,40,12,key='c_nports')
                c_nslabs = st.number_input("N_slabs ",1,6,1,key='c_nslabs')
                c_slab_p = st.number_input("Slab Pitch[mm] ",0.0,15.0,8.0,0.5,key='c_slab_p')
            c1,c2 = st.columns(2)
            with c1:
                cond_fpi = st.slider("FPI  ",0,35,20,key='cfpi_mchx')
                c_mchx_ft = st.number_input("핀두께[mm] ",0.0,0.20,0.06,0.01,key='c_mchx_ft')
            with c2:
                c_fin_mat = st.selectbox("핀재질 ",["Al","Cu"],key='c_mchx_fmat')
                c_tube_mat = st.selectbox("튜브재질 ",["Al","Cu"],key='c_mchx_tmat')

    with st.expander("운전 조건", expanded=True):
        T_amb = st.number_input("T_amb [°C]",0.0,45.0,27.0,0.5,key='T_amb')
        RH_in = st.slider("RH",0.0,0.95,0.70,0.05,key='RH_in')
        CMM = st.number_input("CMM",0.0,30.0,6.75,0.5,key='CMM')
        gap_mode = st.selectbox("Gap Mode",["open","semi","sealed"],index=1,key='gap_mode')
        seal_frac = st.slider("Seal Fraction",0.0,1.0,0.7,0.1,disabled=(gap_mode!="semi"),key='seal_frac')
        c1,c2 = st.columns(2)
        with c1: frame_mat = st.selectbox("프레임 재질",["Al","Steel"],key='frame_mat')
        with c2: A_frame = st.number_input("프레임 단면적 [cm²]",0.0,20.0,4.0,0.5,key='A_frame')

    with st.expander("냉매", expanded=False):
        refrigerant = st.selectbox("냉매",["R410A","R32","R134a","R290","R1234yf"],key='refr')
        T_evap = st.number_input("T_evap [°C]",-10.0,20.0,5.0,0.5,key='T_evap')
        T_cond = st.number_input("T_cond [°C]",0.0,65.0,45.0,0.5,key='T_cond')

    with st.expander("시뮬레이션 범위", expanded=False):
        gap_min = st.number_input("Gap Min [mm]",0,50,5,key='gap_min')
        gap_max = st.number_input("Gap Max [mm]",0,200,100,key='gap_max')
        gap_pts = st.slider("포인트 수",0,50,25,key='gap_pts')

    st.divider()
    run_btn = st.button("▶ 해석 실행", type="primary", use_container_width=True)
    st.divider()
    st.subheader("📊 비교 분석")
    if 'cases' not in st.session_state: st.session_state.cases = []
    add_case = st.button("+ 현재 설정 케이스 추가", use_container_width=True)
    if st.session_state.cases:
        st.caption(f"{len(st.session_state.cases)}개 케이스")
        for i,cc in enumerate(st.session_state.cases): st.text(f"  {i+1}. {cc['label']}")
        if st.button("🗑 전체 삭제"): st.session_state.cases=[]; st.rerun()

# ═══════════════════════════════════════════════════════════
#  Config
# ═══════════════════════════════════════════════════════════
def build_config():
    cfg = dict(
        ref=dict(refrigerant=refrigerant,T_sat_evap=T_evap,T_sat_cond=T_cond),
        gap=dict(T_amb=T_amb,RH_in=RH_in,CMM=CMM,gap_mode=gap_mode,seal_fraction=seal_frac,
                 frame_material=frame_mat,A_frame=A_frame),
        sim=dict(gap_min=gap_min,gap_max=gap_max,gap_points=gap_pts),
    )
    if evap_hx_type=="FT":
        ed = dict(hx_type='FT',W=evap_W/1000,H=evap_H/1000,D=evap_D/1000,
            fpi=evap_fpi,fin_type=evap_fin,tube_layout=evap_layout,
            tube_rows=evap_Nr,tube_cols=evap_Nt,
            fin_thickness=e_fin_t,
            tube_do=e_Do,tube_di=e_Di,
            tube_pitch_t=e_Pt,tube_pitch_l=e_Pl,
            fin_material=e_fin_mat,tube_material=e_tube_mat,
            h_drain=e_drain)
        if evap_fin=='wavy': ed.update(wavy_height=e_wavy_h,wavy_angle=e_wavy_a)
        if evap_fin in ['louvered','louver']: ed.update(ft_louver_pitch=e_lv_p,ft_louver_angle=e_lv_a)
        if evap_fin=='slit': ed.update(slit_num=e_slit_n,slit_height=e_slit_h)
        cfg['evap']=ed
    else:
        cfg['evap']=dict(hx_type='MCHX',W=evap_W/1000,H=evap_H/1000,D=evap_D/1000,
            ch_width=e_ch_w,ch_height=e_ch_h,ch_wall=e_ch_wall,
            louver_pitch=e_lp,louver_angle=e_la,
            n_ports=e_nports,n_slabs=e_nslabs,slab_pitch=e_slab_p,
            fin_pitch=25.4/evap_fpi,fin_thickness=e_mchx_ft,
            fin_material=e_fin_mat,tube_material=e_tube_mat)
    if cond_hx_type=="FT":
        cd = dict(hx_type='FT',W=cond_W/1000,H=cond_H/1000,D=cond_D/1000,
            fpi=cond_fpi,fin_type=cond_fin,tube_layout=cond_layout,
            tube_rows=cond_Nr,tube_cols=cond_Nt,
            fin_thickness=c_fin_t,
            tube_do=c_Do,tube_di=c_Di,
            tube_pitch_t=c_Pt,tube_pitch_l=c_Pl,
            fin_material=c_fin_mat,tube_material=c_tube_mat,
            h_drain=c_drain)
        if cond_fin=='wavy': cd.update(wavy_height=c_wavy_h,wavy_angle=c_wavy_a)
        if cond_fin in ['louvered','louver']: cd.update(ft_louver_pitch=c_lv_p,ft_louver_angle=c_lv_a)
        if cond_fin=='slit': cd.update(slit_num=c_slit_n,slit_height=c_slit_h)
        cfg['cond']=cd
    else:
        cfg['cond']=dict(hx_type='MCHX',W=cond_W/1000,H=cond_H/1000,D=cond_D/1000,
            ch_width=c_ch_w,ch_height=c_ch_h,ch_wall=c_ch_wall,
            louver_pitch=c_lp,louver_angle=c_la,
            n_ports=c_nports,n_slabs=c_nslabs,slab_pitch=c_slab_p,
            fin_pitch=25.4/cond_fpi,fin_thickness=c_mchx_ft,
            fin_material=c_fin_mat,tube_material=c_tube_mat)
    return cfg

if save_json_placeholder:
    cfg_for_save = build_config()
    json_str = json.dumps(cfg_for_save, indent=2, ensure_ascii=False)
    save_json_placeholder.download_button(
        "💾 현재 설정 JSON 저장", json_str,
        "gap_config.json", "application/json",
        use_container_width=True)

# ═══════════════════════════════════════════════════════════
#  해석 엔진
# ═══════════════════════════════════════════════════════════
def build_hx(d, defs):
    if d.get('hx_type')=='MCHX':
        return MCHXSpec(W=d['W'],H=d['H'],D=d['D'],
            ch_width=d.get('ch_width',0.8)*1e-3, ch_height=d.get('ch_height',1.5)*1e-3,
            ch_wall=d.get('ch_wall',0.4)*1e-3, n_ports=int(d.get('n_ports',12)),
            fin_pitch=d.get('fin_pitch',1.0)*1e-3, fin_thickness=d.get('fin_thickness',0.06)*1e-3,
            louver_pitch=d.get('louver_pitch',1.2)*1e-3, louver_angle=d.get('louver_angle',27.0),
            n_slabs=int(d.get('n_slabs',1)), slab_pitch=d.get('slab_pitch',8.0)*1e-3)
    return build_spec({**defs, **d})

def run_analysis(cfg):
    DEF=dict(fin_pitch=1.8,fin_thickness=0.10,fin_height=25.4,tube_do=9.52,tube_di=8.52,
        tube_pitch_t=25.4,tube_pitch_l=22.0,corr_mode="etype",wavy_height=1.5,wavy_angle=17.5,
        ft_louver_pitch=1.7,ft_louver_angle=28.0,slit_num=6,slit_height=1.0,
        fin_material="Al",tube_material="Cu",h_drain=25.0)
    rd=cfg.get('ref',{}); gd=cfg.get('gap',{}); sd=cfg.get('sim',{})
    evap=build_hx(cfg.get('evap',{}),DEF); cond=build_hx(cfg.get('cond',{}),DEF)
    ref=RefrigerantState(**rd)
    gp=GapParams(T_amb=gd['T_amb'],RH_in=gd['RH_in'],CMM=gd['CMM'],
        gap_mode=gd['gap_mode'],seal_fraction=gd.get('seal_fraction',0.7),
        mode='forced',frame_material=gd.get('frame_material','Al'),
        A_frame=gd.get('A_frame',4.0)*1e-4)
    gaps=np.linspace(sd['gap_min'],sd['gap_max'],int(sd['gap_points']))
    geo_e=compute_ft_geometry(evap) if evap.hx_type=='FT' else compute_mchx_geometry(evap)
    geo_c=compute_ft_geometry(cond) if cond.hx_type=='FT' else compute_mchx_geometry(cond)
    V_face=gp.CMM/(60.0*evap.W*evap.H)
    ua_e=compute_UA(evap,geo_e,ref,V_face,'evap'); ua_c=compute_UA(cond,geo_c,ref,V_face,'cond')
    et=f"{evap.fin_type.capitalize()} {evap.tube_layout[0].upper()}" if evap.hx_type=='FT' else 'MCHX'
    tag=f"{et} | {ref.refrigerant} | T_evap={ref.T_sat_evap}°C | CMM={gp.CMM} | {gp.gap_mode}"
    results_a=sweep(gaps,evap,cond,geo_e,geo_c,ua_e,ua_c,ref,gp)
    case_b=dict(name=tag,T_in=gd['T_amb'],RH_in=gd['RH_in'],CMM=gd['CMM'],T_wall=ref.T_sat_evap,
        gap_mm=20,theta_deg=0,eta_coeff=5e-4,w_bridge=0.75,gap_mode=gd['gap_mode'],
        seal_fraction=gd.get('seal_fraction',0.7),hx_type=evap.hx_type,
        tube_layout=getattr(evap,'tube_layout','staggered'))
    result_b=analyze_combined(case_b,gaps,evap,geo_e,cond,geo_c,ua_e,ua_c,ref)
    cp=compute_carry_penalty(result_b,evap,geo_e,gaps,ref.T_sat_cond); apply_carry_penalty(results_a,cp)
    inlet=InletCondition(T_in=gd['T_amb'],RH_in=gd['RH_in'],CMM=gd['CMM'],T_wall_evap=ref.T_sat_evap)
    result_c=sweep_dp(gaps,evap,cond,geo_e,geo_c,inlet)
    co=CarryoverSpec(evap,geo_e) if evap.hx_type=='FT' else MCHXCarryoverSpec(evap,geo_e)
    Wc=we_crit(co)+1e-9; Von=v_onset(co) if evap.hx_type=='FT' else 2.5
    Pr=np.array([monte_carlo(g,V_face,co,0,0.75,N=40,seed=42)['P_reach'] for g in gaps])
    Vs=np.linspace(0.3,max(V_face*2.5,6),60); Ws=np.array([we_ch(v,co)/Wc for v in Vs])
    vb=dict(P_reach_arr=Pr,We_ratio=we_ch(V_face,co)/Wc,V_face=V_face,V_onset=Von,
        q_cond=result_b['cr']['q_cond'],carry_penalty=cp,V_sweep=Vs,We_sweep=Ws,gaps=gaps)
    return dict(gaps=gaps,results_a=results_a,result_b=result_b,result_c=result_c,
        carry_penalty=cp,viz_b=vb,evap=evap,cond=cond,ref=ref,gp=gp,
        geo_e=geo_e,geo_c=geo_c,ua_e=ua_e,ua_c=ua_c,tag=tag,V_face=V_face,
        gap_d=gd,inlet=inlet,label=tag)

def make_df(res):
    rows=[]
    for i,g in enumerate(res['gaps']):
        r=res['results_a'][i]; fl=res['result_b']['flux_arr'][i]; rc=res['result_c']
        rk='SAFE' if fl<FLUX_CAUTION else 'CAUTION' if fl<FLUX_DANGER else 'DANGER' if fl<FLUX_SEVERE else 'SEVERE'
        rows.append(dict(Gap=f"{g:.1f}",Cap=f"{r['cap_ratio']:.1f}",Q_net=f"{r['Q_net']:.0f}",
            Q_sen=f"{r['Q_sensible']:.0f}",Q_lat=f"{r['Q_latent']:.0f}",SHR=f"{r['SHR']:.3f}",
            Dehumid=f"{r['dehumid_rate']:.0f}",dT=f"{r['dT_recir']:.2f}",
            Flux=f"{fl:.2f}",Risk=rk,dP=f"{rc['dp_total'][i]:.1f}"))
    return pd.DataFrame(rows)

# ═══════════════════════════════════════════════════════════
#  실행
# ═══════════════════════════════════════════════════════════
cfg=build_config()
if add_case:
    with st.spinner("케이스 저장..."): st.session_state.cases.append(run_analysis(cfg))
    st.rerun()
if run_btn:
    with st.spinner("해석 중..."): st.session_state.results=run_analysis(cfg)

# ═══════════════════════════════════════════════════════════
#  결과
# ═══════════════════════════════════════════════════════════
if 'results' in st.session_state and st.session_state.results:
  try:
    res=st.session_state.results; gaps=res['gaps']; tag=res['tag']
    i50=np.argmin(np.abs(gaps-50)); r50=res['results_a'][i50]
    c1,c2,c3,c4,c5=st.columns(5)
    c1.metric("Cap(G=50)",f"{r50['cap_ratio']:.1f}%")
    c2.metric("Q_net",f"{r50['Q_net']:.0f}W")
    c3.metric("SHR",f"{r50['SHR']:.3f}")
    c4.metric("ΔP",f"{res['result_c']['dp_total'][i50]:.0f}Pa")
    fl50=res['result_b']['flux_arr'][i50]
    c5.metric("Risk",('SAFE' if fl50<FLUX_CAUTION else 'CAUTION' if fl50<FLUX_DANGER else 'DANGER' if fl50<FLUX_SEVERE else 'SEVERE'))
    st.caption(f"**{tag}** | V={res['V_face']:.2f}m/s | UA_e={res['ua_e']['UA']:.1f}W/K | η_o={res['ua_e']['eta_o']:.3f}")

    # 고화질 렌더링 함수
    def show_fig(fig, dpi=200):
        buf = io.BytesIO()
        fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight', facecolor=fig.get_facecolor())
        buf.seek(0)
        st.image(buf, use_container_width=True)
        plt.close(fig)

    ta,tb,tc,ts,tcmp,td=st.tabs(["🌡️ A 냉방","💧 B 비말","📊 C 압강","🖼️ 개략도","📈 비교","📋 데이터"])
    with ta:
        f=make_single_fig_a(gaps,res['results_a'],res['gp'],res['ref'],tag); show_fig(f)
    with tb:
        f=make_single_fig_b(gaps,res['result_b'],res['viz_b'],tag); show_fig(f)
    with tc:
        f=make_single_fig_c(gaps,res['result_c'],res['inlet'],tag); show_fig(f)
    with ts:
        gi=20 if 20<=gaps[-1] else gaps[len(gaps)//2]; ii=np.argmin(np.abs(gaps-gi))
        f=make_single_schematic(res['evap'],res['cond'],res['geo_e'],res['geo_c'],
            res['ua_e'],res['ua_c'],res['gp'],res['ref'],res['results_a'][ii],res['viz_b'],tag)
        if f:
            show_fig(f)
        elif os.path.exists("gap_schematic.png"):
            st.image("gap_schematic.png",use_container_width=True)
    with tcmp:
        cs=st.session_state.get('cases',[])
        if len(cs)<2: st.info("비교: 사이드바에서 **+ 케이스 추가** 2회 이상")
        else:
            from run_compare import compare_fig_a,compare_fig_b,compare_fig_c
            st.subheader(f"📈 {len(cs)}개 케이스 비교")
            sa,sb,sc=st.tabs(["A","B","C"])
            with sa: compare_fig_a(cs); st.image("compare_a.png",use_container_width=True) if os.path.exists("compare_a.png") else None
            with sb: compare_fig_b(cs); st.image("compare_b.png",use_container_width=True) if os.path.exists("compare_b.png") else None
            with sc: compare_fig_c(cs); st.image("compare_c.png",use_container_width=True) if os.path.exists("compare_c.png") else None
    with td:
        df=make_df(res); st.dataframe(df,use_container_width=True,hide_index=True)
        st.download_button("📥 CSV",df.to_csv(index=False),"gap_results.csv","text/csv")
        st.subheader("설계 조건")
        st.json({"냉매":res['ref'].refrigerant,"T_evap/T_cond":f"{res['ref'].T_sat_evap}/{res['ref'].T_sat_cond}°C",
            "T_amb":f"{res['gap_d']['T_amb']}°C","RH":f"{res['gap_d']['RH_in']:.0%}","CMM":str(res['gap_d']['CMM']),
            "Mode":res['gap_d']['gap_mode'],"V_face":f"{res['V_face']:.2f}m/s",
            "Evap":res['evap'].hx_type,"Cond":res['cond'].hx_type,
            "UA_e":f"{res['ua_e']['UA']:.1f}W/K","UA_c":f"{res['ua_c']['UA']:.1f}W/K"})
  except Exception as e:
    st.error(f"결과 표시 오류: {e}")
    st.info("설정 변경 후 **▶ 해석 실행**을 다시 클릭하세요.")
    if st.button("🔄 결과 초기화"):
        del st.session_state['results']
        st.rerun()
else:
    st.info("👈 사이드바에서 설정 후 **▶ 해석 실행** 클릭")
    st.markdown("""
    ### 시뮬레이터 소개
    - **Module A** — Gap 열손실 → Cap Retention
    - **Module B** — 비말동반 → Risk 등급
    - **Module C** — 압력강하 → 시스템 ΔP
    
    Threlkeld 정석 ε-NTU + Tube-Row Segment. **FT / MCHX** 모두 지원.
    """)

# ═══════════════════════════════════════════════════════════
#  문서 (항상 표시)
# ═══════════════════════════════════════════════════════════
st.divider()
st.header("📚 문서")
doc1,doc2,doc3,doc4=st.tabs(["📘 사용 가이드","📋 평가서","🏗️ 구조 설명서","📐 수학적 모델"])

def _load_doc(name):
    for p in [f"docs/{name}", name]:
        if os.path.exists(p):
            with open(p,'r',encoding='utf-8') as f: return f.read()
    return f"⚠️ `{name}` 파일을 찾을 수 없습니다. `docs/` 폴더를 확인하세요."

def _safe_md(name, tab):
    with tab:
        try:
            txt = _load_doc(name)
            # Streamlit markdown에서 $ 충돌 방지
            st.markdown(txt, unsafe_allow_html=True)
        except Exception as e:
            st.error(f"문서 렌더링 오류: {e}")

_safe_md("DOC_GUIDE.md", doc1)
_safe_md("DOC_EVAL.md", doc2)
_safe_md("DOC_STRUCTURE.md", doc3)
_safe_md("DOC_MATH.md", doc4)