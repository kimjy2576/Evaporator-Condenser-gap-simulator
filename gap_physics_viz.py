"""
Gap 물리 현상 가시화 — 증발기-응축기 간격 열유동 개략도
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Arc, Circle
import matplotlib.patheffects as pe
import os

_NOTO = "/usr/share/fonts/opentype/noto/NotoSansCJK-Regular.ttc"
if os.path.exists(_NOTO):
    from matplotlib import font_manager as fm
    fm.fontManager.addfont(_NOTO)
    import matplotlib; matplotlib.rcParams['font.family'] = fm.FontProperties(fname=_NOTO).get_name()
import matplotlib; matplotlib.rcParams['axes.unicode_minus'] = False

BG     = "#0a0e17"
PANEL  = "#0d1420"
CYAN   = "#00c8ff"
RED    = "#ff5e3a"
YEL    = "#ffc940"
GRN    = "#2adf8a"
PUR    = "#d070ff"
ORA    = "#ff9f43"
DIM    = "#4a6070"
TEXT   = "#c8d8e8"
BORD   = "#1e2d42"
WHITE  = "#e8f0f8"

def _glow(color, lw=2):
    return [pe.withStroke(linewidth=lw+2, foreground=color+'44')]


def draw_gap_physics(params=None):
    """메인 Gap 물리 현상도 (상단) + 열저항 네트워크 (하단)"""

    fig = plt.figure(figsize=(24, 16), facecolor=BG)

    # ── 상단: Gap 물리 개략도 ──────────────────────────────────
    ax = fig.add_axes([0.02, 0.38, 0.96, 0.58])
    ax.set_xlim(-1, 21)
    ax.set_ylim(-1, 11)
    ax.set_facecolor(BG)
    ax.axis('off')

    # ====== 증발기 블록 ======
    evap_x, evap_w = 1.0, 3.5
    evap_y, evap_h = 1.0, 8.5
    evap = FancyBboxPatch((evap_x, evap_y), evap_w, evap_h,
            boxstyle="round,pad=0.15", fc=CYAN+"18", ec=CYAN, lw=2.5)
    ax.add_patch(evap)
    # 핀 표현
    for y in np.linspace(evap_y+0.5, evap_y+evap_h-0.5, 12):
        ax.plot([evap_x+0.2, evap_x+evap_w-0.2], [y, y],
                color=CYAN, lw=0.8, alpha=0.3)
    # 튜브 표현
    for cy in np.linspace(evap_y+1.2, evap_y+evap_h-1.2, 5):
        for cx in [evap_x+1.0, evap_x+2.5]:
            c = Circle((cx, cy), 0.35, fc=CYAN+"33", ec=CYAN, lw=1.2)
            ax.add_patch(c)
    ax.text(evap_x+evap_w/2, evap_y+evap_h+0.25, 'EVAPORATOR',
            ha='center', fontsize=14, color=CYAN, fontweight='bold')
    ax.text(evap_x+evap_w/2, evap_y-0.35,
            '$T_{surf}$ = 5°C\n$T_{sat,evap}$',
            ha='center', fontsize=9, color=CYAN, alpha=0.8)

    # ====== 응축기 블록 ======
    cond_x = 14.5
    cond = FancyBboxPatch((cond_x, evap_y), evap_w, evap_h,
            boxstyle="round,pad=0.15", fc=RED+"18", ec=RED, lw=2.5)
    ax.add_patch(cond)
    for y in np.linspace(evap_y+0.5, evap_y+evap_h-0.5, 12):
        ax.plot([cond_x+0.2, cond_x+evap_w-0.2], [y, y],
                color=RED, lw=0.8, alpha=0.3)
    for cy in np.linspace(evap_y+1.2, evap_y+evap_h-1.2, 5):
        for cx in [cond_x+1.0, cond_x+2.5]:
            c = Circle((cx, cy), 0.35, fc=RED+"33", ec=RED, lw=1.2)
            ax.add_patch(c)
    ax.text(cond_x+evap_w/2, evap_y+evap_h+0.25, 'CONDENSER',
            ha='center', fontsize=14, color=RED, fontweight='bold')
    ax.text(cond_x+evap_w/2, evap_y-0.35,
            '$T_{surf}$ = 45°C\n$T_{sat,cond}$',
            ha='center', fontsize=9, color=RED, alpha=0.8)

    # ====== GAP 영역 ======
    gap_x1 = evap_x + evap_w
    gap_x2 = cond_x
    gap_cx = (gap_x1 + gap_x2) / 2

    # Gap 경계 점선
    ax.plot([gap_x1, gap_x1], [0.3, 10.2], color=BORD, ls='--', lw=1, alpha=0.5)
    ax.plot([gap_x2, gap_x2], [0.3, 10.2], color=BORD, ls='--', lw=1, alpha=0.5)
    ax.annotate('', xy=(gap_x2-0.1, 10.3), xytext=(gap_x1+0.1, 10.3),
                arrowprops=dict(arrowstyle='<->', color=YEL, lw=2))
    ax.text(gap_cx, 10.55, 'Gap G [mm]', ha='center', fontsize=12,
            color=YEL, fontweight='bold')

    # Gap 배경
    gap_bg = FancyBboxPatch((gap_x1+0.1, evap_y+0.1), gap_x2-gap_x1-0.2, evap_h-0.2,
            boxstyle="round,pad=0.1", fc=YEL+"06", ec='none')
    ax.add_patch(gap_bg)

    # ====== ① 재순환 (Recirculation) ======
    # 상단 외부 재순환 화살표 (cond 출구 → evap 입구)
    style = "Simple,tail_width=1.5,head_width=8,head_length=6"
    # 응축기 출구 상방
    arr_recir = FancyArrowPatch((cond_x+evap_w+0.3, 8.0), (evap_x-0.3, 8.0),
        connectionstyle="arc3,rad=0.4", arrowstyle=style,
        fc=ORA, ec=ORA, lw=1.5, alpha=0.8)
    ax.add_patch(arr_recir)
    ax.text(gap_cx, 10.0, '① 외부 재순환 (Recirculation)',
            ha='center', fontsize=10, color=ORA, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc=BG, ec=ORA, alpha=0.8))
    ax.text(gap_cx, 9.4,
            '응축기 배출 고온 공기 → 증발기 입구 되돌아옴\n'
            '$\\Delta T_{recir}$ ∝ (1 - seal_fraction) / Gap',
            ha='center', fontsize=8, color=ORA, alpha=0.8)

    # ====== ② 혼합 (Mixing) ======
    # Gap 중앙부 와류
    for yc in [3.5, 5.5, 7.5]:
        theta = np.linspace(0, 2*np.pi, 50)
        r = 0.5 + 0.15*np.sin(3*theta)
        xs = gap_cx - 1.5 + r*np.cos(theta)
        ys = yc + r*np.sin(theta)*0.6
        ax.plot(xs, ys, color=PUR, lw=1, alpha=0.4)
        # 와류 화살표
        ax.annotate('', xy=(gap_cx-1.5+0.5, yc+0.35),
                    xytext=(gap_cx-1.5+0.3, yc+0.4),
                    arrowprops=dict(arrowstyle='->', color=PUR, lw=1.5))

    ax.text(gap_cx-1.5, 2.3, '② 내부 혼합',
            ha='center', fontsize=10, color=PUR, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc=BG, ec=PUR, alpha=0.8))
    ax.text(gap_cx-1.5, 1.7,
            'Gap 내 냉공기-고온공기\njat entrainment',
            ha='center', fontsize=7.5, color=PUR, alpha=0.8)

    # ====== 주 공기 흐름 ======
    # 증발기 입구 (좌 → 증발기)
    for y in [3, 5, 7]:
        ax.annotate('', xy=(evap_x, y), xytext=(-0.5, y),
                    arrowprops=dict(arrowstyle='->', color=CYAN, lw=2.5, alpha=0.7))
    ax.text(-0.7, 5.3, '입구 공기\n$T_{amb}$, RH',
            ha='center', fontsize=9, color=TEXT, alpha=0.8)

    # 증발기 → Gap
    for y in [3, 5, 7]:
        ax.annotate('', xy=(gap_x1+0.5, y), xytext=(gap_x1-0.3, y),
                    arrowprops=dict(arrowstyle='->', color=GRN, lw=2, alpha=0.6))

    # Gap → 응축기
    for y in [3, 5, 7]:
        ax.annotate('', xy=(gap_x2+0.3, y), xytext=(gap_x2-0.5, y),
                    arrowprops=dict(arrowstyle='->', color=YEL, lw=2, alpha=0.6))

    # 응축기 출구
    for y in [3, 5, 7]:
        ax.annotate('', xy=(cond_x+evap_w+1.0, y), xytext=(cond_x+evap_w+0.2, y),
                    arrowprops=dict(arrowstyle='->', color=RED, lw=2.5, alpha=0.7))
    ax.text(cond_x+evap_w+1.2, 5.3, '출구 공기\n$T_{out}$ > $T_{amb}$',
            ha='center', fontsize=9, color=TEXT, alpha=0.8)

    # ====== ③ 복사 (Radiation) ======
    # 파형 화살표 (cond → evap)
    xr = np.linspace(cond_x-0.2, gap_x1+0.2, 80)
    yr_base = 8.5
    yr = yr_base + 0.15*np.sin(12*np.pi*(xr-xr[0])/(xr[-1]-xr[0]))
    ax.plot(xr, yr, color=RED, lw=2, alpha=0.7,
            path_effects=_glow(RED, 1))
    ax.annotate('', xy=(gap_x1+0.2, yr_base), xytext=(gap_x1+1.0, yr_base),
                arrowprops=dict(arrowstyle='->', color=RED, lw=2))
    ax.text(gap_cx+1.5, 8.9, '③ 복사 $Q_{rad}$',
            ha='center', fontsize=10, color=RED, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc=BG, ec=RED, alpha=0.8))
    ax.text(gap_cx+1.5, 8.3,
            '$q_{rad} = \\sigma F_{12} (T_{cond}^4 - T_{evap}^4)$\n'
            'Gap↓ → $F_{12}$↑',
            ha='center', fontsize=7.5, color=RED, alpha=0.8)

    # ====== ④ 전도 (Conduction through frame) ======
    # 프레임 (상하단)
    frame_y_top = evap_y + evap_h
    frame_y_bot = evap_y
    for fy in [frame_y_bot, frame_y_top]:
        frame = FancyBboxPatch((gap_x1-0.2, fy-0.15), gap_x2-gap_x1+0.4, 0.3,
                boxstyle="round,pad=0.05", fc=DIM+"55", ec=DIM, lw=1.5)
        ax.add_patch(frame)
    ax.text(gap_cx, frame_y_top+0.5, '④ 프레임 전도 $Q_{cond}$',
            ha='center', fontsize=9, color=YEL, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.15', fc=BG, ec=YEL, alpha=0.7))
    # 전도 화살표
    ax.annotate('', xy=(gap_x1+0.5, frame_y_top+0.05),
                xytext=(gap_x2-0.5, frame_y_top+0.05),
                arrowprops=dict(arrowstyle='<->', color=YEL, lw=1.5, ls='--'))
    ax.text(gap_cx, frame_y_bot-0.5,
            'Frame\n$R_{cond} = L/k$',
            ha='center', fontsize=7.5, color=DIM, alpha=0.8)

    # ====== ⑤ 비말동반 (Carryover) ======
    # 물방울 궤적 (evap 표면 → gap → cond)
    np.random.seed(42)
    for i in range(8):
        y0 = np.random.uniform(2.5, 7.5)
        x_traj = np.linspace(gap_x1+0.2, gap_x1 + np.random.uniform(2, gap_x2-gap_x1-0.5), 20)
        y_traj = y0 + np.cumsum(np.random.normal(0, 0.04, 20)) - 0.02*np.arange(20)
        reached = x_traj[-1] > cond_x - 0.5
        clr = RED if reached else GRN
        ax.plot(x_traj, y_traj, color=clr, lw=0.8, alpha=0.5)
        # 물방울
        sz = np.random.uniform(15, 40)
        ax.scatter(x_traj[-1], y_traj[-1], s=sz, c=clr, alpha=0.7, marker='o', edgecolors='white', lw=0.3)

    ax.text(gap_cx+1.5, 2.3, '⑤ 비말동반 (Carryover)',
            ha='center', fontsize=10, color=CYAN, fontweight='bold',
            bbox=dict(boxstyle='round,pad=0.2', fc=BG, ec=CYAN, alpha=0.8))
    ax.text(gap_cx+1.5, 1.5,
            '응결수 액적 → Gap 비산 → 응축기 도달\n'
            'We > $We_{crit}$ → 이탈  |  P_reach ∝ 1/Gap\n'
            '응축기 재증발 → $Q_{carry}$ = $\\dot{m}_{carry}$ × $h_{fg}$',
            ha='center', fontsize=7.5, color=CYAN, alpha=0.8)

    # 증발기 표면 응결수 표시
    for y in [2.5, 4.0, 5.5, 7.0, 8.2]:
        ax.scatter(gap_x1+0.05, y, s=25, c='#66ccff', alpha=0.8, marker='o')
    ax.text(gap_x1+0.6, 1.3, '응결수\n$q_{cond}$',
            fontsize=7, color='#66ccff', ha='center', alpha=0.8)

    # ====== ⑥ 드레인 ======
    # 하단 드레인
    drain_y = evap_y + 0.15
    for dx in np.linspace(gap_x1+0.5, gap_x2-0.5, 5):
        ax.annotate('', xy=(dx, drain_y-0.5), xytext=(dx, drain_y+0.3),
                    arrowprops=dict(arrowstyle='->', color='#66ccff', lw=1, alpha=0.5))
    ax.text(gap_cx, 0.1, '드레인 (중력 배수)', fontsize=8, color='#66ccff',
            ha='center', alpha=0.6)

    # ====== 온도 프로파일 (우측) ======
    ax_t = fig.add_axes([0.88, 0.42, 0.10, 0.50])
    ax_t.set_facecolor(BG)
    temps = np.array([5, 15, 25, 35, 45])
    labels = ['$T_{evap}$=5°C', '', '$T_{amb}$=25°C', '', '$T_{cond}$=45°C']
    colors_t = [CYAN, DIM, GRN, DIM, RED]
    for t, l, c in zip(temps, labels, colors_t):
        ax_t.barh(t, 1, height=3, color=c, alpha=0.3)
        if l:
            ax_t.text(0.5, t, l, ha='center', fontsize=8, color=c, fontweight='bold')
    ax_t.set_xlim(0, 1)
    ax_t.set_ylim(0, 50)
    ax_t.set_ylabel('Temperature [°C]', color=TEXT, fontsize=9)
    ax_t.tick_params(colors=TEXT, labelsize=8)
    ax_t.set_xticks([])
    for s in ax_t.spines.values(): s.set_color(BORD)
    ax_t.set_title('T profile', color=TEXT, fontsize=9)

    # ── 하단: 열저항 네트워크 ──────────────────────────────────
    ax2 = fig.add_axes([0.04, 0.04, 0.92, 0.30])
    ax2.set_xlim(0, 20)
    ax2.set_ylim(0, 6)
    ax2.set_facecolor(PANEL)
    ax2.axis('off')

    ax2.text(10, 5.7, '열저항 네트워크 & 에너지 밸런스',
             ha='center', fontsize=13, color=WHITE, fontweight='bold')

    # 노드
    nodes = {
        'T_evap':  (1.5, 3.5, '$T_{evap}$\n5°C', CYAN),
        'T_in':    (4.5, 3.5, '$T_{in,eff}$', GRN),
        'Q_evap':  (7.5, 3.5, '$Q_{evap}$\nε-NTU', CYAN),
        'Q_net':   (10.5, 3.5, '$Q_{net}$', GRN),
        'T_cond':  (18.5, 3.5, '$T_{cond}$\n45°C', RED),
    }
    # Loss nodes (branch down)
    loss_nodes = {
        'Q_rad':   (13.0, 1.2, '$Q_{rad}$\n복사', RED),
        'Q_cond':  (15.0, 1.2, '$Q_{cond}$\n전도', YEL),
        'Q_mix':   (17.0, 1.2, '$Q_{mix}$\n혼합', PUR),
        'Q_carry': (19.0, 1.2, '$Q_{carry}$\n비말동반', CYAN),
    }
    # Recirculation (branch up from T_in)
    recir_nodes = {
        'dT_recir': (4.5, 5.5, '$\\Delta T_{recir}$\n재순환', ORA),
    }

    def _draw_node(ax, x, y, label, color, r=0.6):
        c = Circle((x, y), r, fc=color+'22', ec=color, lw=2)
        ax.add_patch(c)
        ax.text(x, y, label, ha='center', va='center', fontsize=7.5,
                color=color, fontweight='bold')

    def _draw_arrow(ax, x1, y1, x2, y2, color, label='', lw=1.5):
        ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                    arrowprops=dict(arrowstyle='->', color=color, lw=lw))
        if label:
            mx, my = (x1+x2)/2, (y1+y2)/2
            ax.text(mx, my+0.3, label, ha='center', fontsize=7, color=color)

    for k, (x, y, l, c) in nodes.items():
        _draw_node(ax2, x, y, l, c)
    for k, (x, y, l, c) in loss_nodes.items():
        _draw_node(ax2, x, y, l, c, r=0.5)
    for k, (x, y, l, c) in recir_nodes.items():
        _draw_node(ax2, x, y, l, c, r=0.5)

    # 주 흐름 연결
    _draw_arrow(ax2, 2.1, 3.5, 3.9, 3.5, TEXT, 'ε-NTU')
    _draw_arrow(ax2, 5.1, 3.5, 6.9, 3.5, TEXT, '$C_{air}(T_{in}-T_{evap})$')
    _draw_arrow(ax2, 8.1, 3.5, 9.9, 3.5, TEXT, '- losses')

    # 에너지 밸런스 수식
    ax2.text(10.5, 2.2,
             '$Q_{net} = Q_{useful} - Q_{rad} - Q_{cond} - Q_{mix} - Q_{carry}$',
             ha='center', fontsize=10, color=WHITE, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.3', fc=PANEL, ec=BORD, lw=1.5))

    # Loss 분기
    _draw_arrow(ax2, 10.5, 2.9, 13.0, 1.7, RED)
    _draw_arrow(ax2, 10.5, 2.9, 15.0, 1.7, YEL)
    _draw_arrow(ax2, 10.5, 2.9, 17.0, 1.7, PUR)
    _draw_arrow(ax2, 10.5, 2.9, 19.0, 1.7, CYAN)

    # 재순환 분기
    _draw_arrow(ax2, 4.5, 5.0, 4.5, 4.1, ORA)
    ax2.text(3.0, 5.5, 'open/semi\n→ $T_{in}$ ↑', fontsize=7, color=ORA, ha='center')

    # T_cond → Q_rad/Q_cond 연결
    _draw_arrow(ax2, 18.5, 2.9, 13.0, 1.7, RED, '')
    _draw_arrow(ax2, 18.5, 2.9, 15.0, 1.7, YEL, '')

    # Cap 정의
    ax2.text(10.5, 0.3,
             '$Cap = Q_{net} / Q_{ref} \\times 100\\%$    |    '
             '$Q_{ref}$: Gap → ∞ (손실 없는 이상 조건)',
             ha='center', fontsize=9, color=DIM)

    # Gap 의존성 요약 (좌측)
    ax2.text(1.0, 0.6, 'Gap↓ 효과:', fontsize=9, color=YEL, fontweight='bold')
    ax2.text(1.0, 0.15,
             '$\\Delta T_{recir}$↑  |  $Q_{rad}$↑  |  $Q_{cond}$↑  |  '
             '$Q_{mix}$↑  |  $P_{reach}$↑  |  $\\Delta P_{gap}$↑',
             fontsize=8, color=DIM)

    fig.suptitle('증발기-응축기 Gap 열유동 물리 현상도',
                 fontsize=18, color=WHITE, fontweight='bold', y=0.98)

    fig.savefig('gap_physics.png', dpi=150, bbox_inches='tight', facecolor=BG)
    plt.close(fig)
    print("  → gap_physics.png")


if __name__ == '__main__':
    draw_gap_physics()
