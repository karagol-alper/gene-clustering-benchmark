"""Shared figure style: validated categorical palette, recessive chrome."""
import matplotlib as mpl

INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#898781"
GRID = "#e1e0d9"
BASE = "#c3c2b7"
SURFACE = "#ffffff"
# categorical slots in fixed order (validated palette, light mode)
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = (
    "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948")
SERIES = [BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED]
PATH, BEN = ORANGE, BLUE            # class identity used in every figure
SEQ = ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]

MM = 1 / 25.4
FULL_W = 174 * MM                    # full journal width


def apply():
    mpl.rcParams.update({
        "font.family": "Liberation Sans", "font.size": 7.5,
        "axes.titlesize": 8.5, "axes.titleweight": "bold", "axes.titlelocation": "left",
        "axes.labelsize": 7.5, "axes.labelcolor": INK, "axes.edgecolor": BASE,
        "axes.linewidth": 0.8, "axes.grid": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "xtick.color": INK2, "ytick.color": INK2, "xtick.labelsize": 7, "ytick.labelsize": 7,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.major.width": 0.6,
        "ytick.major.width": 0.6, "legend.fontsize": 6.8, "legend.frameon": False,
        "figure.dpi": 150, "savefig.dpi": 600, "savefig.bbox": "tight",
        "savefig.pad_inches": 0.03, "lines.linewidth": 1.6, "lines.markersize": 4.5,
        "pdf.fonttype": 42, "ps.fonttype": 42, "text.color": INK,
    })


def panel(ax, letter, title, x=-0.02, y=1.0):
    ax.set_title(f"{letter}  {title}", loc="left", pad=6, fontsize=8.5, fontweight="bold",
                 color=INK, x=x, y=y)


def hgrid(ax):
    ax.yaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)


def vgrid(ax):
    ax.xaxis.grid(True, color=GRID, lw=0.6)
    ax.set_axisbelow(True)
