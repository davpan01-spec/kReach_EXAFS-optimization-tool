"""Plot style shared by all notebooks."""
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext

C_SIG, C_NOISE, C_INK, C_MUTED = "#2a78d6", "#eb6834", "#1f1f1e", "#8a8984"
C_SCHEME = ["#2a78d6", "#eb6834", "#1baf7a"]          # categorical slots 1-3 (validated all-pairs)
K_LABEL = r"$k$ (Å$^{-1}$)"


def set_style():
    plt.rcParams.update({
        "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
        "font.family": "sans-serif", "font.size": 10, "mathtext.fontset": "dejavusans",
        "axes.linewidth": 0.8, "axes.labelsize": 10.5, "axes.titlesize": 10.5,
        "axes.edgecolor": C_INK, "axes.labelcolor": C_INK, "text.color": C_INK,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 4, "ytick.major.size": 4, "xtick.minor.size": 2, "ytick.minor.size": 2,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.color": C_INK, "ytick.color": C_INK,
        "lines.linewidth": 1.5, "legend.frameon": False, "legend.fontsize": 9,
        "axes.grid": False,
    })


def panel_label(ax, s):
    ax.text(0.015, 0.97, s, transform=ax.transAxes, fontweight="bold", va="top", ha="left")


def sci(x):
    """1e4 -> $10^{4}$, 2.5e4 -> $2.5\\times10^{4}$"""
    m, e = f"{x:.1e}".split("e"); m = m.rstrip("0").rstrip(".")
    return rf"$10^{{{int(e)}}}$" if m == "1" else rf"${m}\times10^{{{int(e)}}}$"


def log_xaxis(ax):
    ax.set_xscale("log")
    ax.xaxis.set_major_locator(LogLocator(base=10))
    ax.xaxis.set_major_formatter(LogFormatterMathtext())
