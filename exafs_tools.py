"""Compatibility module: the helpers now live in the exafs_optimizer package.

Kept so that `from exafs_tools import ...` (used by EXAFS_MC_FIT.ipynb) keeps working.
"""
from exafs_optimizer import __version__
from exafs_optimizer.style import (C_SIG, C_NOISE, C_INK, C_MUTED, C_SCHEME, K_LABEL,
                                   set_style, panel_label, sci, log_xaxis)
from exafs_optimizer.feff import simulate_feff, plot_feff
