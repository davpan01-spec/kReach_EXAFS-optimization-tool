"""EXAFS Optimizer: plan the counting statistics of an EXAFS measurement before beamtime.

    from exafs_optimizer import Config, Planner
    plan = Planner(Config(cif_file="examples/UO2_COD9009049.cif", count_rate=1.6e3))
    plan.feff_simulation()        # section 3
    plan.counting_statistics()    # section 4
    ...
"""
__version__ = "0.2.1-alpha"

from .config import Config
from .planner import Planner
from .feff import simulate_feff, plot_feff
from .style import (C_SIG, C_NOISE, C_INK, C_MUTED, C_SCHEME, K_LABEL,
                    set_style, panel_label, sci, log_xaxis)
