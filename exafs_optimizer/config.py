"""All user inputs of EXAFS_OPTIMIZER.ipynb in one object."""
from dataclasses import dataclass, field
import numpy as np


@dataclass
class Config:
    """Control parameters. The section numbers refer to EXAFS_OPTIMIZER.ipynb."""

    # ---- structure and FEFF simulation (section 3) ----
    cif_file: str = "examples/UO2_COD9009049.cif"   # relative to the notebook, or absolute
    absorber: str = "U"            # absorbing element symbol
    edge: str = "L3"               # absorption edge
    absorber_site: int = 1         # crystallographic site of the absorber (if several)
    cluster_size: float = 9.0      # FEFF cluster radius (Å), must be > R_cut
    R_cut: float = 7               # keep paths with Reff < R_cut (Å)
    use_feff8: bool = True         # True: feff8l, False: feff6l
    S02: float = 1                 # path parameters applied to all paths
    E0: float = 0.0
    SIGMA2: float = 0.005          # Debye-Waller factor (Å^2)
    kmin: float = 3.0              # FT range (Å^-1)
    kmax: float = 15.0
    kweight: int = 2               # k-weight, used everywhere

    # ---- measurement (sections 4-8) ----
    count_rate: float = 1.6e3      # counts/s per data point, background subtracted
    t_point: float = 30            # counting time per data point (s)
    edge_step: float = 1.25        # edge jump (transmission) or edge fraction of counts (fluorescence)
    E_start: float = 30.0          # start of the EXAFS scan above the edge (eV)
    K_scan_max: float = 13         # end of the scan (Å^-1), at most the simulated grid end (20 Å^-1)
    dk_exp: float = 0.05           # k-spacing of the measured data (Å^-1)

    # ---- detection test (sections 4-8) ----
    win: float = 1.5               # analysis window width in k (Å^-1)
    alpha: float = 0.003           # false-alarm probability
    power: float = 0.95            # required detection probability per window

    # ---- k_res vs counting time per point (section 5) ----
    T_range: np.ndarray = field(default_factory=lambda: np.logspace(-2, 4, 121))   # s per point
    edge_steps: list = field(default_factory=lambda: [1])                         # one curve each

    # ---- R-space comparison and simulated scans (sections 6, 9) ----
    seed: int = 0                  # random seed of the noisy realisations

    # ---- k-weighted counting time (sections 7-8) ----
    time_exps: list = field(default_factory=lambda: [0, 1, 2, 3])                 # n in t0 (k/k_start)^n
    T0_range: np.ndarray = field(default_factory=lambda: np.logspace(-3, 2, 101))  # t0 (s)
    n_grid: np.ndarray = field(default_factory=lambda: np.arange(0, 4.001, 0.01))  # exponents scanned
    k_target: float = 12.0         # (a) required k_res (Å^-1)
    T_budget: float = 3.5          # (b) total counting time available (h)
    t0_min: object = "auto"        # floor on t0 at k_start (s), summed over all scans, or "auto"
    snr_lowk: float = 5            # "auto": rms(chi)/sigma_chi per point in the first window

    # ---- splitting into scans (section 9) ----
    plan_from: str = "a"           # "a": plan of section 8(a), "b": of section 8(b)
    margin: float = 1.25           # safety factor on the counting time
    n_scans: int = 20              # number of scans (averaged afterwards)
    t_dead_point: float = 1        # dead time per point in every scan (s)
    t_overhead_scan: float = 120   # dead time per scan (s)
