"""Computations, printouts and plots behind each section of EXAFS_OPTIMIZER.ipynb.

Each public method of Planner is one notebook section. A section computes what it needs
from the earlier ones automatically, so the sections can also be run out of order.
"""
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterMathtext
from scipy.stats import chi2, ncx2
from larch import Group
from larch.xafs import etok, ktoe, xftf

from .config import Config
from .style import C_SIG, C_NOISE, C_INK, C_MUTED, C_SCHEME, K_LABEL, panel_label, sci, log_xaxis
from .feff import simulate_feff, plot_feff


class Planner:
    """Holds a Config and the results of each section."""

    def __init__(self, cfg=None, **overrides):
        self.cfg = cfg if cfg is not None else Config()
        for key, value in overrides.items():
            if not hasattr(self.cfg, key):
                raise AttributeError(f"unknown control parameter: {key}")
            setattr(self.cfg, key, value)
        self.sim = None
        self._grid = self._opt = False

    # =====================================================================
    # computations shared between sections
    # =====================================================================
    def _ensure_sim(self, verbose=False):
        if self.sim is None:
            c = self.cfg
            self.sim, self.paths = simulate_feff(
                c.cif_file, c.absorber, edge=c.edge, absorber_site=c.absorber_site,
                cluster_size=c.cluster_size, R_cut=c.R_cut, use_feff8=c.use_feff8,
                s02=c.S02, e0=c.E0, sigma2=c.SIGMA2, verbose=verbose)
            self.kchi = self.sim.chi * self.sim.k**c.kweight      # k-weighted chi(k)

    def _ensure_grid(self):
        """Scan grid, detection-test constants and the chosen measurement of section 4."""
        self._ensure_sim()
        if self._grid:
            return
        c, sim = self.cfg, self.sim
        # scan grid: E0 + E_start to K_scan_max in steps of dk_exp
        self.k_start = etok(c.E_start)                # k = 0.5123 sqrt(E) Å^-1; 2.81 Å^-1 for 30 eV
        n_pts = int(np.floor((c.K_scan_max - self.k_start) / c.dk_exp + 1e-9)) + 1
        self.karr = self.k_start + c.dk_exp * np.arange(n_pts)        # last point <= K_scan_max
        self.npoints = self.karr.size
        self.in_scan = (sim.k >= self.k_start) & (sim.k <= c.K_scan_max)  # no data outside the scan
        # detection test
        self.M = max(int(round(c.win / c.dk_exp)), 1)                # measured points per window
        self.T_thr = chi2.ppf(1 - c.alpha, self.M)
        self.dk_sim = sim.k[1] - sim.k[0]
        self.n_win = max(int(round(c.win / self.dk_sim)), 1)          # simulation points per window
        # chosen measurement
        self.N_counts = c.count_rate * c.t_point
        self.sigma = 1/np.sqrt(self.N_counts)/c.edge_step             # 1-sigma noise in chi per point
        self.p_detect = self.detection_probability(self.N_counts)
        self.k_resolve = self.k_res_of(self.N_counts)
        self.t_total = self.npoints * c.t_point
        self._grid = True

    def detection_probability(self, N_k, es=None):
        """P(detect) of the window centred on each sim.k, for N_k counts per point (scalar or per sim.k).
        The sum over measured points is computed on the simulation grid, rescaled by dk_sim/dk_exp."""
        es = self.cfg.edge_step if es is None else es
        snr2 = np.where(self.in_scan, self.sim.chi**2 * N_k * es**2, 0)   # chi^2 / sigma^2
        lam = np.convolve(snr2, np.ones(self.n_win), mode="same") * (self.dk_sim / self.cfg.dk_exp)
        return ncx2.sf(self.T_thr, self.M, lam)

    def k_res_of(self, N_k, es=None):
        """Largest k in the scan whose window is detected with probability >= power (nan if none)."""
        ok = (self.detection_probability(N_k, es) >= self.cfg.power) & self.in_scan
        return self.sim.k[np.where(ok)[0][-1]] if ok.any() else np.nan

    def k_res_for(self, t0, n):
        """k_res for t(k) = t0 (k/k_start)^n."""
        return self.k_res_of(self.cfg.count_rate * t0 * (self.sim.k / self.k_start)**n)

    def total_time(self, t0, n):
        """Counting time (s) over the EXAFS scan for t(k) = t0 (k/k_start)^n."""
        return (t0 * (self.karr / self.k_start)**n).sum()

    def k_res_fine(self, N_k):
        """k_res interpolated linearly between the last resolved simulation point and the next one."""
        p = self.detection_probability(N_k)
        ok = (p >= self.cfg.power) & self.in_scan
        if not ok.any():
            return np.nan
        i = np.where(ok)[0][-1]
        if i + 1 < p.size and self.in_scan[i + 1]:
            return self.sim.k[i] + self.dk_sim * (p[i] - self.cfg.power) / (p[i] - p[i + 1])
        return self.sim.k[i]                             # resolved up to the end of the scan

    def t0_for_target(self, k_tgt, n, lo=1e-6, hi=1e6, n_iter=60):
        """Smallest t0 (s) with k_res >= k_tgt for exponent n (nan if not reachable up to hi)."""
        if not self.k_res_for(hi, n) >= k_tgt - 1e-9:
            return np.nan
        for _ in range(n_iter):                          # bisection in log(t0)
            mid = np.sqrt(lo * hi)
            if self.k_res_for(mid, n) >= k_tgt - 1e-9:
                hi = mid
            else:
                lo = mid
        return hi

    def _ensure_optimum(self):
        """Both optimisation problems of section 8."""
        self._ensure_grid()
        if self._opt:
            return
        c, sim = self.cfg, self.sim
        assert c.k_target <= c.K_scan_max, "k_target must lie inside the scan (<= K_scan_max)"
        # floor on t0 (s): signal/noise snr_lowk per point in the first window of the scan
        self.chi_lowk = np.sqrt(np.mean(sim.chi[self.in_scan & (sim.k <= self.k_start + c.win)]**2))
        self.t0_floor = ((c.snr_lowk / (self.chi_lowk * c.edge_step))**2 / c.count_rate
                         if c.t0_min == "auto" else c.t0_min)
        # (a) target range -> minimum total time
        t0_a = np.array([self.t0_for_target(c.k_target, n) for n in c.n_grid])
        self.t0_a = np.maximum(t0_a, self.t0_floor)      # floor on the low-k counting time (nan stays nan)
        self.time_a = np.array([self.total_time(t0, n) for t0, n in zip(self.t0_a, c.n_grid)])
        self.ia = np.nanargmin(self.time_a)
        self.i0 = np.argmin(np.abs(c.n_grid))            # n = 0 reference (constant time per point)
        # (b) time budget -> maximum k_res
        t0_b = c.T_budget * 3600 / np.array([((self.karr / self.k_start)**n).sum() for n in c.n_grid])
        self.t0_b = np.where(t0_b >= self.t0_floor, t0_b, np.nan)  # too steep: low k below the floor
        self.kres_b = np.array([self.k_res_for(t0, n) for t0, n in zip(self.t0_b, c.n_grid)])
        if np.all(np.isnan(self.kres_b)):                # every n gives t0 < floor: budget too small
            raise ValueError(
                f"T_budget = {c.T_budget:g} h is below the t0 floor even at n = 0 "
                f"(needs >= {self.t0_floor * self.karr.size / 3600:.3g} h): raise T_budget or lower snr_lowk / t0_min")
        self.best_b = np.isclose(self.kres_b, np.nanmax(self.kres_b))
        self.ib = np.where(self.best_b)[0][0]
        self._opt = True

    # =====================================================================
    # section 3
    # =====================================================================
    def feff_simulation(self):
        """FEFF run, list of paths, chi(k) and |chi(R)|."""
        c = self.cfg
        self.sim = None; self._grid = self._opt = False   # re-run FEFF and everything after it
        self._ensure_sim(verbose=True)
        plot_feff(self.sim, self.paths, c.kmin, c.kmax, c.kweight, c.R_cut,
                  title=f"FEFF simulation · {os.path.basename(c.cif_file)} · {c.absorber} {c.edge}-edge · "
                        rf"{len(self.paths)} paths < {c.R_cut:g} Å · $\sigma^2$ = {c.SIGMA2} Å$^2$")

    # =====================================================================
    # section 4
    # =====================================================================
    def counting_statistics(self):
        """Noise of the chosen measurement, detection probability and k_res."""
        self._ensure_grid()
        c, sim, kchi = self.cfg, self.sim, self.kchi
        k_start, k_resolve, sigma, N_counts = self.k_start, self.k_resolve, self.sigma, self.N_counts
        t_total = self.t_total
        print(f"scan {k_start:.2f}-{c.K_scan_max} Å^-1 (E0 + {c.E_start:g} eV), {self.npoints} points of {c.dk_exp} Å^-1")
        print(f"{c.count_rate:.3g} counts/s x {c.t_point:g} s/point = {N_counts:.1e} counts/point, sigma_chi = {sigma:.2e}")
        print(f"total counting time: {t_total:.0f} s = {t_total/60:.1f} min = {t_total/3600:.2f} h")
        print(f"max resolvable k (P_detect >= {c.power}, false alarm {c.alpha}, {self.M} points per window): {k_resolve:.2f} Å^-1")

        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.5, 6), sharex=True,
                                       gridspec_kw={"height_ratios": [1.6, 1], "hspace": 0.08})
        noise = sigma * sim.k**c.kweight
        for ax in (ax1, ax2):
            ax.axvspan(0, k_start, color=C_MUTED, alpha=0.10, lw=0)
            ax.axvspan(c.K_scan_max, sim.k.max(), color=C_MUTED, alpha=0.10, lw=0)
            if not np.isnan(k_resolve):
                ax.axvline(k_resolve, color=C_INK, ls="--", lw=1)
        if not np.isnan(k_resolve):
            ax1.text(k_resolve, 0.97, rf"  $k_{{res}}$ = {k_resolve:.2f} Å$^{{-1}}$",
                     transform=ax1.get_xaxis_transform(), ha="left", va="top", fontsize=9.5, fontweight="bold")

        ax1.fill_between(sim.k, noise, -noise, color=C_NOISE, alpha=0.25, lw=0, label=r"noise, $\pm1\sigma$")
        ax1.plot(sim.k, kchi, color=C_SIG, lw=1.3, label="simulated EXAFS")
        ax1.axhline(0, color=C_MUTED, lw=0.6)
        lim = 1.15 * max(np.abs(kchi).max(), noise[sim.k <= np.nan_to_num(k_resolve)].max(initial=0))
        ax1.set_ylim(-lim, lim); ax1.set_xlim(0, sim.k.max())
        ax1.set_ylabel(rf"$k^{{{c.kweight}}}\chi(k)$ (Å$^{{-{c.kweight}}}$)")
        ax1.legend(loc="lower left")
        ax1.set_title(f"{sci(c.count_rate)} counts/s × {c.t_point:g} s = {sci(N_counts)} counts/pt · "
                      rf"edge step = {c.edge_step:g} · $\sigma_\chi$ = {sci(sigma)} · $\Delta k$ = {c.dk_exp:g} Å$^{{-1}}$"
                      "\nshaded: outside the scan", fontsize=9, loc="left")
        panel_label(ax1, "(a)")

        ax2.plot(sim.k, np.where(self.in_scan, self.p_detect, np.nan), color=C_INK, lw=1.4)
        ax2.axhline(c.power, color=C_MUTED, ls=":", lw=1)
        ax2.text(0.99, c.power - 0.02, f"required power {c.power:g}", transform=ax2.get_yaxis_transform(),
                 ha="right", va="top", fontsize=8.5, color=C_MUTED)
        ax2.set_ylim(-0.04, 1.2); ax2.set_yticks([0, 0.25, 0.5, 0.75, 1])
        ax2.set_xlabel(K_LABEL)
        ax2.set_ylabel("detection probability")
        ax2.text(0.99, 0.55, f"{c.win:g} Å$^{{-1}}$ window ({self.M} pts)\nfalse-alarm rate α = {c.alpha:g}",
                 transform=ax2.transAxes, ha="right", va="center", fontsize=8.5, color=C_MUTED)
        panel_label(ax2, "(b)")
        plt.show()

    # =====================================================================
    # section 5
    # =====================================================================
    def kres_vs_time(self):
        """k_res vs constant counting time per point, one curve per edge step."""
        self._ensure_grid()
        c = self.cfg
        # edge step is an ordered quantity -> one hue, light (small step) to dark (large step)
        colors = [C_SIG] if len(c.edge_steps) == 1 else plt.cm.Blues(np.linspace(0.35, 1.0, len(c.edge_steps)))

        fig, ax = plt.subplots(figsize=(7.5, 4.6))
        for es, col in zip(c.edge_steps, colors):
            kr = np.array([self.k_res_of(c.count_rate * t, es) for t in c.T_range])
            ax.step(c.T_range, kr, where="mid", color=col, lw=1.8, label=f"{es:g}")
        ax.plot(c.t_point, self.k_resolve, "o", color=C_INK, ms=5, zorder=3)
        ax.annotate(f"  section 4: {c.t_point:g} s/pt", (c.t_point, self.k_resolve), fontsize=8.5, va="top")

        ax.axhline(c.K_scan_max, color=C_MUTED, ls="--", lw=1)
        ax.text(c.T_range[0] * 1.3, c.K_scan_max, "end of scan", va="bottom", fontsize=8.5, color=C_MUTED)
        log_xaxis(ax)
        ax.set_xlim(c.T_range[0], c.T_range[-1]); ax.set_ylim(0, c.K_scan_max * 1.1)
        ax.set_xlabel("counting time per data point, $t$ (s)")
        ax.set_ylabel(r"max. resolvable $k$ (Å$^{-1}$)")
        ax.grid(True, which="major", color=C_MUTED, alpha=0.25, lw=0.6)
        top = ax.secondary_xaxis("top", functions=(lambda t: t * c.count_rate, lambda N: N / c.count_rate))
        top.xaxis.set_major_formatter(LogFormatterMathtext())
        top.set_xlabel(f"counts per data point, $N$ (at {sci(c.count_rate)} counts/s)")
        ax.tick_params(axis="x", which="both", top=False)
        if len(c.edge_steps) > 1:                        # one series: the note names it, no legend
            ax.legend(title="edge step", loc="upper left", bbox_to_anchor=(0, 0.92), title_fontsize=9)
        es_txt = f"edge step = {c.edge_steps[0]:g}\n" if len(c.edge_steps) == 1 else ""
        ax.text(0.98, 0.04, es_txt + rf"P$_{{detect}}\geq$ {c.power:g}, α = {c.alpha:g}" + "\n"
                rf"{c.win:g} Å$^{{-1}}$ window, $\Delta k$ = {c.dk_exp:g} Å$^{{-1}}$",
                transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5, color=C_MUTED)
        plt.show()

    # =====================================================================
    # section 6
    # =====================================================================
    def r_space(self):
        """|chi(R)| over the full range vs the range up to k_res, noiseless and noisy."""
        self._ensure_grid()
        c, sim, kchi, k_resolve = self.cfg, self.sim, self.kchi, self.k_resolve
        kmin, kmax, kweight = c.kmin, c.kmax, c.kweight
        if np.isnan(k_resolve) or k_resolve < kmin + 2:
            print(f"k_res = {k_resolve:.2f} Å^-1 leaves less than 2 Å^-1 above kmin = {kmin}: no meaningful FT")
            return
        k_hi = min(k_resolve, kmax)
        rng = np.random.default_rng(c.seed)
        sigma_sim = self.sigma * np.sqrt(c.dk_exp / self.dk_sim)   # same information content on the sim grid
        ft = {"theory": (sim.chi, kmax),
              "noiseless": (sim.chi, k_hi),
              "noisy": (sim.chi + np.where(self.in_scan, rng.normal(0, sigma_sim, sim.k.size), 0), k_hi)}
        groups = {}
        for name, (chi_, k2) in ft.items():
            g = Group(k=sim.k, chi=chi_)
            xftf(g, kmin=kmin, kmax=k2, kweight=kweight, dk=1, window="hanning")
            groups[name] = g
        dR = {name: np.pi / (2 * (k2 - kmin)) for name, (_, k2) in ft.items()}
        print(f"R resolution pi/(2 dk): {dR['theory']:.3f} Å for k = [{kmin}, {kmax}], "
              f"{dR['noiseless']:.3f} Å for k = [{kmin}, {k_hi:.2f}]")

        style = {"theory": dict(color=C_SIG, lw=1.6, label=f"theory, $k$ = {kmin:g}–{kmax:g} Å$^{{-1}}$"),
                 "noiseless": dict(color=C_NOISE, lw=1.6, label=f"$k$ = {kmin:g}–{k_hi:.1f} Å$^{{-1}}$ ($k_{{res}}$)"),
                 "noisy": dict(color=C_MUTED, lw=0.9, label=f"same range, with noise ({c.t_point:g} s/pt)")}
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"wspace": 0.28})
        ax1.plot(sim.k, np.where(self.in_scan, groups["noisy"].chi * sim.k**kweight, np.nan), **style["noisy"])
        ax1.plot(sim.k, kchi, color=C_SIG, lw=1.3)
        for name, ls in (("theory", "--"), ("noiseless", "-")):
            ax1.plot(sim.k, groups[name].kwin * np.abs(kchi).max() * 1.1, color=style[name]["color"], lw=0.9, ls=ls)
        ax1.axhline(0, color=C_MUTED, lw=0.6)
        ax1.set_xlim(0, sim.k.max()); ax1.set_ylim(-1.3*np.abs(kchi).max(), 1.3*np.abs(kchi).max())
        ax1.set_xlabel(K_LABEL); ax1.set_ylabel(rf"$k^{{{kweight}}}\chi(k)$ (Å$^{{-{kweight}}}$)")
        ax1.text(0.98, 0.03, "lines on top: FT windows", transform=ax1.transAxes, ha="right", va="bottom",
                 fontsize=8, color=C_MUTED)
        panel_label(ax1, "(a)")

        for name in ("noisy", "theory", "noiseless"):
            ax2.plot(groups[name].r, groups[name].chir_mag, **style[name])
        ax2.set_xlim(0, c.R_cut + 1); ax2.set_ylim(0, None)
        ax2.set_xlabel(r"$R$ (Å)"); ax2.set_ylabel(rf"$|\chi(R)|$ (Å$^{{-{kweight+1}}}$)")
        ax2.legend(loc="upper right")
        panel_label(ax2, "(b)")
        fig.suptitle(rf"$|\chi(R)|$ with data up to $k_{{res}}$ · {sci(self.N_counts)} counts/pt · "
                     rf"$\delta R$ = {dR['theory']:.2f} → {dR['noiseless']:.2f} Å", fontsize=10, y=1.0)
        plt.show()

    # =====================================================================
    # section 7
    # =====================================================================
    def kweighted_time(self):
        """k_res vs t0 and vs total time for t(k) = t0 (k/k_start)^n, n in time_exps."""
        self._ensure_grid()
        c, k_start = self.cfg, self.k_start
        kres = {n: np.array([self.k_res_for(t0, n) for t0 in c.T0_range]) for n in c.time_exps}
        ttot = {n: np.array([self.total_time(t0, n) for t0 in c.T0_range]) for n in c.time_exps}

        print(f"{c.count_rate:.3g} counts/s, edge step = {c.edge_step}, scan {k_start:.2f}-{c.K_scan_max} Å^-1 "
              f"(E0 + {c.E_start:g} eV) in {c.dk_exp} steps, {self.npoints} points")
        print(f"{'T0 (s)':>8}" + "".join(f"{'n='+str(n)+' k_res':>12}{'time (h)':>10}" for n in c.time_exps))
        for t0 in [1e-3, 1e-2, 1e-1, 1, 10, 100]:
            i = np.argmin(np.abs(c.T0_range - t0))
            print(f"{t0:>8g}" + "".join(f"{kres[n][i]:>12.2f}{ttot[n][i]/3600:>10.3g}" for n in c.time_exps))

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True, gridspec_kw={"wspace": 0.08})
        for n, col in zip(c.time_exps, C_SCHEME):
            ax1.step(c.T0_range, kres[n], where="mid", color=col, lw=1.8, label=rf"$n$ = {n}")
            ax2.step(ttot[n]/3600, kres[n], where="mid", color=col, lw=1.8, label=rf"$n$ = {n}")
        for ax in (ax1, ax2):
            log_xaxis(ax)
            ax.axhline(c.K_scan_max, color=C_MUTED, ls="--", lw=1)
            ax.grid(True, which="major", color=C_MUTED, alpha=0.25, lw=0.6)
        ax2.text(0.97, c.K_scan_max, "end of scan", ha="right", transform=ax2.get_yaxis_transform(), va="bottom", fontsize=8.5, color=C_MUTED)
        ax1.set_xlim(c.T0_range[0], c.T0_range[-1]); ax1.set_ylim(0, c.K_scan_max*1.08)
        ax1.set_xlabel(r"counting time per point at $k_{start}$, $t_0$ (s)")
        ax1.set_ylabel(r"max. resolvable $k$ (Å$^{-1}$)")
        ax1.legend(title=r"$t(k) = t_0\,(k/k_{start})^{n}$", loc="lower right", title_fontsize=9)
        panel_label(ax1, "(a)")
        ax2.set_xlabel("total scan time (h)")
        ax2.text(0.98, 0.04, "same data, plotted vs total time:\nfair comparison of the schemes",
                 transform=ax2.transAxes, ha="right", va="bottom", fontsize=8.5, color=C_MUTED)
        panel_label(ax2, "(b)")
        fig.suptitle(rf"$R_0$ = {sci(c.count_rate)} counts/s · edge step = {c.edge_step:g} · $k_{{start}}$ = {k_start:.2f} Å$^{{-1}}$ ($E_0$+{c.E_start:g} eV) · "
                     rf"$\Delta k$ = {c.dk_exp:g} Å$^{{-1}}$ · P$_{{detect}}\geq$ {c.power:g}, α = {c.alpha:g}",
                     fontsize=10, y=0.99)
        plt.show()

    # =====================================================================
    # section 8
    # =====================================================================
    def optimal_exponent(self):
        """(a) fastest n for k_target, (b) largest k_res for T_budget."""
        self._ensure_optimum()
        c, n_grid = self.cfg, self.cfg.n_grid
        t0_a, time_a, ia, i0 = self.t0_a, self.time_a, self.ia, self.i0
        t0_b, kres_b, best_b, ib = self.t0_b, self.kres_b, self.best_b, self.ib
        print(f"t0 floor: {self.t0_floor:.3g} s at k_start, summed over all scans"
              + (f" (auto: rms chi = {self.chi_lowk:.3g} in {self.k_start:.2f}-{self.k_start + c.win:.2f} Å^-1, S/N = {c.snr_lowk:g})"
                 if c.t0_min == "auto" else ""))
        print(f"(a) k_res >= {c.k_target:g} Å^-1: fastest n = {n_grid[ia]:.2f}, t0 = {t0_a[ia]:.3g} s, "
              f"total {time_a[ia]/3600:.3g} h  (n = 0: {time_a[i0]/3600:.3g} h, "
              f"{100*(1 - time_a[ia]/time_a[i0]):.0f}% saved)")
        nb_txt = (f"n = {n_grid[best_b].min():.2f}-{n_grid[best_b].max():.2f} (all equal)" if best_b.sum() > 1
                  else f"n = {n_grid[ib]:.2f}")
        print(f"(b) {c.T_budget:g} h budget: best {nb_txt}, t0 = {t0_b[ib]:.3g} s, k_res = {kres_b[ib]:.2f} Å^-1  "
              f"(n = 0: {kres_b[i0]:.2f} Å^-1)")
        for lbl, i in (("(a)", ia), ("(b)", ib)):
            if i in (0, n_grid.size - 1):
                print(f"    {lbl}: the optimum n = {n_grid[i]:.2f} is at the edge of n_grid: widen n_grid or raise t0_min")
        if np.nanmax(kres_b) >= c.K_scan_max - self.dk_sim:
            print("    k_res reaches the end of the scan: extend K_scan_max or reduce T_budget to discriminate")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"wspace": 0.3})
        ax1.plot(n_grid, time_a/3600, color=C_SIG, lw=1.6)
        ax1.plot(n_grid[ia], time_a[ia]/3600, "*", color=C_NOISE, ms=12, zorder=3,
                 label=f"optimum $n$ = {n_grid[ia]:.2f}: {time_a[ia]/3600:.3g} h")
        ax1.set_ylim(0, 1.1*np.nanmax(time_a)/3600)
        ax1.set_xlabel(r"time exponent $n$"); ax1.set_ylabel("total scan time (h)")
        ax1.set_title(rf"(a) time to reach $k_{{res}}$ = {c.k_target:g} Å$^{{-1}}$", fontsize=9.5, loc="left")
        ax1.legend(loc="lower left")
        ax1.grid(True, which="major", color=C_MUTED, alpha=0.25, lw=0.6)

        ax2.step(n_grid, kres_b, where="mid", color=C_SIG, lw=1.6)
        nb_lo, nb_hi = n_grid[best_b].min(), n_grid[best_b].max()
        ax2.axvspan(nb_lo - 0.025, nb_hi + 0.025, color=C_NOISE, alpha=0.15, lw=0,
                    label=rf"optimum: $k_{{res}}$ = {kres_b[ib]:.2f} Å$^{{-1}}$, $n$ = {nb_lo:.2f}–{nb_hi:.2f}"
                          if best_b.sum() > 1 else rf"optimum: $k_{{res}}$ = {kres_b[ib]:.2f} Å$^{{-1}}$, $n$ = {nb_lo:.2f}")
        ax2.axhline(c.K_scan_max, color=C_MUTED, ls="--", lw=1)
        ax2.text(0.02, c.K_scan_max, "end of scan", transform=ax2.get_yaxis_transform(), va="bottom",
                 fontsize=8.5, color=C_MUTED)
        ax2.set_ylim(np.nanmin(kres_b) - 1, c.K_scan_max + 0.8)
        ax2.set_xlabel(r"time exponent $n$"); ax2.set_ylabel(r"max. resolvable $k$ (Å$^{-1}$)")
        ax2.set_title(rf"(b) $k_{{res}}$ for a {c.T_budget:g} h budget", fontsize=9.5, loc="left")
        ax2.legend(loc="lower center")
        ax2.grid(True, which="major", color=C_MUTED, alpha=0.25, lw=0.6)
        fig.suptitle(rf"$t(k) = t_0\,(k/k_{{start}})^n$ · {sci(c.count_rate)} counts/s · edge step = {c.edge_step:g} · "
                     rf"{self.npoints} points, {self.k_start:.2f}–{c.K_scan_max:g} Å$^{{-1}}$", fontsize=10, y=1.02)
        plt.show()

    # =====================================================================
    # section 9
    # =====================================================================
    def split_scans(self):
        """Instrument settings for n_scans scans, dead time, and a simulated measurement."""
        self._ensure_optimum()
        c, sim, n_grid = self.cfg, self.sim, self.cfg.n_grid
        k_start, karr, npoints, in_scan = self.k_start, self.karr, self.npoints, self.in_scan
        count_rate, n_scans, kweight = c.count_rate, c.n_scans, c.kweight

        # ---- scheme from section 8 ----
        if c.plan_from == "a":
            n_p, t0_p = n_grid[self.ia], self.t0_a[self.ia]
        else:                                            # 8(b): ties on the grid broken with the interpolated k_res
            tied = np.where(self.best_b)[0]
            i_p = tied[np.argmax([self.k_res_fine(count_rate * self.t0_b[i] * (sim.k / k_start)**n_grid[i]) for i in tied])]
            n_p, t0_p = n_grid[i_p], self.t0_b[i_p]
        w = (karr / k_start)**n_p                        # t(k_i) / t0 on the scan grid
        t0_tot = c.margin * t0_p                         # t0 summed over all scans
        k_plan = self.k_res_fine(count_rate * t0_tot * (sim.k / k_start)**n_p)

        def plan(ns):
            """(t0 per scan, counting time, dead time) in s for the measurement split into ns scans."""
            return t0_tot / ns, t0_tot * w.sum(), ns * (npoints * c.t_dead_point + c.t_overhead_scan)

        t0_scan, t_count, t_dead = plan(n_scans)
        print(f"section 8({c.plan_from}): n = {n_p:.2f}, t0 = {t0_p:.3g} s (floor {self.t0_floor:.3g} s); x margin {c.margin:g} -> t0 = {t0_tot:.3g} s in total, "
              f"k_res = {k_plan:.2f} Å^-1")
        if c.plan_from == "a" and c.K_scan_max < c.k_target + c.win / 2:
            print(f"  note: K_scan_max < k_target + win/2 = {c.k_target + c.win/2:g} Å^-1: the last windows are truncated")
        print(); print("instrument settings, per scan:")
        print(f"  scan range:  k = {k_start:.2f} (E0 + {c.E_start:g} eV) to K_scan_max = {c.K_scan_max:g} Å^-1 "
              f"(E0 + {ktoe(c.K_scan_max):.0f} eV), dk = {c.dk_exp:g} Å^-1, {npoints} points")
        print(f"  time/point:  t(k) = {t0_scan:.3g} s x (k/{k_start:.2f})^{n_p:.2f}, "
              f"from {t0_scan:.3g} s at k_start to {t0_scan*w[-1]:.3g} s at {karr[-1]:.2f} Å^-1")
        print(f"  scans:       {n_scans}, each {t_count/n_scans/3600:.2f} h counting + {t_dead/n_scans/60:.1f} min dead time "
              f"= {(t_count + t_dead)/n_scans/3600:.2f} h")
        print(f"  total:       {t_count/3600:.1f} h counting + {t_dead/3600:.1f} h dead time = {(t_count + t_dead)/3600:.1f} h")
        print(f"  analysis:    use k up to k_res = {k_plan:.2f} Å^-1 (FT/fit kmax)")
        k_useful = min(k_plan + c.win / 2, c.K_scan_max)  # data beyond this is not resolved
        f_waste = w[karr > k_useful].sum() / w.sum()
        if f_waste > 0.05:
            print(f"  WARNING: {100*f_waste:.0f}% of the counting time ({f_waste*t_count/3600:.1f} h) goes to k > {k_useful:.2f} Å^-1, "
                  f"beyond k_res + win/2: set K_scan_max to about {k_useful:.2f} Å^-1 and re-run sections 4-9")
        if t0_scan < c.t_dead_point:
            print(f"  note: at low k the counting time per point ({t0_scan:.2g} s) is shorter than the dead time "
                  f"({c.t_dead_point:g} s): fewer, longer scans waste less time")

        opts = sorted({1, 2, 5, 10, 20, 50, n_scans})
        print(f"\n{'scans':>6}{'t0/scan (s)':>13}{'h per scan':>12}{'counting (h)':>14}{'dead (h)':>10}{'total (h)':>11}")
        for m in opts:
            a, b, d = plan(m)
            print(f"{m:>6}{a:>13.3g}{(b + d)/m/3600:>12.2f}{b/3600:>14.1f}{d/3600:>10.1f}{(b + d)/3600:>11.1f}"
                  + ("   <- n_scans" if m == n_scans else ""))

        # ---- simulated measurement: n_scans noisy scans and their average ----
        tk = t0_tot / n_scans * (np.maximum(sim.k, k_start) / k_start)**n_p        # s per point, one scan
        sig_pt = np.where(in_scan, 1 / np.sqrt(count_rate * tk) / c.edge_step, 0)  # sigma_chi per measured point
        sig_sim = sig_pt * np.sqrt(c.dk_exp / self.dk_sim)                         # same information on the sim grid
        rng = np.random.default_rng(c.seed)
        chi_scans = np.where(in_scan, sim.chi + rng.normal(size=(n_scans, sim.k.size)) * sig_sim, 0)
        chi_avg = chi_scans.mean(axis=0)
        chi_true = np.where(in_scan, sim.chi, 0)

        k_lo, k_hi = max(c.kmin, k_start), min(k_plan, c.K_scan_max)              # FT over the resolved range
        def ft(chi_):
            g = Group(k=sim.k, chi=chi_)
            xftf(g, kmin=k_lo, kmax=k_hi, kweight=kweight, dk=1, window="hanning")
            return g
        g_true, g_avg, g_one = ft(chi_true), ft(chi_avg), ft(chi_scans[0])
        g_noise = ft(chi_avg - chi_true)                                           # noise of the average alone
        sig_avg = sig_pt / np.sqrt(n_scans)                                        # final noise per point in chi
        noise_R = np.sqrt(np.mean(g_noise.chir_mag[g_noise.r <= c.R_cut + 1]**2))  # rms noise level in |chi(R)|
        print(f"\nfinal noise per point in chi (average of {n_scans} scans): {sig_avg[in_scan][0]:.2g} at k_start, "
              f"{np.interp(k_hi, sim.k, sig_avg):.2g} at k_res = {k_hi:.2f} Å^-1 (one scan: x{np.sqrt(n_scans):.2g}); "
              f"rms noise in |chi(R)|: {noise_R:.2g}, peak of |chi(R)|: {g_true.chir_mag.max():.3g}")

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.5, 3.9), gridspec_kw={"wspace": 0.28})
        kw = sim.k**kweight
        band = np.where(in_scan, sig_avg * kw, np.nan)
        ax1.fill_between(sim.k, -band, band, color=C_NOISE, alpha=0.25, lw=0, zorder=0,
                         label=rf"final noise $\pm1\sigma$: $\sigma_\chi$ = {sig_avg[in_scan][0]:.1e} → {np.interp(k_hi, sim.k, sig_avg):.1e}")
        ax1.plot(sim.k, np.where(in_scan, chi_scans[0] * kw, np.nan), color=C_MUTED, lw=0.6, alpha=0.8, label="one scan")
        ax1.plot(sim.k, np.where(in_scan, chi_avg * kw, np.nan), color=C_NOISE, lw=1.0, label=f"average of {n_scans} scans")
        ax1.plot(sim.k, np.where(in_scan, sim.chi * kw, np.nan), color=C_SIG, lw=1.3, label="noiseless")
        ax1.axvline(k_hi, color=C_INK, ls="--", lw=1)
        ax1.text(k_hi, 0.97, rf"  $k_{{res}}$ = {k_hi:.2f}", transform=ax1.get_xaxis_transform(), ha="left", va="top", fontsize=8.5)
        lim = 1.4 * np.abs(sim.chi * kw)[in_scan].max()
        ax1.set_xlim(0, c.K_scan_max + 0.5); ax1.set_ylim(-lim, lim)
        ax1.axhline(0, color=C_MUTED, lw=0.6)
        ax1.set_xlabel(K_LABEL); ax1.set_ylabel(rf"$k^{{{kweight}}}\chi(k)$ (Å$^{{-{kweight}}}$)")
        ax1.legend(loc="lower left", fontsize=8)
        panel_label(ax1, "(a)")

        ax2.plot(g_one.r, g_one.chir_mag, color=C_MUTED, lw=0.9, label="one scan")
        ax2.plot(g_avg.r, g_avg.chir_mag, color=C_NOISE, lw=1.4, label=f"average of {n_scans} scans")
        ax2.plot(g_true.r, g_true.chir_mag, color=C_SIG, lw=1.4, ls="--", label="noiseless")
        ax2.fill_between(g_noise.r, 0, g_noise.chir_mag, color=C_NOISE, alpha=0.25, lw=0,
                         label=f"noise of the average (rms {noise_R:.2g})")
        ax2.set_xlim(0, c.R_cut + 1); ax2.set_ylim(0, None)
        ax2.set_xlabel(r"$R$ (Å)"); ax2.set_ylabel(rf"$|\chi(R)|$ (Å$^{{-{kweight+1}}}$)")
        ax2.legend(loc="upper right", fontsize=8)
        panel_label(ax2, "(b)")
        fig.suptitle(rf"{n_scans} scans of $t(k)$ = {t0_tot/n_scans:.3g} s $(k/k_{{start}})^{{{n_p:.2f}}}$ · "
                     rf"FT {k_lo:.2f}–{k_hi:.2f} Å$^{{-1}}$, $k^{{{kweight}}}$, Hanning, dk = 1 · seed {c.seed}", fontsize=10, y=1.0)
        plt.show()
