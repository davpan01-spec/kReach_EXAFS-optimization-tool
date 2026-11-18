"""FEFF simulation from a CIF, and the plot of the simulated signal."""
import os, glob, warnings
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import hilbert
from larch import Group
from larch.xrd.structure2feff import structure2feffinp
from larch.xafs import feff8l, feff6l, feffpath, ff2chi, xftf

from .style import C_SIG, C_INK, C_MUTED, K_LABEL, panel_label

warnings.filterwarnings("ignore", category=DeprecationWarning, module=r"larch\.")  # pymatgen dict-interface noise

def simulate_feff(cif_file, absorber, edge="K", absorber_site=1, cluster_size=7.0, R_cut=5.0,
                  use_feff8=True, s02=1.0, e0=0.0, sigma2=0.005, verbose=True):
    """CIF -> feff.inp -> FEFF run -> paths with Reff < R_cut -> summed chi(k).

    FEFF runs in '<cif name>_feff' next to the CIF. Returns (sim, paths): sim is a larch
    Group with k and chi on larch's default grid (dk = 0.05 Å^-1), paths the list of FeffPath.
    """
    if not os.path.isfile(cif_file):
        raise FileNotFoundError(f"CIF file not found: {cif_file}")

    feff_folder = os.path.splitext(os.path.abspath(cif_file))[0] + "_feff"
    os.makedirs(feff_folder, exist_ok=True)
    with open(cif_file, "r") as f:
        cif_text = f.read()            # pass text: larch's filename branch mis-detects the .cif format
    feff_inp = structure2feffinp(cif_text, absorber, edge=edge, cluster_size=cluster_size,
                                 absorber_site=absorber_site, version8=use_feff8)
    if feff_inp.startswith("# could not read"):
        raise ValueError(f"pymatgen could not parse CIF: {cif_file}")
    with open(os.path.join(feff_folder, "feff.inp"), "w") as f:
        f.write(feff_inp)

    runner = feff8l if use_feff8 else feff6l
    runner(folder=feff_folder, feffinp="feff.inp", verbose=False)

    paths = []
    for fname in sorted(glob.glob(os.path.join(feff_folder, "feff*.dat"))):
        path = feffpath(fname, s02=s02, e0=e0, sigma2=sigma2)
        if path.reff < R_cut:
            paths.append(path)
    if not paths:
        raise RuntimeError(f"No FEFF paths with Reff < {R_cut} Å in {feff_folder}")

    if verbose:
        print(f"{len(paths)} paths with Reff < {R_cut} Å:")
        print(f"{'file':<14}{'Reff(Å)':>9}{'degen':>8}{'nleg':>6}  geometry")
        for path in paths:
            geom = " -> ".join(g[0] for g in path.geom) + f" -> {path.geom[0][0]}"
            print(f"{os.path.basename(path.filename):<14}{path.reff:>9.4f}{path.degen:>8.1f}{path.nleg:>6d}  {geom}")

    sim = Group()
    ff2chi(paths, group=sim)
    return sim, paths


def plot_feff(sim, paths, kmin, kmax, kweight, R_cut, title=""):
    """(a) k-weighted chi(k) with envelope and FT window, (b) chi(R) with path Reff ticks."""
    xftf(sim, kmin=kmin, kmax=kmax, kweight=kweight, dk=1, window="hanning")
    kchi = sim.chi * sim.k**kweight
    kchi_env = np.abs(hilbert(kchi))                   # amplitude envelope
    win_k = sim.kwin * kchi_env.max() * 1.15          # FT window, scaled for display
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 3.8), gridspec_kw={"wspace": 0.28})

    ax1.fill_between(sim.k, kchi_env, -kchi_env, color=C_SIG, alpha=0.13, lw=0, label="amplitude envelope")
    ax1.plot(sim.k, win_k, color=C_MUTED, lw=0.9, ls="--", label=f"FT window ({kmin:g}–{kmax:g} Å$^{{-1}}$)")
    ax1.plot(sim.k, kchi, color=C_SIG, lw=1.3, label=r"$\Sigma$ paths")
    ax1.axhline(0, color=C_MUTED, lw=0.6)
    ax1.set_xlim(0, sim.k.max())
    ax1.set_xlabel(K_LABEL); ax1.set_ylabel(rf"$k^{{{kweight}}}\chi(k)$ (Å$^{{-{kweight}}}$)")
    ax1.legend(loc="lower right")
    panel_label(ax1, "(a)")

    ax2.plot(sim.r, sim.chir_re, color=C_MUTED, lw=0.9, label=r"Re[$\chi(R)$]")
    ax2.plot(sim.r, sim.chir_mag, color=C_SIG, lw=1.6, label=r"$|\chi(R)|$")
    ax2.axhline(0, color=C_MUTED, lw=0.6)
    ymax = sim.chir_mag[sim.r <= R_cut + 1].max()
    for path in paths:                                 # Reff of each included path
        ax2.plot([path.reff, path.reff], [-0.12*ymax, -0.04*ymax], color=C_INK, lw=1)
    ax2.axvline(R_cut, color=C_INK, ls=":", lw=1)
    ax2.text(R_cut - 0.05, 1.08*ymax, rf"$R_{{cut}}$ = {R_cut:g} Å", ha="right", va="top", fontsize=9)
    ax2.set_xlim(0, R_cut + 1); ax2.set_ylim(-1.0*ymax, 1.12*ymax)
    ax2.set_xlabel(r"$R$ (Å)"); ax2.set_ylabel(rf"$\chi(R)$ (Å$^{{-{kweight+1}}}$)")
    ax2.legend(loc="upper right", bbox_to_anchor=(1, 0.9))
    ax2.text(0.985, 0.03, "ticks: path $R_{eff}$\n(no phase correction)", transform=ax2.transAxes,
             ha="right", va="bottom", fontsize=8, color=C_MUTED)
    panel_label(ax2, "(b)")
    fig.suptitle(title, fontsize=10, y=1.0)
    plt.show()
