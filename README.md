# EXAFS Optimizer

**Version 0.2.1-alpha** (test release, see [CHANGELOG.md](CHANGELOG.md))

Plan the counting statistics of an EXAFS measurement before beamtime. From a CIF file, a count rate and a counting-time scheme, the tool estimates:

- how far in $k$ the signal can be resolved above counting noise;
- what the Fourier transform looks like with only that $k$ range;
- which $k$-weighted counting time $t(k)=t_0(k/k_{start})^n$ reaches a target $k$ fastest, or reaches the largest $k$ within a time budget;
- how to split that plan into scans, and what the split costs in dead time.

A companion notebook for the Monte Carlo uncertainty of the fitted structural parameters vs $k_{max}$ is *coming in a later release*.

The method, an example and the limitations are described in [DOCUMENTATION.md](DOCUMENTATION.md).

## Files
| File | Content |
|---|---|
| `EXAFS_OPTIMIZER.ipynb` | Main workflow: control parameters, section descriptions and outputs (FEFF simulation, detection test, counting-time optimisation, splitting into scans) |
| `exafs_optimizer/` | Backend: `config.py` (all parameters and defaults), `planner.py` (computations and plots, one method per notebook section), `feff.py` (FEFF simulation from a CIF), `style.py` (plot style) |
| `exafs_tools.py` | Compatibility module that re-exports the FEFF and plot helpers for older notebooks |
| `examples/` | Public example structures from COD (UO₂, UF₄, UCl₄, β-Na₂UO₄), see [examples/README.md](examples/README.md) |
| `DOCUMENTATION.md` | Theory, example results, advantages and limitations |
| `requirements.txt` | Tested package versions |
| `EXAFS_MC_FIT.ipynb` | Monte Carlo fit: parameter uncertainty vs $k_{max}$ (*coming in a later release*, not included) |

## Installation
```
pip install -r requirements.txt
```
This installs the tested versions, including `xraylarch` 2026.3.1. `xraylarch` includes FEFF6/FEFF8 and pymatgen. Python 3.13 was used for testing.

## Usage
1. Open `EXAFS_OPTIMIZER.ipynb`. By default it uses `examples/UO2_COD9009049.cif`. To use your own structure, set `cif_file` to its path, either relative to the notebook or absolute.
2. Edit the **Control parameters** cell (absorber, edge, count rate, time per point, scan range, …).
3. Run the notebook from the top. FEFF writes its files to `<cif name>_feff/` next to the CIF.

The backend can also be used without the notebook:
```python
from exafs_optimizer import Planner
plan = Planner(cif_file="examples/UO2_COD9009049.cif", count_rate=1.6e3, k_target=12)
plan.optimal_exponent()       # section 8; any parameter not given keeps its default
```

## Citation
If you use this code, please cite: *[paper reference to be added]*. If you use one of the example structures, also cite its original article (listed in [examples/README.md](examples/README.md)).

## License
MIT, see [LICENSE](LICENSE). The example CIF files come from the Crystallography Open Database and keep their own terms (see their headers).
