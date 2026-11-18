# Changelog

## 0.2.1-alpha (test release)

### Added
- `examples/`: public room-temperature structures from COD (UO₂, UF₄, UCl₄, β-Na₂UO₄). The notebook runs out of the box with `examples/UO2_COD9009049.cif`.
- `requirements.txt` with the tested versions (xraylarch 2026.3.1, Python 3.13).
- MIT `LICENSE`.
- `exafs_optimizer.__version__`, printed by the setup cell.
- Section 9: splitting the plan into scans, with dead time per point and per scan.
- Simulated noise floor of the average of `n_scans` scans.

### Changed
- All code moved out of the notebook into the `exafs_optimizer/` package (`Config` for the parameters, `Planner` with one method per section). The notebook now contains only the control parameters, the section descriptions and the outputs. Each section computes what it needs from earlier ones, so sections can also be run out of order. `exafs_tools.py` remains as a compatibility module.
- `DOCUMENTATION.md` §2 now describes the default notebook setup (UO₂), so testers can compare their numbers directly.
- The default `cif_file` is a relative path to the example structure.
- The notebook is shipped without saved outputs.

### Fixed
- The count rate was printed with `:.0e`, so 1.6×10³ appeared as `2e+03`. It now prints with 3 significant figures.
- The scan grid could end just above `K_scan_max` (13.01 instead of 13 Å⁻¹), counting the point towards the time but not the signal. The last point is now ≤ `K_scan_max`.
- The section 3 plot title showed the full CIF path. It now shows the file name only.
- A minimizer bug in the helper.

### Not included
- `EXAFS_MC_FIT.ipynb` (Monte Carlo parameter uncertainty) is coming in a later release.

### What to check as a tester
1. **Installation.** Does `pip install -r requirements.txt` work on your system (OS, Python version)? Does FEFF8 run?
2. **Default run.** Run `EXAFS_OPTIMIZER.ipynb` from the top without changes. Do your printed numbers match DOCUMENTATION.md §2 (k_res = 10.35 Å⁻¹ at 30 s/point; n = 2.84 and 3.1 h for k_target = 12 Å⁻¹)?
3. **Your own structure.** Try a CIF and edge you know from experiment. Is k_res plausible compared with your real data? Is it roughly an upper bound, as expected?
4. **Other examples.** UF₄ has two U sites (`absorber_site`). Does the FEFF step work for all four example files?
5. **Clarity.** Are the control parameters, printouts and plots understandable without asking? Is anything in DOCUMENTATION.md unclear or wrong?
6. **Edge cases.** Very low count rates, a `k_target` close to `K_scan_max`, a `T_budget` that is too small. Do you get clear messages instead of errors?

Please report problems with the OS, Python version and the full error message.
