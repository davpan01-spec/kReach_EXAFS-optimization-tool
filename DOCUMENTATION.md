# EXAFS Optimizer: method, advantages and limitations

EXAFS Optimizer is a pre-beamtime planning tool. For a known or assumed crystal structure, it estimates the largest photoelectron wavenumber $k_{res}$ at which the EXAFS signal can still be distinguished from counting noise, for a given count rate and counting-time scheme. It also finds the $k$-weighted counting-time scheme that reaches a target $k$ in the shortest time. A companion Monte Carlo notebook, which will give the resulting uncertainties of the fitted structural parameters, is coming in a later release (§1.6). The code uses Larch [1], FEFF [2,3] and pymatgen [4].

## 1. Theory

### 1.1 Simulated signal
The CIF is converted to a FEFF input with pymatgen and Larch. FEFF6 or FEFF8 then computes the scattering paths, and all paths with $R_{eff} < R_{cut}$ are summed in the standard EXAFS equation:

$$\chi(k)=\sum_j \frac{S_0^2 N_j f_j(k)}{k R_j^2}\, e^{-2\sigma^2 k^2} e^{-2R_j/\lambda(k)} \sin\!\left[2kR_j+\delta_j(k)\right],$$

using a common $S_0^2$, $\Delta E_0$ and Debye–Waller factor $\sigma^2$ for all paths. The result is $\chi(k)$ on a grid with $\delta k_{sim}$ = 0.05 Å⁻¹.

### 1.2 Noise model
The measurement provides a count rate $R_0$ (counts s⁻¹ per data point, background subtracted) and a counting time $t$ per point, so each point collects $N = R_0 t$ counts. With Poisson statistics, the noise in the absorption $\mu x$ is approximately $1/\sqrt{N}$. Normalising by the edge step $\Delta\mu x$ gives the noise in $\chi$:

$$\sigma_\chi = \frac{1}{\sqrt{N}\,\Delta\mu x}.$$

In fluorescence, $\Delta\mu x$ is replaced by the fraction of detected counts that come from the absorber edge. Noise is assumed to be Gaussian and uncorrelated between data points spaced by $\delta k_{exp}$ [5].

### 1.3 Detection test and $k_{res}$
For a window of width $w$ in $k$ (about one oscillation of the main shell) containing $M = w/\delta k_{exp}$ measured points, we use the energy-detector statistic [6]

$$T=\sum_{i=1}^{M}\frac{\chi_{obs}(k_i)^2}{\sigma_\chi^2}.$$

For pure noise, $T$ follows a $\chi^2_M$ distribution. With signal, it follows a non-central $\chi^2_M(\lambda)$ with $\lambda=\sum_i \chi(k_i)^2/\sigma_\chi^2$. The threshold $T_{thr}$ is the $(1-\alpha)$ quantile of $\chi^2_M$, where $\alpha$ is the false-alarm rate (default 0.003, i.e. 3σ). The detection probability, or power, is

$$P_D(k) = P\!\left[\chi^2_M(\lambda(k)) > T_{thr}\right],$$

evaluated for the window centred at each $k$. The sum $\lambda$ is computed on the finer simulation grid and rescaled by $\delta k_{sim}/\delta k_{exp}$. The window is counted as resolved when $P_D \geq P_{req}$ (default 0.95), and $k_{res}$ is the largest resolved $k$ inside the scan, $[k_{start}, k_{scan,max}]$. Because both the error types and the window are stated explicitly, $k_{res}$ is a reproducible, structure-specific criterion rather than a rule of thumb.

### 1.4 Counting-time schemes
Most of the EXAFS amplitude decays at high $k$, so the counting time per point is increased with $k$:

$$t(k)=t_0\left(\frac{k}{k_{start}}\right)^n,\qquad N(k)=R_0\,t(k),\qquad T_{tot}=\sum_i t(k_i).$$

$n = 0$ is constant time per point. With $n = 2w$, the noise in $k^w\chi(k)$ is constant over the scan. The detection test of §1.3 is applied with the $k$-dependent $\sigma_\chi(k)$. The notebook then solves two optimisation problems over $n$:

- **(a) Target range.** Minimise $T_{tot}$ subject to $k_{res} \ge k_{target}$. $k_{res}$ never decreases with $t_0$, so the smallest sufficient $t_0$ for each $n$ is found by bisection.
- **(b) Time budget.** Maximise $k_{res}$ subject to $T_{tot} = T_{budget}$, i.e. $t_0 = T_{budget}/\sum_i (k_i/k_{start})^n$.

**Low-$k$ floor.** $k_{res}$ tests only the high-$k$ end, so without a constraint both problems favour ever steeper $n$, moving all the time to high $k$. A floor $t_0 \ge t_{0,min}$ on the total counting time per point at $k_{start}$ keeps the low-$k$ data usable for fitting. By default (`"auto"`) it follows from the simulated signal: $t_{0,min} = \left(\mathrm{SNR}_{low} / (\chi_{rms}\,\Delta\mu x)\right)^2 / R_0$, where $\chi_{rms}$ is the rms of $\chi$ in the first window $[k_{start}, k_{start}+w]$. The noise per point at the start of the scan is then $\chi_{rms}/\mathrm{SNR}_{low}$ (notebook default $\mathrm{SNR}_{low}$ = 5, `snr_lowk`). In (a), $t_0$ is raised to the floor; in (b), exponents that would fall below it are excluded. The optimum then lies inside the scanned range of $n$.

**Splitting into scans.** The scheme of (a) or (b), with its counting time multiplied by a safety margin, is divided over $N_s$ averaged scans of $t_0/N_s$ each. With counting noise only, this gives the same counts per point, and so the same noise and $k_{res}$, as one scan with $t_0$. The split changes only the dead time: $N_s\,(N_{pts}\,t_{dead,pt} + t_{dead,scan})$. A margin is advisable because $k_{res}$ from (a) sits exactly at the detection limit, and $k_{res}$ jumps between groups of oscillations as the counting time changes.

### 1.5 Consequence in $R$ space
$\chi(k)$ is Fourier transformed (Hanning window, $k$-weight $w$) over the full range $[k_{min}, k_{max}]$ and over the resolvable range $[k_{min}, k_{res}]$, both noiseless and for one noisy realisation. Comparing the two transforms shows how the shorter range lowers the $R$ resolution, $\delta R \approx \pi/(2\Delta k)$, and how many independent points remain [7].

### 1.6 Monte Carlo parameter uncertainty (`EXAFS_MC_FIT.ipynb`, coming in a later release)
*This notebook is not part of release 0.2.1-alpha. The method is described here for completeness.*

Noisy copies of the simulated $\chi(k)$ are fitted in $R$ space with `feffit`, using the same paths and four global parameters: $S_0^2$, $\Delta E_0$, an isotropic expansion $\alpha$ with $\Delta R_j=\alpha R_{eff,j}$, and $\sigma^2$. This is repeated for several values of $k_{max}$. The standard deviation of the fitted values over $n_{MC}$ realisations is the actual 1σ uncertainty. It is compared with the uncertainty reported by `feffit`, and its own sampling error is $\approx s/\sqrt{2(n_{MC}-1)}$.

## 2. Example: UO₂, U L₃ edge
This is the default setup of `EXAFS_OPTIMIZER.ipynb`, using `examples/UO2_COD9009049.cif` and the default control parameters. It uses FEFF8 with 38 paths up to $R_{cut}$ = 7 Å, $\sigma^2$ = 0.005 Å², $R_0 = 1.6\times10^3$ counts s⁻¹, $\Delta\mu x$ = 1.25, and a scan from $E_0$+30 eV (2.81 Å⁻¹) to 13 Å⁻¹ in steps of 0.05 Å⁻¹ (204 points). The detection test uses a 1.5 Å⁻¹ window, α = 0.003 and $P_{req}$ = 0.95, and the low-$k$ floor uses $\mathrm{SNR}_{low}$ = 5.

- **Constant time.** At 30 s per point ($4.8\times10^4$ counts, 1.7 h), the signal is resolved up to $k_{res}$ = 10.35 Å⁻¹.
- **Low-$k$ floor.** $\chi_{rms}$ = 0.069 in the first window gives $t_{0,min}$ = 2.1 s.
- **Target range of 12 Å⁻¹.** The optimal exponent is $n$ = 2.84 with $t_0$ = 2.1 s, which needs 3.1 h of counting compared with 7.2 h at constant time (57 % less).
- **3.5 h budget.** $k_{res}$ = 12.35 Å⁻¹ for any $n$ between 2.81 and 2.94, compared with 10.80 Å⁻¹ at $n$ = 0.
- **Splitting into scans.** The plan of (a) with a margin of 1.25 gives $k_{res}$ = 12.43 Å⁻¹. Split into 20 scans with 1 s dead time per point and 2 min per scan, this takes 3.9 h of counting plus 1.8 h of dead time (5.7 h). Over 2–5 scans, the dead time would be only 0.2–0.5 h.

The results are deterministic except for the noisy realisations in sections 6 and 9 (random seed `seed`).

## 3. Advantages
- **Before the beamtime.** Counting statistics and scan strategy can be planned, and beamtime requested, from a structure file alone.
- **Specific to the structure.** The criterion uses the actual simulated signal, including its damping, scatterers and interference between shells, rather than a generic $k$ range.
- **Explicit statistics.** $k_{res}$ is defined by a stated false-alarm rate, required power and window, so the result is reproducible and can be compared between samples and beamlines.
- **Quantitative optimisation.** The tool returns the time exponent and $t_0$ for a target $k$ or a time budget, together with the time saved compared with constant counting.
- **Fast.** One evaluation of $k_{res}$ takes about 1 ms. Full scans over the parameters run in seconds. Only the Monte Carlo fit (later release) takes longer, about 1 min for 250 fits.
- **Two complementary views.** The detection test says *where* the signal is present. The Monte Carlo fit (later release) says *what it is worth* for the structural parameters.
- **Open and standard.** The code is built on Larch, FEFF and pymatgen, and all inputs are collected in one parameter cell.

## 4. Limitations
- **Counting noise only.** Systematic errors are not included: monochromator glitches, harmonics, sample inhomogeneity and thickness effects, self-absorption, detector dead time and non-linearity, $I_0$ noise, and errors from background subtraction. In practice these often limit $k$ before counting statistics do, so $k_{res}$ is an upper bound.
- **The model must be right.** The results depend on the assumed structure and, most strongly, on $\sigma^2$, which sets the high-$k$ damping and depends on temperature and disorder. A single global $\sigma^2$ is used for all paths. Paths beyond $R_{cut}$ and any FEFF inaccuracies in amplitude and phase are neglected.
- **Detection is not fittability.** A window can pass the test while still adding little information to a fit. The Monte Carlo notebook (later release) addresses this, but only for the four-parameter model with fixed $k$-weight.
- **Depends on the analysis choices.** $k_{res}$ depends on the window width, $\alpha$ and $P_{req}$. It is resolved only to the 0.05 Å⁻¹ simulation grid, and windows near the end of the scan are truncated.
- **Simplified noise.** Points are taken as independent (no oversampling or rebinning correlations). Fluorescence background is captured only through the effective edge step.
- **Counting time only.** Dead time per point (monochromator motion, readout) is not included. It adds the same constant for every $n$, so it does not change the optimum in (a), but it reduces the effective budget in (b). Section 9 adds the dead time per point and per scan when the measurement is split into scans.
- **Idealised Monte Carlo (later release).** The fit model is exact (no model error), the noise is white, and the time per point is constant.

## References
1. M. Newville, *J. Phys.: Conf. Ser.* **430**, 012007 (2013).
2. J. J. Rehr, J. J. Kas, F. D. Vila, M. P. Prange, K. Jorissen, *Phys. Chem. Chem. Phys.* **12**, 5503 (2010).
3. A. L. Ankudinov, B. Ravel, J. J. Rehr, S. D. Conradson, *Phys. Rev. B* **58**, 7565 (1998).
4. S. P. Ong *et al.*, *Comput. Mater. Sci.* **68**, 314 (2013).
5. M. Newville, B. I. Boyanov, D. E. Sayers, *J. Synchrotron Rad.* **6**, 264 (1999).
6. H. Urkowitz, *Proc. IEEE* **55**, 523 (1967).
7. E. A. Stern, *Phys. Rev. B* **48**, 9825 (1993).
