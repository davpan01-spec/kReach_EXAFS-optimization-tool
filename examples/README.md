# Example structures

Room-temperature crystal structures from the [Crystallography Open Database](https://www.crystallography.net/cod/) (COD). They are included so that the notebook runs without any external files. Each file keeps its original COD header. When you use one of them, cite the original article.

| File | Compound | Space group | COD ID | Source |
|---|---|---|---|---|
| `UO2_COD9009049.cif` | UO₂ (uraninite) | Fm-3m | [9009049](https://www.crystallography.net/cod/9009049.html) | R. W. G. Wyckoff, *Crystal Structures*, 2nd ed., vol. 1, p. 239 (Interscience, 1963), via AMCSD |
| `UF4_COD2310432.cif` | UF₄ | C2/c | [2310432](https://www.crystallography.net/cod/2310432.html) | A. C. Larson, R. B. Roof Jr., D. T. Cromer, *Acta Cryst.* **17**, 555 (1964) |
| `UCl4_COD1530722.cif` | UCl₄ | I4₁/amd | [1530722](https://www.crystallography.net/cod/1530722.html) | T. Schleid, G. Meyer, L. R. Morss, *J. Less-Common Met.* **132**, 69 (1987) |
| `Na2UO4_beta_COD1530200.cif` | β-Na₂UO₄ | Pccn | [1530200](https://www.crystallography.net/cod/1530200.html) | L. M. Kovba, *Radiokhimiya* **13**, 309 (1971) |

Notes:
- `UO2_COD9009049.cif` is the default `cif_file` of `EXAFS_OPTIMIZER.ipynb`.
- **Na₂UO₄:** COD contains only the β phase (quenched high-temperature phase), not α-Na₂UO₄, which is stable at room temperature. The Kovba entry is used because it is fully ordered. COD 1001339 (Gasperin 1985) uses partial occupancies, which FEFF cannot use directly.
- UF₄ has two U sites. Select one with `absorber_site`.
- FEFF writes its files to `<cif name>_feff/` next to the CIF. These folders are ignored by git.
