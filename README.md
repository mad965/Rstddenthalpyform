# Rstdenthalpyform

R package to retrieve the standard enthalpy of formation ΔfH° (and related quantities: Gibbs energy of formation ΔfG°, entropy S°, heat capacity Cp°) of chemical compounds from **three independent thermochemical databases**, with their **source**, **reference conditions**, and **provenance**.

> Version 0.1.0. The package provides access to the databases, a multi-source search function, and a linear regression model to predict ΔfH° from elemental composition. Missing values remain missing: no value is ever invented.

## Data sources

| Source | Records | ΔfH° | ΔfG° | S | Cp | Notes |
| --- | --- | --- | --- | --- | --- | --- |
| **ATcT** (Active Thermochemical Tables, v1.228) | 3 455 | ✅ | ❌ | ❌ | ❌ | Strict ΔfH°, mostly gas phase |
| **RMG** (Reaction Mechanism Generator database) | 24 920 | ✅ | ❌ | ✅ | ❌ | Consolidated entries, gas phase |
| **CHNOSZ** (OBIGT, via `subcrt()`) | 3 007 | ✅ | ✅ | ✅ | ✅ | All phases, computed via `subcrt()` |

The three sources are kept **separate** (no merging) in `data/ATcT.rda`, `data/RMG.rda`, and `data/CHNOSZ.rda`. Each file contains a data frame with the **same 18-column schema**, allowing uniform querying without mixing provenances.

### Units and conventions

- `delta_fH` and `delta_fG`: **kJ/mol**.
- `S` and `Cp`: **J/(mol·K)**.
- `phase`: `"g"` (gas), `"l"` (liquid), `"cr"` (crystal), `"aq"` (aqueous).
- ATcT and RMG provide a strict ΔfH°. CHNOSZ provides an **apparent** enthalpy of formation (Benson-Helgeson convention) computed by `subcrt()`: it is not strictly identical to ATcT/RMG, and the `source` column allows this to be traced.
- Missing values are kept as `NA` and never replaced by zero.

## Installation

```r
# from the package directory
install.packages("devtools")          # once
devtools::document()                  # generates help files (man/) from roxygen comments
devtools::install()
or, from the unzipped archive: devtools::install_local("Rstdenthalpyform").

Usage
Multi-source search
r
library(Rstdenthalpyform)

# All sources for methane (gas)
get_enthalpy("CH4", phase = "g")

# A single source
get_enthalpy("CH4", source = "ATcT", phase = "g")

# Partial, case-insensitive search
get_enthalpy("methane", exact = FALSE)
get_enthalpy() returns one row per observation found (source × phase × temperature), with the columns: source, name, formula, CAS, phase, compound_nature, delta_fH, delta_fG, S, Cp, temperature, uncertainty, unit, source_id, reference.

Predicting ΔfH° by linear regression
The package includes a linear model of the form:

ΔfH∘=∑XβX⋅nXΔ 
 
where
nX is the number of atoms of element 
X in the molecular formula. One model is fitted per source, using the elements most frequent in that source.

r
load("data/ATcT.rda"); load("data/RMG.rda"); load("data/CHNOSZ.rda")

models <- fit_all_models(atct_harmonized, rmg_harmonized, chnosz_harmonized)

# Coefficients (kJ/mol per atom)
compare_model_coefficients(models)

# Predict for a new compound
predict_enthalpy(models$ATcT, "C2H6")
Predicted vs. true values
r
plot_predicted_vs_true(models$ATcT, atct_harmonized, "ATcT")
Points are coloured by the period of the main element in the formula (periodic table row) and shaded by the physical phase. A dashed line shows the ideal 
y=x.

Code organisation
text
R/
  get_enthalpy.R          multi-source search (ATcT, RMG, CHNOSZ)
  prediction_model.R      regression model, prediction, plotting
  chemical_identification.R  formula parsing, species search
  import_data.R           OBIGT reading and normalisation
  families.R              chemical families
  database.R              database access, references, duplicates
  prediction.R            result object, related properties
  method_selection.R      ThermoPredictor, get_formation_enthalpy()
  validation.R            integrity and consistency checks
  units.R                 unit conversions
Data organisation
text
data/
  ATcT.rda        # 3 455 rows, 18 columns
  RMG.rda         # 24 920 rows, 18 columns
  CHNOSZ.rda      # 3 007 rows, 18 columns
inst/extdata/
  ATcT/           # consolidated CSV
  RMG/            # consolidated CSV
  OBIGT/          # CHNOSZ raw CSVs + consolidated CSV
data-raw/
  consolidate_atct.R
  consolidate_rmg.R
  consolidate_chnosz.R
  make_rda.R
  utils_harmonization.R
The scripts in data-raw/ (excluded from the installed package) document how each source was imported, cleaned, unit-converted, and classified. They are reproducible.

Chemical families
list_families() returns the requested families. Classification is rule-based, using the formula, charge, phase, and origin file. Two families are not determinable with these data and remain empty rather than guessed: organometallics (metal–carbon bond) and biomolecules (CHNOSZ's biotic file was merged into organic_aq). Radicals cannot be identified either: the "Ions and radicals" family contains only ions.

Data consistency check
check_gibbs_consistency() recomputes ΔfG° = ΔfH° − T (S − Σ S°elements) for every entry that has G, H, and S, using elemental entropies supplied by CHNOSZ. Deviations above the tolerance (5 kJ/mol by default) are flagged: they indicate entries to examine, not necessarily errors.

Known limitations
Liquid water: in OBIGT, the water entry (liquid phase) has no H value in the CSV. subcrt() nevertheless computes it, and get_enthalpy("H2O", phase = "l") returns the computed value. Ice (ice) and steam (steam) also have a value.

CHNOSZ: 54 species (aqueous arsenides, tripeptides) have no delta_fH computable by subcrt(); they remain NA.

Inconsistent entries: among 2 853 CHNOSZ entries with G, H, and S, 66 deviate by more than 5 kJ/mol from the recomputed ΔfG° = ΔfH° − T ΔfS° (see check_gibbs_consistency()). They are flagged, not corrected.

No uncertainties in OBIGT; ATcT and RMG provide them when available.

A single value per species in ATcT and CHNOSZ, several possible conformers in RMG: all rows are returned, and the source_id column distinguishes them.

Updating the data
See the scripts in data-raw/:

consolidate_atct.R: reads ATcT_enthalpies.csv, harmonises, writes ATcT_consolidated.csv.

consolidate_rmg.R: reads rmg_entries.csv, harmonises, writes RMG_consolidated.csv.

consolidate_chnosz.R: reads OBIGT CSVs, computes via subcrt(), writes CHNOSZ_consolidated.csv.

make_rda.R: converts the 3 consolidated CSVs into .rda files in data/.

Then update inst/extdata/OBIGT/SOURCE.txt with the CHNOSZ version used.

License and credits
GPL-3 (as CHNOSZ, whose data are reused). If you use these data, please cite:

ATcT: Ruscic, B. (2024). Active Thermochemical Tables, version 1.228.

RMG: Gao, C. W., et al. (2016). RMG-database.

CHNOSZ: Dick, J. M. (2019). CHNOSZ: Thermodynamic calculations and diagrams for geochemistry. Frontiers in Earth Science, 7, 152. doi:10.3389/feart.2019.00180.

References for each value are in the ref1 and ref2 columns (see get_references()).

