# Xi/Alpha fsaverage10k Pipeline

This repository converts Brainstorm/Matlab cortical maps to `fsaverage10k`, computes age-group/global averages, and runs Hansen-style receptor dominance analysis.

## Before Running

Set your own project root first.  
Do **not** use machine-specific paths from other users.

Recommended convention:

- `<PROJECT_ROOT>/result`
- `<PROJECT_ROOT>/averagemaps`
- `<PROJECT_ROOT>/xialpha_parcellate`
- `<PROJECT_ROOT>/xialpha_regression`
- `<PROJECT_ROOT>/xialpha_CCA`
- `<PROJECT_ROOT>/xialpha_hansen`
- `<PROJECT_ROOT>/figure`
- `<PROJECT_ROOT>/hansen_receptors-main`

In each script, update root constants before running:

- `matlabTopy.py`: input/output roots used for conversion.
- `averagemaps.py`: `RESULT_ROOT`, `META_JSON_PATH`, `OUTPUT_ROOT`.
- `xialpha_Schaefer100.py`: `RESULT_ROOT`, `META_JSON_PATH`, `OUTPUT_ROOT`.
- `xialpha_regression.py`: `INPUT_ROOT`, `OUTPUT_ROOT`.
- `xialpha_CCA.py`: `REGRESSION_ROOT`, `REFERENCE_DIR`, `OUTPUT_ROOT`, `FIGURE_ROOT`.
- `run_xialpha_hansen_dominance.py`: `REPO_ROOT`, `AVERAGE_MAPS_ROOT`, `OUTPUTS_ROOT`, `REFERENCE_DIR`, `HANSEN_ROOT`, `CACHE_ROOT`.

## Dependencies

- Python 3.10+
- Main packages:
  - `numpy`
  - `scipy`
  - `nibabel`
  - `nilearn`
  - `neuromaps`
  - `matplotlib`
  - `seaborn`
  - `statsmodels`
  - `scikit-learn`
- Required toolkit:
  - `mat2py_maps`: https://github.com/ronald1129/mat2py_maps
- Hansen reference/code:
  - https://github.com/netneurolab/hansen_receptors/tree/main/code

## Pipeline Steps

### 1) Convert Brainstorm fsaverage8k to fsaverage10k `.gii`

Run:

```powershell
python <PROJECT_ROOT>/matlabTopy.py
```

Output location (typical):

- `<PROJECT_ROOT>/result/<SubjectID>/<Alpha_or_Xi>/<Feature>/maps/*.shape.gii`
- (optional, depending on your conversion settings): `figures/`, `reports/`

### 2) Build age-group and global average maps

Run:

```powershell
python <PROJECT_ROOT>/averagemaps.py
```

Output location:

- `<PROJECT_ROOT>/averagemaps/age0_20/.../maps/*.shape.gii`
- `<PROJECT_ROOT>/averagemaps/age20_40/.../maps/*.shape.gii`
- `<PROJECT_ROOT>/averagemaps/age40_60/.../maps/*.shape.gii`
- `<PROJECT_ROOT>/averagemaps/age60_80/.../maps/*.shape.gii`
- `<PROJECT_ROOT>/averagemaps/age80_100/.../maps/*.shape.gii`
- `<PROJECT_ROOT>/averagemaps/ageall/.../maps/*.shape.gii`
- `<PROJECT_ROOT>/averagemaps/age_groups_summary.json`
- `<PROJECT_ROOT>/averagemaps/averagemaps_summary.json`

### 3) Subject-level Schaefer100 parcellation for Xi/Alpha (age-ready table)

Run:

```powershell
python <PROJECT_ROOT>/xialpha_Schaefer100.py
```

Output location:

- `<PROJECT_ROOT>/xialpha_parcellate/xi/Power/Xi_estimate_Power_schaefer100.csv`
- `<PROJECT_ROOT>/xialpha_parcellate/xi/Width/Xi_estimate_Width_schaefer100.csv`
- `<PROJECT_ROOT>/xialpha_parcellate/xi/Exponent/Xi_estimate_Exponent_schaefer100.csv`
- `<PROJECT_ROOT>/xialpha_parcellate/alpha/Power/Alpha_estimate_Power_schaefer100.csv`
- `<PROJECT_ROOT>/xialpha_parcellate/alpha/Width/Alpha_estimate_Width_schaefer100.csv`
- `<PROJECT_ROOT>/xialpha_parcellate/alpha/Exponent/Alpha_estimate_Exponent_schaefer100.csv`
- `<PROJECT_ROOT>/xialpha_parcellate/alpha/PAF/Alpha_estimate_PAF_schaefer100.csv`
- `<PROJECT_ROOT>/xialpha_parcellate/run_summary.json`

Each CSV format:

- Column 1: `subject_name`
- Column 2: `age`
- Columns 3..102: `region1..region100`

### 4) Age regression (linear + nonlinear) to get age-effect maps

Run:

```powershell
python <PROJECT_ROOT>/xialpha_regression.py
```

Output location:

- Linear outputs: `<PROJECT_ROOT>/xialpha_regression/linear/<xi_or_alpha>/<feature>/`
- Nonlinear outputs: `<PROJECT_ROOT>/xialpha_regression/nonlinear/<method>/<xi_or_alpha>/<feature>/`
- Global summary: `<PROJECT_ROOT>/xialpha_regression/run_summary.json`

Main age-effect files:

- Linear: `beta_age_map.csv`
- Nonlinear (`quadratic`, `cubic`, `spline_df4`): `delta_pred_map.csv`

### 5) CCA between Xi/Alpha age-effects and receptor maps

Run:

```powershell
python <PROJECT_ROOT>/xialpha_CCA.py
```

Output location:

- Numeric outputs:
  - `<PROJECT_ROOT>/xialpha_CCA/linear/GLM/`
  - `<PROJECT_ROOT>/xialpha_CCA/nonlinear/quadratic/`
  - `<PROJECT_ROOT>/xialpha_CCA/nonlinear/cubic/`
  - `<PROJECT_ROOT>/xialpha_CCA/nonlinear/spline_df4/`
- CCA figures:
  - `<PROJECT_ROOT>/figure/CCA/linear/GLM/cca_overview.png`
  - `<PROJECT_ROOT>/figure/CCA/nonlinear/<method>/cca_overview.png`

### 6) Run Hansen-style dominance analysis (optional complementary analysis)

Single age bin:

```powershell
python <PROJECT_ROOT>/run_xialpha_hansen_dominance.py --age-bin age0_20 --n-spins 10000
```

Global weighted:

```powershell
python <PROJECT_ROOT>/run_xialpha_hansen_dominance.py --global-weighted --n-spins 10000
```

Output location:

- `<PROJECT_ROOT>/xialpha_hansen/<age>_hansen_dominance/`
- `<PROJECT_ROOT>/xialpha_hansen/global_weighted_hansen_dominance/`

## Function-Level Output Map

### `averagemaps.py`

- `main()`:
  - writes average maps (`*.shape.gii`) under `<PROJECT_ROOT>/averagemaps/<age_group>/<source>/<feature>/maps/`
  - writes per-feature summary: `.../<feature>/summary.json`
  - writes global summaries:
    - `<PROJECT_ROOT>/averagemaps/age_groups_summary.json`
    - `<PROJECT_ROOT>/averagemaps/averagemaps_summary.json`

### `run_xialpha_hansen_dominance.py`

- `ensure_reference_dir()`:
  - writes/copies Hansen reference files into `<PROJECT_ROOT>/hansen_reference/`
- `main()`:
  - writes full analysis outputs into:
    - `<PROJECT_ROOT>/xialpha_hansen/<age>_hansen_dominance/` or
    - `<PROJECT_ROOT>/xialpha_hansen/global_weighted_hansen_dominance/`
  - outputs include figure(s) and analysis artifacts generated by `run_hansen_dominance_analysis(...)`

### `xialpha_Schaefer100.py`

- `main()`:
  - reads all subject Xi/Alpha fsaverage10k maps under `<PROJECT_ROOT>/result/`
  - applies Xi/Alpha preprocessing rules (`signed_log1p`, alpha anchor weighting, xi positive mask)
  - parcellates to Schaefer100
  - writes 7 feature CSVs under `<PROJECT_ROOT>/xialpha_parcellate/<xi_or_alpha>/<feature>/`
  - writes:
    - per-feature `summary.json`
    - `<PROJECT_ROOT>/xialpha_parcellate/run_summary.json`
    - `<PROJECT_ROOT>/xialpha_parcellate/errors.json` (if any)

### `xialpha_regression.py`

- `main()`:
  - loads each feature CSV from `<PROJECT_ROOT>/xialpha_parcellate/`
  - fits:
    - linear GLM/OLS
    - nonlinear models: `quadratic`, `cubic`, `spline_df4`
  - writes linear results under `<PROJECT_ROOT>/xialpha_regression/linear/...`
  - writes nonlinear results under `<PROJECT_ROOT>/xialpha_regression/nonlinear/<method>/...`
  - writes `<PROJECT_ROOT>/xialpha_regression/run_summary.json`

### `xialpha_CCA.py`

- `main()`:
  - builds Xi/Alpha age-effect matrix (`100 x 7`) from regression outputs
  - loads receptor matrix (`100 x 19`) from `<PROJECT_ROOT>/hansen_reference/`
  - runs CCA for:
    - `linear/GLM`
    - `nonlinear/quadratic`
    - `nonlinear/cubic`
    - `nonlinear/spline_df4`
  - writes numeric CCA outputs under `<PROJECT_ROOT>/xialpha_CCA/<linear_or_nonlinear>/<method>/`
  - writes figures under `<PROJECT_ROOT>/figure/CCA/<linear_or_nonlinear>/<method>/`
  - writes `<PROJECT_ROOT>/xialpha_CCA/run_summary.json`

## Batch Run All Age Groups

```powershell
$SCRIPT = "<PROJECT_ROOT>/run_xialpha_hansen_dominance.py"

Get-ChildItem "<PROJECT_ROOT>/averagemaps" -Directory |
Where-Object { $_.Name -like "age*" -and $_.Name -ne "ageall" } |
ForEach-Object {
    python $SCRIPT --age-bin $_.Name --n-spins 10000
}
```

## Notes

- First run may be slower because `neuromaps` downloads atlas/transform caches.
- Use smaller `--n-spins` (for example `100`) for quick smoke tests.
