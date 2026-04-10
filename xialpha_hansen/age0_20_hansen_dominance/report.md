# age0_20 xialpha Hansen-style dominance analysis

- Figure: `E:\b2f10k\xialpha_hansen\age0_20_hansen_dominance\age0_20_hansen_dominance.png`
- SVG: `E:\b2f10k\xialpha_hansen\age0_20_hansen_dominance\age0_20_hansen_dominance.svg`
- Spins: `10000`
- Significance stars: `spin_p < 0.05`
- Train split: `0.75`
- Random seed: `1234`

## Metadata
- analysis_scope: `age_bin`
- age_bin: `age0_20`
- subject_count: `837`
- raw_maps_root: `E:\b2f10k\averagemaps\age0_20`
- dataset_key: `age0_20`

## Notes
- A signed log1p transform was applied to all xi/alpha maps before weighting, plotting, parcellation, and regression.
- Alpha-family maps (AA, AB, AE, AF) were weighted by thresholded normalized raw AA anchors at 0.05.
- Xi-family maps (XA, XB, XE) were masked only by raw XA > 0 without normalization.
- All predictor and target columns are max-scaled before z-scoring and regression.
- Off-cortex vertices are excluded from both the figure and the parcel-wise analysis.
- The Schaefer100 atlas is projected to native fsaverage10k using neuromaps and cached locally.
- Feature display order is fixed to: XA, XB, XE, AA, AB, AE, AF.

## Surface NaN summary
- XA: off-cortex NaNs excluded from analysis/plotting, L=888, R=881
- XB: off-cortex NaNs excluded from analysis/plotting, L=888, R=881
- XE: off-cortex NaNs excluded from analysis/plotting, L=888, R=881
- AA: off-cortex NaNs excluded from analysis/plotting, L=888, R=881
- AB: off-cortex NaNs excluded from analysis/plotting, L=888, R=881
- AE: off-cortex NaNs excluded from analysis/plotting, L=888, R=881
- AF: off-cortex NaNs excluded from analysis/plotting, L=888, R=881

## Model summary
- XA: adjusted R^2=0.6496, null mean=0.5442, null p95=0.7283, spin p=0.1914, FDR q=0.2233, train corr mean=0.8495, test corr mean=0.3998
- XB: adjusted R^2=0.8456, null mean=0.7571, null p95=0.8950, spin p=0.2535, FDR q=0.2535, train corr mean=0.9324, test corr mean=0.6156
- XE: adjusted R^2=0.7470, null mean=0.4949, null p95=0.7334, spin p=0.0301, FDR q=0.04214, train corr mean=0.9162, test corr mean=0.4783
- AA: adjusted R^2=0.4438, null mean=0.0236, null p95=0.1819, spin p=0.0021, FDR q=0.003675, train corr mean=0.7374, test corr mean=0.4967
- AB: adjusted R^2=0.5719, null mean=0.1340, null p95=0.3025, spin p=0.0003, FDR q=0.0006999, train corr mean=0.8191, test corr mean=0.5522
- AE: adjusted R^2=0.5854, null mean=0.1444, null p95=0.3178, spin p=0.0002, FDR q=0.0006999, train corr mean=0.8243, test corr mean=0.5814
- AF: adjusted R^2=0.5801, null mean=0.1415, null p95=0.3106, spin p=0.0002, FDR q=0.0006999, train corr mean=0.8222, test corr mean=0.5695