# age20_40 xialpha Hansen-style dominance analysis

- Figure: `E:\b2f10k\xialpha_hansen\age20_40_hansen_dominance\age20_40_hansen_dominance.png`
- SVG: `E:\b2f10k\xialpha_hansen\age20_40_hansen_dominance\age20_40_hansen_dominance.svg`
- Spins: `10000`
- Significance stars: `spin_p < 0.05`
- Train split: `0.75`
- Random seed: `1234`

## Metadata
- analysis_scope: `age_bin`
- age_bin: `age20_40`
- subject_count: `736`
- raw_maps_root: `E:\b2f10k\averagemaps\age20_40`
- dataset_key: `age20_40`

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
- XA: adjusted R^2=0.6120, null mean=0.4122, null p95=0.5708, spin p=0.0157, FDR q=0.0157, train corr mean=0.8499, test corr mean=0.3718
- XB: adjusted R^2=0.8416, null mean=0.5535, null p95=0.7050, spin p=9.999e-05, FDR q=0.0006999, train corr mean=0.9388, test corr mean=0.7091
- XE: adjusted R^2=0.7555, null mean=0.4247, null p95=0.6470, spin p=0.0005, FDR q=0.001575, train corr mean=0.9152, test corr mean=0.5368
- AA: adjusted R^2=0.4193, null mean=0.0311, null p95=0.1873, spin p=0.0027, FDR q=0.00315, train corr mean=0.7373, test corr mean=0.4611
- AB: adjusted R^2=0.5302, null mean=0.1405, null p95=0.3054, spin p=0.0012, FDR q=0.00168, train corr mean=0.8065, test corr mean=0.4824
- AE: adjusted R^2=0.5346, null mean=0.1569, null p95=0.3309, spin p=0.0008999, FDR q=0.001575, train corr mean=0.8047, test corr mean=0.4942
- AF: adjusted R^2=0.5317, null mean=0.1488, null p95=0.3155, spin p=0.0008999, FDR q=0.001575, train corr mean=0.8058, test corr mean=0.4868