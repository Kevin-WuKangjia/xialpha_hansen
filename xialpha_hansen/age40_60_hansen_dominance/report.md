# age40_60 xialpha Hansen-style dominance analysis

- Figure: `E:\b2f10k\xialpha_hansen\age40_60_hansen_dominance\age40_60_hansen_dominance.png`
- SVG: `E:\b2f10k\xialpha_hansen\age40_60_hansen_dominance\age40_60_hansen_dominance.svg`
- Spins: `10000`
- Significance stars: `spin_p < 0.05`
- Train split: `0.75`
- Random seed: `1234`

## Metadata
- analysis_scope: `age_bin`
- age_bin: `age40_60`
- subject_count: `120`
- raw_maps_root: `E:\b2f10k\averagemaps\age40_60`
- dataset_key: `age40_60`

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
- XA: adjusted R^2=0.6680, null mean=0.2753, null p95=0.4163, spin p=9.999e-05, FDR q=0.00035, train corr mean=0.8736, test corr mean=0.5053
- XB: adjusted R^2=0.7407, null mean=0.4812, null p95=0.6918, spin p=0.009799, FDR q=0.009799, train corr mean=0.8952, test corr mean=0.7659
- XE: adjusted R^2=0.8058, null mean=0.5096, null p95=0.7169, spin p=0.0014, FDR q=0.001633, train corr mean=0.9210, test corr mean=0.7802
- AA: adjusted R^2=0.4408, null mean=0.0492, null p95=0.2385, spin p=0.0009999, FDR q=0.0014, train corr mean=0.7474, test corr mean=0.4500
- AB: adjusted R^2=0.6220, null mean=0.1392, null p95=0.3227, spin p=0.0002, FDR q=0.00035, train corr mean=0.8513, test corr mean=0.5462
- AE: adjusted R^2=0.6319, null mean=0.1550, null p95=0.3411, spin p=0.0002, FDR q=0.00035, train corr mean=0.8553, test corr mean=0.5590
- AF: adjusted R^2=0.6257, null mean=0.1485, null p95=0.3333, spin p=0.0002, FDR q=0.00035, train corr mean=0.8525, test corr mean=0.5574