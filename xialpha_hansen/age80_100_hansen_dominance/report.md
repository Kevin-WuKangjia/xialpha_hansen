# age80_100 xialpha Hansen-style dominance analysis

- Figure: `E:\b2f10k\xialpha_hansen\age80_100_hansen_dominance\age80_100_hansen_dominance.png`
- SVG: `E:\b2f10k\xialpha_hansen\age80_100_hansen_dominance\age80_100_hansen_dominance.svg`
- Spins: `10000`
- Significance stars: `spin_p < 0.05`
- Train split: `0.75`
- Random seed: `1234`

## Metadata
- analysis_scope: `age_bin`
- age_bin: `age80_100`
- subject_count: `19`
- raw_maps_root: `E:\b2f10k\averagemaps\age80_100`
- dataset_key: `age80_100`

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
- XA: adjusted R^2=0.5720, null mean=0.4846, null p95=0.7319, spin p=0.2659, FDR q=0.3102, train corr mean=0.8087, test corr mean=0.4244
- XB: adjusted R^2=0.5628, null mean=0.4883, null p95=0.6407, spin p=0.3151, FDR q=0.3151, train corr mean=0.8209, test corr mean=0.4139
- XE: adjusted R^2=0.6497, null mean=0.4750, null p95=0.7180, spin p=0.1991, FDR q=0.2787, train corr mean=0.8641, test corr mean=0.2347
- AA: adjusted R^2=0.2900, null mean=0.0348, null p95=0.2638, spin p=0.0359, FDR q=0.06282, train corr mean=0.6718, test corr mean=0.4435
- AB: adjusted R^2=0.4011, null mean=0.0771, null p95=0.2351, spin p=0.0005999, FDR q=0.0014, train corr mean=0.7486, test corr mean=0.3975
- AE: adjusted R^2=0.4280, null mean=0.0835, null p95=0.2274, spin p=9.999e-05, FDR q=0.00035, train corr mean=0.7609, test corr mean=0.4171
- AF: adjusted R^2=0.4167, null mean=0.0776, null p95=0.2251, spin p=9.999e-05, FDR q=0.00035, train corr mean=0.7578, test corr mean=0.3951