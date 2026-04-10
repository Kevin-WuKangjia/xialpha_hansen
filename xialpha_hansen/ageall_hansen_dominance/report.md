# ageall xialpha Hansen-style dominance analysis

- Figure: `E:\b2f10k\xialpha_hansen\ageall_hansen_dominance\ageall_hansen_dominance.png`
- SVG: `E:\b2f10k\xialpha_hansen\ageall_hansen_dominance\ageall_hansen_dominance.svg`
- Spins: `10000`
- Significance stars: `spin_p < 0.05`
- Train split: `0.75`
- Random seed: `1234`

## Metadata
- analysis_scope: `age_bin`
- age_bin: `ageall`
- subject_count: `0`
- raw_maps_root: `E:\b2f10k\averagemaps\ageall`
- dataset_key: `ageall`

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
- XA: adjusted R^2=0.6429, null mean=0.4514, null p95=0.6175, spin p=0.0263, FDR q=0.0263, train corr mean=0.8556, test corr mean=0.4189
- XB: adjusted R^2=0.8832, null mean=0.7273, null p95=0.8587, spin p=0.008299, FDR q=0.009682, train corr mean=0.9477, test corr mean=0.7240
- XE: adjusted R^2=0.7192, null mean=0.4408, null p95=0.6670, spin p=0.005799, FDR q=0.008119, train corr mean=0.8997, test corr mean=0.5188
- AA: adjusted R^2=0.4353, null mean=0.0292, null p95=0.1872, spin p=0.0023, FDR q=0.004025, train corr mean=0.7409, test corr mean=0.4859
- AB: adjusted R^2=0.5604, null mean=0.1402, null p95=0.3077, spin p=0.0005, FDR q=0.001167, train corr mean=0.8188, test corr mean=0.5230
- AE: adjusted R^2=0.5689, null mean=0.1537, null p95=0.3285, spin p=0.0004, FDR q=0.001167, train corr mean=0.8209, test corr mean=0.5495
- AF: adjusted R^2=0.5639, null mean=0.1484, null p95=0.3180, spin p=0.0005, FDR q=0.001167, train corr mean=0.8198, test corr mean=0.5354