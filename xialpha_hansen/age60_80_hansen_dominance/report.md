# age60_80 xialpha Hansen-style dominance analysis

- Figure: `E:\b2f10k\xialpha_hansen\age60_80_hansen_dominance\age60_80_hansen_dominance.png`
- SVG: `E:\b2f10k\xialpha_hansen\age60_80_hansen_dominance\age60_80_hansen_dominance.svg`
- Spins: `10000`
- Significance stars: `spin_p < 0.05`
- Train split: `0.75`
- Random seed: `1234`

## Metadata
- analysis_scope: `age_bin`
- age_bin: `age60_80`
- subject_count: `253`
- raw_maps_root: `E:\b2f10k\averagemaps\age60_80`
- dataset_key: `age60_80`

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
- XA: adjusted R^2=0.5233, null mean=0.4294, null p95=0.6692, spin p=0.2488, FDR q=0.2488, train corr mean=0.7995, test corr mean=0.5110
- XB: adjusted R^2=0.7271, null mean=0.5464, null p95=0.7661, spin p=0.1531, FDR q=0.1786, train corr mean=0.8901, test corr mean=0.6284
- XE: adjusted R^2=0.7026, null mean=0.5107, null p95=0.7351, spin p=0.1274, FDR q=0.1783, train corr mean=0.8876, test corr mean=0.5173
- AA: adjusted R^2=0.3957, null mean=0.0160, null p95=0.1738, spin p=0.0025, FDR q=0.004375, train corr mean=0.7287, test corr mean=0.3741
- AB: adjusted R^2=0.5398, null mean=0.1060, null p95=0.2736, spin p=0.0002, FDR q=0.0004666, train corr mean=0.8148, test corr mean=0.3993
- AE: adjusted R^2=0.5546, null mean=0.1298, null p95=0.3036, spin p=0.0002, FDR q=0.0004666, train corr mean=0.8191, test corr mean=0.4342
- AF: adjusted R^2=0.5488, null mean=0.1191, null p95=0.2895, spin p=0.0002, FDR q=0.0004666, train corr mean=0.8175, test corr mean=0.4166