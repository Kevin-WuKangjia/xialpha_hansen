# -*- coding: utf-8 -*-
"""
Figure 4: Dynamics
"""

import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations
from scipy.stats import zscore, pearsonr, ttest_ind
import seaborn as sns
from matplotlib.colors import ListedColormap
import matplotlib.gridspec as gridspec
from scipy.spatial.distance import squareform, pdist
from scipy.spatial import cKDTree
from sklearn.linear_model import LinearRegression
from statsmodels.stats.multitest import multipletests
from neuromaps.nulls.spins import gen_spinsamples
from neuromaps.datasets import fetch_atlas
from nilearn import plotting
import nibabel as nib
import os

def get_reg_r_sq(X, y):
    lin_reg = LinearRegression()
    lin_reg.fit(X, y)
    yhat = lin_reg.predict(X)
    SS_Residual = sum((y - yhat) ** 2)
    SS_Total = sum((y - np.mean(y)) ** 2)
    r_squared = 1 - (float(SS_Residual)) / SS_Total
    adjusted_r_squared = 1 - (1 - r_squared) * \
        (len(y) - 1) / (len(y) - X.shape[1] - 1)
    return adjusted_r_squared


def get_dominance_stats(X, y):
    """
    Minimal replacement for netneurotools.stats.get_dominance_stats.
    """
    def remove_ret(tpl, elem):
        lst = list(tpl)
        lst.remove(elem)
        return tuple(lst)

    n_predictor = X.shape[-1]
    predictor_combs = [list(combinations(range(n_predictor), i))
                       for i in range(1, n_predictor + 1)]

    model_r_sq = {}
    for len_group in predictor_combs:
        for idx_tuple in len_group:
            model_r_sq[idx_tuple] = get_reg_r_sq(X[:, idx_tuple], y)

    model_metrics = {}
    individual_dominance = np.array([model_r_sq[(i_pred,)]
                                     for i_pred in range(n_predictor)]).reshape(1, -1)
    model_metrics["individual_dominance"] = individual_dominance

    partial_dominance = [[] for _ in range(n_predictor - 1)]
    for i_len in range(n_predictor - 1):
        i_len_combs = list(combinations(range(n_predictor), i_len + 2))
        for j_node in range(n_predictor):
            j_node_sel = [v for v in i_len_combs if j_node in v]
            reduced_list = [remove_ret(comb, j_node) for comb in j_node_sel]
            diff_values = [model_r_sq[j_node_sel[i]] - model_r_sq[reduced_list[i]]
                           for i in range(len(reduced_list))]
            partial_dominance[i_len].append(float(np.mean(diff_values)))

    partial_dominance = np.array(partial_dominance)
    model_metrics["partial_dominance"] = partial_dominance

    total_dominance = np.mean(np.r_[individual_dominance, partial_dominance], axis=0)
    model_metrics["total_dominance"] = total_dominance
    model_metrics["full_r_sq"] = model_r_sq[tuple(range(n_predictor))]
    return model_metrics, model_r_sq


def cv_slr_distance_dependent(X, y, coords, train_pct=.75, metric='rsq'):
    '''
    cross validates linear regression model using distance-dependent method.
    X = n x p matrix of input variables
    y = n x 1 matrix of output variable
    coords = n x 3 coordinates of each observation
    train_pct (between 0 and 1), percent of observations in training set
    metric = {'rsq', 'corr'}
    '''

    P = squareform(pdist(coords, metric="euclidean"))
    train_metric = []
    test_metric = []

    for i in range(len(y)):
        distances = P[i, :]  # for every node
        idx = np.argsort(distances)

        train_idx = idx[:int(np.floor(train_pct * len(coords)))]
        test_idx = idx[int(np.floor(train_pct * len(coords))):]

        mdl = LinearRegression()
        mdl.fit(X[train_idx, :], y[train_idx])
        if metric == 'rsq':
            # get r^2 of train set
            train_metric.append(get_reg_r_sq(X[train_idx, :], y[train_idx]))

        elif metric == 'corr':
            rho, _ = pearsonr(mdl.predict(X[train_idx, :]), y[train_idx])
            train_metric.append(rho)

        yhat = mdl.predict(X[test_idx, :])
        if metric == 'rsq':
            # get r^2 of test set
            SS_Residual = sum((y[test_idx] - yhat) ** 2)
            SS_Total = sum((y[test_idx] - np.mean(y[test_idx])) ** 2)
            r_squared = 1 - (float(SS_Residual)) / SS_Total
            adjusted_r_squared = 1-(1-r_squared)*((len(y[test_idx]) - 1) /
                                                  (len(y[test_idx]) -
                                                   X.shape[1]-1))
            test_metric.append(adjusted_r_squared)

        elif metric == 'corr':
            rho, _ = pearsonr(yhat, y[test_idx])
            test_metric.append(rho)

    return train_metric, test_metric


def get_reg_r_pval(X, y, spins, nspins):
    emp = get_reg_r_sq(X, y)
    null = np.zeros((nspins, ))
    for s in range(nspins):
        null[s] = get_reg_r_sq(X[spins[:, s], :], y)
    return (1 + sum(null > emp))/(nspins + 1)


def get_perm_p(emp, null):
    return (1 + sum(abs(null - np.mean(null))
                    > abs(emp - np.mean(null)))) / (len(null) + 1)


def find_schaefer_annot_paths():
    candidates = [
        'C:/Users/Administrator/nnt-data/atl-schaefer2018/fsaverage',
        'C:/Users/Administrator/nilearn_data/schafer_surface',
    ]
    fname_l = 'atl-Schaefer2018_space-fsaverage_hemi-L_desc-100Parcels7Networks_deterministic.annot'
    fname_r = 'atl-Schaefer2018_space-fsaverage_hemi-R_desc-100Parcels7Networks_deterministic.annot'
    for root in candidates:
        l_path = os.path.join(root, fname_l)
        r_path = os.path.join(root, fname_r)
        if os.path.exists(l_path) and os.path.exists(r_path):
            return l_path, r_path
    raise FileNotFoundError(f'Cannot find Schaefer annot files. Tried: {candidates}')


def load_surf_coords(path_):
    return np.asarray(nib.load(path_).darrays[0].data, dtype=float)


def build_fsaverage10k_schaefer100_labels():
    annot_l, annot_r = find_schaefer_annot_paths()
    labels_l_164k, _, _ = nib.freesurfer.read_annot(annot_l)
    labels_r_164k, _, _ = nib.freesurfer.read_annot(annot_r)

    fsavg_164k = fetch_atlas('fsaverage', density='164k')
    fsavg_10k = fetch_atlas('fsaverage', density='10k')

    sph_l_164k = load_surf_coords(fsavg_164k['sphere'].L)
    sph_r_164k = load_surf_coords(fsavg_164k['sphere'].R)
    sph_l_10k = load_surf_coords(fsavg_10k['sphere'].L)
    sph_r_10k = load_surf_coords(fsavg_10k['sphere'].R)

    tree_l = cKDTree(sph_l_164k)
    tree_r = cKDTree(sph_r_164k)
    _, idx_l = tree_l.query(sph_l_10k, k=1)
    _, idx_r = tree_r.query(sph_r_10k, k=1)

    labels_l_10k = labels_l_164k[idx_l].astype(np.int32)
    labels_r_10k = labels_r_164k[idx_r].astype(np.int32)
    return labels_l_10k, labels_r_10k


def parcel100_to_fsaverage10k_surface(parcel_vec, labels_l_10k, labels_r_10k):
    left = np.zeros(labels_l_10k.shape[0], dtype=float)
    right = np.zeros(labels_r_10k.shape[0], dtype=float)
    for pid in range(1, 51):
        left[labels_l_10k == pid] = parcel_vec[pid - 1]
        right[labels_r_10k == pid] = parcel_vec[50 + pid - 1]
    return left, right


def plot_combined_panel(power, power_band, dominance_plot, receptor_names, cmap_seq):
    atlas = fetch_atlas("fsaverage", density="10k")
    surf_l = atlas["inflated"].L
    surf_r = atlas["inflated"].R
    labels_l_10k, labels_r_10k = build_fsaverage10k_schaefer100_labels()

    nfeat = len(power_band)
    bar_vals = np.sum(dominance_plot, axis=1)
    heat_data = np.zeros_like(dominance_plot)
    np.divide(
        dominance_plot,
        np.sum(dominance_plot, axis=1, keepdims=True),
        out=heat_data,
        where=np.sum(dominance_plot, axis=1, keepdims=True) > 0,
    )

    fig = plt.figure(figsize=(18, 2.4 * nfeat))
    gs = gridspec.GridSpec(
        nrows=nfeat,
        ncols=4,
        width_ratios=[1.2, 1.2, 1.0, 2.0],
        wspace=0.15,
        hspace=0.05,
    )

    for i, feat in enumerate(power_band):
        left_vals, right_vals = parcel100_to_fsaverage10k_surface(power[:, i], labels_l_10k, labels_r_10k)
        both = np.r_[left_vals, right_vals]
        vmin = np.nanpercentile(both, 2)
        vmax = np.nanpercentile(both, 98)
        if np.isclose(vmin, vmax):
            vmax = vmin + 1e-6

        ax_l = fig.add_subplot(gs[i, 0], projection='3d')
        plotting.plot_surf_stat_map(
            surf_l, left_vals, hemi='left', view='lateral',
            axes=ax_l, colorbar=False, bg_on_data=False,
            threshold=None, vmin=vmin, vmax=vmax,
        )
        ax_l.set_axis_off()
        ax_l.text2D(-0.08, 0.5, feat, transform=ax_l.transAxes, va='center', ha='right', fontsize=9)
        if i == 0:
            ax_l.set_title("Left hemisphere", fontsize=10, pad=4)

        ax_r = fig.add_subplot(gs[i, 1], projection='3d')
        plotting.plot_surf_stat_map(
            surf_r, right_vals, hemi='right', view='lateral',
            axes=ax_r, colorbar=False, bg_on_data=False,
            threshold=None, vmin=vmin, vmax=vmax,
        )
        ax_r.set_axis_off()
        if i == 0:
            ax_r.set_title("Right hemisphere", fontsize=10, pad=4)

    ax_bar = fig.add_subplot(gs[:, 2])
    y = np.arange(nfeat)
    ax_bar.barh(y, bar_vals, color='#4C72B0')
    ax_bar.set_yticks(y)
    ax_bar.set_yticklabels(power_band, fontsize=9)
    ax_bar.invert_yaxis()
    ax_bar.set_xlabel('R^2', fontsize=10)
    ax_bar.set_title('Model fit by feature', fontsize=11)
    ax_bar.grid(axis='x', alpha=0.25)

    ax_heat = fig.add_subplot(gs[:, 3])
    sns.heatmap(
        heat_data,
        xticklabels=receptor_names,
        yticklabels=power_band,
        cmap=cmap_seq,
        linewidths=.3,
        ax=ax_heat,
        cbar_kws={"label": "Normalized dominance"},
    )
    ax_heat.set_xlabel('Receptor data', fontsize=10)
    ax_heat.set_ylabel('')
    ax_heat.set_title('Receptor dominance', fontsize=11)
    ax_heat.tick_params(axis='x', labelrotation=90, labelsize=8)
    ax_heat.tick_params(axis='y', labelsize=9)

    fig.tight_layout()
    fig.savefig(pathresult + 'figures/dominance_panel_3part.png', dpi=300, bbox_inches='tight')
    plt.close(fig)

"""
set-up
"""

path = 'E:/b2f10k/hansen_receptors-main/'
pathresult = 'E:/b2f10k/xialpha_hansen/age0_20/'
scale = 'scale100'
apply_fdr_zero = False
os.makedirs(pathresult + 'results', exist_ok=True)
os.makedirs(pathresult + 'figures', exist_ok=True)

coords = np.genfromtxt(path+'data/schaefer/coordinates/Schaefer_100_centres.txt')[:, 1:]
nnodes = coords.shape[0]
hemiid = np.zeros((nnodes, ))
hemiid[:int(nnodes/2)] = 1
nspins = 100
spins = gen_spinsamples(coords, hemiid, n_rotate=nspins, seed=1234)
if isinstance(spins, tuple):
    spins = spins[0]

# load MEG power
power = np.genfromtxt('E:/b2f10k/xialpha_parcellate/age0_20/xialpha_parcellated_schaefer100.csv', delimiter=',', skip_header=1)
power_band = ["Alpha_estimate_Exponent",
  "Alpha_PAF",
  "Alpha_Power",
  "Alpha_Width",
  "Xi_Exponent",
  "Xi_Power",
  "Xi_Width"]
if power.shape[1] != len(power_band):
    raise ValueError(f"power columns ({power.shape[1]}) != power_band length ({len(power_band)})")

# load the receptor data
receptor_data = np.genfromtxt(path+'results/receptor_data_'+scale+'.csv', delimiter=',')
receptor_names = np.load(path+'data/receptor_names_pet.npy')

# colourmaps
cmap = np.genfromtxt(path+'data/colourmap.csv', delimiter=',')
cmap_div = ListedColormap(cmap)
cmap_seq = ListedColormap(cmap[128:, :])


"""
Dominance analysis
"""

model_metrics = dict([])
train_metric = np.zeros([nnodes, len(power_band)])
test_metric = np.zeros(train_metric.shape)
model_pval = np.zeros((len(power_band), ))

for i in range(len(power_band)):
    print(i)
    m, _ = get_dominance_stats(zscore(receptor_data),
                               zscore(power[:, i]))
    model_metrics[power_band[i]] = m
    # cross validate the model
    train_metric[:, i], test_metric[:, i] = \
        cv_slr_distance_dependent(zscore(receptor_data),
                                  zscore(power[:, i]),
                                  coords, .75,
                                  metric='corr')
    # get model pval
    model_pval[i] = get_reg_r_pval(zscore(receptor_data),
                                   zscore(power[:, i]), 
                                   spins, nspins)

dominance = np.zeros((len(power_band), len(receptor_names)))

for i in range(len(model_metrics)):
    tmp = model_metrics[power_band[i]]
    dominance[i, :] = tmp["total_dominance"]
dominance_raw = dominance.copy()
np.save(pathresult+'results/dominance_power.npy', dominance_raw)
np.save(pathresult+'results/power_cv_train.npy', train_metric)
np.save(pathresult+'results/power_cv_test.npy', test_metric)

model_pval_raw = model_pval.copy()
model_pval = multipletests(model_pval, method='fdr_bh')[1]
if apply_fdr_zero:
    dominance[np.where(model_pval >= 0.05)[0], :] = 0
dominance_sig = dominance.copy()
np.save(pathresult+'results/model_pval_raw.npy', model_pval_raw)
np.save(pathresult+'results/model_pval_fdr.npy', model_pval)
np.save(pathresult+'results/dominance_power_sig.npy', dominance_sig)
np.save(pathresult+'results/dominance_power_plot.npy', dominance_sig if apply_fdr_zero else dominance_raw)
print("model_pval_raw:", model_pval_raw)
print("model_pval_fdr:", model_pval)
print("apply_fdr_zero:", apply_fdr_zero)
print("dominance_raw_row_sum:", np.sum(dominance_raw, axis=1))
print("dominance_sig_row_sum:", np.sum(dominance_sig, axis=1))
dominance_plot = dominance_sig if apply_fdr_zero else dominance_raw
plot_combined_panel(power, power_band, dominance_plot, receptor_names, cmap_seq)

# plot cross validation
fig, (ax1, ax2) = plt.subplots(2, 1)
sns.violinplot(data=train_metric, ax=ax1)
sns.violinplot(data=test_metric, ax=ax2)
ax1.set(ylabel='train set correlation', ylim=(-1, 1))
ax2.set_xticklabels(power_band, rotation=90)
ax2.set(ylabel='test set correlation', ylim=(-1, 1))
plt.tight_layout()
plt.savefig(pathresult+'figures/violin_crossval_power.png', dpi=300)

# compare dominance across receptor classes
exc = ['5HT2a', '5HT4', '5HT6', 'D1', 'mGluR5', 'A4B2', 'M1', 'NMDA']
inh = ['5HT1a', '5HT1b', 'CB1', 'D2', 'GABAa', 'H3', 'MOR']
mami = ['5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', '5HTT', 'D1',
        'D2', 'DAT', 'H3', 'NET']
nmami = list(set(receptor_names) - set(mami))
metab = ['5HT1a', '5HT1b', '5HT2a', '5HT4', '5HT6', 'CB1', 'D1',
         'D2', 'H3', 'M1', 'mGluR5', 'MOR']
iono = ['A4B2', 'GABAa', 'NMDA']
gspath = ['5HT4', '5HT6', 'D1']
gipath = ['CB1', 'D2', 'H3', '5HT1a', '5HT1b', 'MOR']
gqpath = ['5HT2a', 'mGluR5', 'M1']

i_exc = np.array([list(receptor_names).index(i) for i in exc])
i_inh = np.array([list(receptor_names).index(i) for i in inh])
i_mami = np.array([list(receptor_names).index(i) for i in mami])
i_nmami = np.array([list(receptor_names).index(i) for i in nmami])
i_metab =  np.array([list(receptor_names).index(i) for i in metab])
i_iono = np.array([list(receptor_names).index(i) for i in iono])
i_gs = np.array([list(receptor_names).index(i) for i in gspath])
i_gi = np.array([list(receptor_names).index(i) for i in gipath])
i_gq = np.array([list(receptor_names).index(i) for i in gqpath])

classes = [[i_exc, i_inh], [i_mami, i_nmami],
           [i_metab, i_iono], [i_gs, i_gi, i_gq]]
class_names = [['exc', 'inh'], ['monoamine', 'not'],
               ['metabotropic', 'ionotropic'], ['gs', 'gi', 'gq']]
plt.ion()
fig, axs = plt.subplots(1, 4, figsize=(15, 3))
axs = axs.ravel()
for i in range(len(classes)):
    print(class_names[i])
    d = [dominance[:, classes[i][j]].flatten() for j in range(len(classes[i]))]
    print(ttest_ind(d[0], d[1]))
    if len(d) > 2:
        print(ttest_ind(d[0], d[2]))
        print(ttest_ind(d[1], d[2]))
    sns.violinplot(data=d, inner=None, color=".8", ax=axs[i])
    sns.stripplot(data=d, ax=axs[i])
    axs[i].set_xticklabels(class_names[i])
    axs[i].set_ylabel('dominance (power)')
plt.tight_layout()
plt.savefig(pathresult+'figures/stripplot_power_rclasses.png', dpi=300)
