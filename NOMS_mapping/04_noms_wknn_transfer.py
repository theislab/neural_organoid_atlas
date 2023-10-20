import warnings

warnings.filterwarnings("ignore")

import os
import sys

import numpy as np
import pandas as pd

import scvi
from scvi.model.utils import mde

import matplotlib.pyplot as plt

import scanpy as sc


os.chdir("/home/fleckj/projects/atlas/scripts")
import importlib
from interfaces import vae_interface
import wknn

wknn = importlib.reload(wknn)

#### Read stuff ####
REF_PATH = "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_common_hv2k.h5ad"

ATLAS_PATH = "/home/fleckj/projects/atlas/scratch/v3/230510_phase3_hv2k_braun.h5ad"

QUERY_PATH = "/home/fleckj/scratch/data/public_datasets/brain_organoids/AminPasca2023brx/230605_pasca_all_v1_common_hv2k.h5ad"

adata_braun = sc.read(REF_PATH)
print(adata_braun)

adata_pasca = sc.read(QUERY_PATH)
print(adata_pasca)


DATA_DIR = "/home/fleckj/projects/atlas/data/results/pasca_scarches/v4_artur/"

pasca_braun_latent = pd.read_csv(
    f"{DATA_DIR}/pasca_map_braun_latent.csv",
    index_col=0,
)

cell_names = pasca_braun_latent.index
cell_names_clean = cell_names.str.replace("-\\d$", "")

pasca_idx = cell_names_clean.isin(adata_pasca.obs_names)
pasca_adata_idx = adata_pasca.obs_names.isin(cell_names_clean)

adata_pasca = adata_pasca[pasca_adata_idx, :]

latent = pasca_braun_latent.values
latent_pasca = latent[pasca_idx, :]
latent_braun = latent[~pasca_idx, :]


#### Compute WKNN ####
wknn_scanvi_q2r, adjs_scanvi_q2r = wknn.get_wknn(
    ref=latent_braun,
    query=latent_pasca,
    k=100,
    query2ref=True,
    ref2query=False,
    weighting_scheme="jaccard_square",
    return_adjs=True,
)

# Pickle
import pickle

with open(
    f"{DATA_DIR}/pasca_braun_wknn_q2r.pkl",
    "wb",
) as f:
    pickle.dump(wknn_scanvi_q2r, f)


#### Label transfer ####
### Read new braun metadata ###
braun_meta = pd.read_csv(
    "/home/fleckj/projects/atlas/data/metadata/0628_braun_meta.tsv.gz",
    sep="\t",
    index_col=0,
)

braun_meta["oa_log_num_wknn_scanvi_q2r_ds_max"] = braun_meta[
    "log_num_wknn_scanvi_q2r_ds_max"
]
braun_meta["oa_observed"] = braun_meta["observed"]

# select columns not in the old metadata
braun_meta_ = braun_meta.loc[:, ~braun_meta.columns.isin(adata_braun.obs.columns)]
# add new metadata
adata_braun.obs = pd.concat([adata_braun.obs, braun_meta_], axis=1)


with open(
    f"{DATA_DIR}/pasca_braun_wknn_q2r.pkl",
    "rb",
) as f:
    wknn_scanvi_q2r = pickle.load(f)

wknn_use = wknn_scanvi_q2r

scores_class = wknn.transfer_labels(adata_braun, adata_pasca, wknn_use, "CellClass")
scores_region = wknn.transfer_labels(adata_braun, adata_pasca, wknn_use, "Subregion")
scores_cluster = wknn.transfer_labels(adata_braun, adata_pasca, wknn_use, "Clusters")
scores_neurontype = wknn.transfer_labels(
    adata_braun, adata_pasca, wknn_use, "neuron_ntt_label_only"
)
scores_ntregion = wknn.transfer_labels(
    adata_braun, adata_pasca, wknn_use, "neuron_ntt_label_only"
)

adata_pasca.obs["wknn_class"] = scores_class["best_label"]
adata_pasca.obs["wknn_region"] = scores_region["best_label"]
adata_pasca.obs["wknn_clusters"] = scores_cluster["best_label"]
adata_pasca.obs["wknn_neuron_type"] = scores_neurontype["best_label"]
adata_pasca.obs["wknn_region_neuron_type"] = scores_ntregion["best_label"]

# Save new adata
adata_pasca.write(
    "/home/fleckj/projects/atlas/scratch/v3/230605_pasca_all_v1_common_hv2k_wknn.h5ad"
)


#### Get braun umap ####
adata_braun.obsm["X_scvi"] = latent_braun
sc.pp.neighbors(adata_braun, use_rep="X_scvi", method="rapids")
sc.tl.umap(adata_braun, method="rapids")
adata_braun.obsm["X_umap_scvi"] = adata_braun.obsm["X_umap"].copy()

p = sc.pl.embedding(
    adata_braun, basis="X_umap_scvi", color=["Donor", "Subregion"], show=False
)
p[0].figure.savefig(
    "/home/fleckj/projects/atlas/plots/pasca_patscreen/braun_umap.png",
    dpi=300,
)

### Read ZHs umap ###
braun_umap = np.load(
    "/home/fleckj/projects/atlas/data/metadata/embed_umap_scanvi.npy", allow_pickle=True
)
adata_braun.obsm["X_umap_scvi"] = braun_umap

p = sc.pl.embedding(
    adata_braun, basis="X_umap_scvi", color=["Donor", "Subregion"], show=False
)
p[0].figure.savefig(
    "/home/fleckj/projects/atlas/plots/pasca_patscreen/braun_umap_zh.png",
    dpi=300,
)


# Get trans probs to plot global presence scores
trans_prob_ref = wknn.get_transition_prob_mat(adata_braun.obsm["X_scvi"], k=50)
adata_braun.obs["log_num_wknn_scanvi_q2r"] = np.log1p(
    wknn.random_walk_with_restart(
        init=np.array(wknn_scanvi_q2r.sum(axis=0)).flatten(),
        transition_prob=trans_prob_ref,
        alpha=0.1,
    )
)

sc.pl.embedding(
    adata_braun,
    basis="X_umap_scvi",
    color=["log_num_wknn_scanvi_q2r"],
    color_map="RdBu",
    vmax=np.percentile(adata_braun.obs["log_num_wknn_scanvi_q2r"], 95),
    frameon=False,
    show=False,
)
plt.savefig(
    "/home/fleckj/projects/atlas/plots/pasca_patscreen/pasca_braun_num_wknn_umap.png",
    dpi=300,
)

# Compute scores per condition
adata_pasca_full = sc.read_h5ad(
    "/home/fleckj/scratch/data/public_datasets/brain_organoids/AminPasca2023brx/230605_pasca_all_01.h5ad"
)

# Sum wknn over datasets
wknn_scanvi_q2r_per_datasets = [
    np.array(
        wknn_scanvi_q2r[adata_pasca.obs["condition"] == x, :].sum(axis=0)
    ).flatten()
    for x in adata_pasca.obs["condition"].cat.categories
]

# Smooth with random walk
wknn_scanvi_q2r_per_datasets_sm = [
    wknn.random_walk_with_restart(init=x, transition_prob=trans_prob_ref, alpha=0.1)
    for x in wknn_scanvi_q2r_per_datasets
]

# Turn to metadata
df_wknn_scanvi_q2r = pd.DataFrame(
    np.concatenate(wknn_scanvi_q2r_per_datasets_sm, axis=1),
    columns=adata_pasca.obs["condition"].cat.categories,
    index=adata_braun.obs_names,
)
df_wknn_scanvi_q2r_norm = (
    df_wknn_scanvi_q2r.apply(lambda x: np.log1p(x), axis=0)
    .apply(lambda x: np.clip(x, np.percentile(x, 1), np.percentile(x, 99)))
    .apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
)

df_wknn_scanvi_q2r_norm.to_csv(
    f"{DATA_DIR}/0711_braun_wknn_scanvi_q2r_per_datasets_sm.tsv.gz",
    sep="\t",
    compression="gzip",
)

fig, axs = plt.subplots(5, 10, figsize=(30, 11))
for i in range(df_wknn_scanvi_q2r_norm.shape[1]):
    idx_row = i // 10
    idx_col = i % 10
    adata_braun.obs["cov_score"] = df_wknn_scanvi_q2r_norm.iloc[:, i]
    sc.pl.embedding(
        adata_braun,
        ax=axs[idx_row, idx_col],
        basis="X_umap_scvi",
        color=["cov_score"],
        color_map="Greys",
        title=[adata_pasca.obs["condition"].cat.categories[i]],
        vmin=0,
        vmax=1,
        show=False,
        frameon=False,
    )

plt.savefig(
    "/home/fleckj/projects/atlas/plots/pasca_patscreen/pasca_braun_covscore_per_condition_umap.png"
)

#### Get max cov score ####
adata_braun.obs["log_num_wknn_scanvi_q2r_ds_max"] = df_wknn_scanvi_q2r_norm.max(1)
# What's detected in hBO?
adata_braun.obs["observed"] = (
    adata_braun.obs["log_num_wknn_scanvi_q2r_ds_max"] > 0.2
).astype("category")

avg_cov_cl_ref = np.array(
    [
        adata_braun.obs["log_num_wknn_scanvi_q2r_ds_max"][
            adata_braun.obs.Clusters == x
        ].mean()
        for x in np.sort(adata_braun.obs.Clusters.unique())
    ]
)
adata_braun.obs["is_in_cov_cl"] = np.isin(
    adata_braun.obs["Clusters"], np.where(avg_cov_cl_ref > 0.2)[0]
)
adata_braun.obs["is_in_cov_cl"] = adata_braun.obs["is_in_cov_cl"].astype("category")

# Plot max coverage as umap
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
sc.pl.embedding(
    adata_braun,
    ax=ax[0],
    basis="X_umap_scvi",
    color=["log_num_wknn_scanvi_q2r_ds_max"],
    color_map="RdBu",
    title="Max coverage (m)",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)
sc.pl.embedding(
    adata_braun,
    ax=ax[1],
    basis="X_umap_scvi",
    color=["observed"],
    title="Paired in hBO (m > 0.2)",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)
sc.pl.embedding(
    adata_braun,
    ax=ax[2],
    basis="X_umap_scvi",
    color=["is_in_cov_cl"],
    title="Paired in hBO (avg-cl > 0.2)",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)
plt.savefig(
    "/home/fleckj/projects/atlas/plots/pasca_patscreen/pasca_braun_max_covscore_umap.png"
)


#### Read ZHs presence scores and compare ####
oa_scores = pd.read_csv(
    "/home/fleckj/projects/atlas/data/presence_scores/presence_score_at_Braun.tsv.gz",
    sep="\t",
    index_col=0,
)

oa_scores_arr = np.array(oa_scores)

adata_braun.obs["max_cov_diff"] = (
    adata_braun.obs["log_num_wknn_scanvi_q2r_ds_max"]
    - adata_braun.obs["oa_log_num_wknn_scanvi_q2r_ds_max"]
)


# Clip at 0 to highlight improvements
adata_braun.obs["pasca_cov_gains"] = np.clip(adata_braun.obs["max_cov_diff"], 0, np.inf)


# Plot max coverage as umap
fig, ax = plt.subplots(1, 4, figsize=(20, 4))
sc.pl.embedding(
    adata_braun,
    ax=ax[0],
    basis="X_umap_scvi",
    color=["log_num_wknn_scanvi_q2r_ds_max"],
    color_map="RdBu",
    title="Max coverage pasca",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)
sc.pl.embedding(
    adata_braun,
    ax=ax[1],
    basis="X_umap_scvi",
    color=["oa_log_num_wknn_scanvi_q2r_ds_max"],
    color_map="RdBu",
    title="Max coverage OA",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)
sc.pl.embedding(
    adata_braun,
    ax=ax[2],
    basis="X_umap_scvi",
    color=["max_cov_diff"],
    color_map="RdBu",
    title="Cov diff pasca - OA",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)
sc.pl.embedding(
    adata_braun,
    ax=ax[3],
    basis="X_umap_scvi",
    color=["pasca_cov_gains"],
    color_map="Blues",
    title="Cov gains in pasca over OA",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)

plt.savefig(
    "/home/fleckj/projects/atlas/plots/pasca_patscreen/pasca_braun_diff_covscore_umap.png"
)

# Compile braun meta with umap coords
braun_umap = adata_braun.obsm["X_umap_scvi"]
braun_umap = pd.DataFrame(braun_umap, index=adata_braun.obs.index)
braun_umap.columns = ["UMAP_1", "UMAP_2"]

braun_meta = pd.concat([adata_braun.obs, braun_umap], axis=1)
braun_meta.to_csv(
    "/home/fleckj/projects/atlas/data/metadata/0711_braun_meta.tsv.gz",
    sep="\t",
    index=True,
    compression="gzip",
)

# Plot pasca to check what umap it is
sc.pl.embedding(
    adata_pasca,
    basis="X_umap",
    color=["condition"],
    frameon=False,
    size=1.5,
    sort_order=False,
    show=False,
)
plt.savefig("/home/fleckj/projects/atlas/plots/pasca_patscreen/pasca_umap.png", dpi=600)

# Compile same for pasca
pasca_umap = adata_pasca.obsm["X_umap"]
pasca_umap = pd.DataFrame(pasca_umap, index=adata_pasca.obs.index)
pasca_umap.columns = ["UMAP_1", "UMAP_2"]

pasca_meta = pd.concat([adata_pasca.obs, pasca_umap], axis=1)
pasca_meta.to_csv(
    "/home/fleckj/projects/atlas/data/metadata/0711_pasca_meta.tsv.gz",
    sep="\t",
    index=True,
    compression="gzip",
)

pasca_q2r_per_dataset = df_wknn_scanvi_q2r_norm.values
pasca_scores_adata = sc.AnnData(
    X=pasca_q2r_per_dataset.T,
    obs=pd.DataFrame(
        index=df_wknn_scanvi_q2r_norm.columns
    )
)
oa_q2r_per_dataset = oa_scores_arr
oa_scores_adata = sc.AnnData(
    X=oa_q2r_per_dataset.T,
    obs=pd.DataFrame(
        index=oa_scores.columns,
    )
)   

np.save(
    "/home/fleckj/projects/atlas/data/presence_scores/pasca_score_per_dataset.npy",
    pasca_q2r_per_dataset,
)
pasca_scores_adata.write(
    "/home/fleckj/projects/atlas/data/presence_scores/pasca_score_per_dataset.h5ad"
)

np.save(
    "/home/fleckj/projects/atlas/data/presence_scores/oa_score_per_dataset.npy",
    oa_q2r_per_dataset,
)
oa_scores_adata.write(
    "/home/fleckj/projects/atlas/data/presence_scores/oa_score_per_dataset.h5ad"
)

adata_braun_full = sc.read_h5ad(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_common.h5ad"
)

adata_braun_subs = sc.pp.subsample(
    adata_braun_full,
    n_obs=200000,
    copy=True,
    random_state=0,
)

adata_braun_subs.write(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_subs200k.h5ad"
)