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

adata_oa = sc.read(ATLAS_PATH)
print(adata_oa)

adata_pasca = sc.read(QUERY_PATH)
print(adata_pasca)

# Filter Yuki_STR
adata_pasca = adata_pasca[~adata_pasca.obs["orig.ident"].isin(["Yuki_STR"]), :]

DATA_DIR = "/home/fleckj/projects/atlas/data/results/pasca_scarches/v4_OA_artur/"

wknn_scanvi_q2r = np.load(
    os.path.join(DATA_DIR, "pasca_to_organoid_wknn_scpoli_q2r.npy"), allow_pickle=True
)

# This seems to be an np.array wrapped around a sparse matrix, so we need to
# extract the sparse matrix
wknn_scanvi_q2r = wknn_scanvi_q2r.item()

# Get trans probs to plot global presence scores
trans_prob_ref = wknn.get_transition_prob_mat(adata_oa.obsm["X_scpoli"], k=50)
adata_oa.obs["log_num_wknn_scanvi_q2r"] = np.log1p(
    wknn.random_walk_with_restart(
        init=np.array(wknn_scanvi_q2r.sum(axis=0)).flatten(),
        transition_prob=trans_prob_ref,
        alpha=0.1,
    )
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
    index=adata_oa.obs_names,
)
df_wknn_scanvi_q2r_norm = (
    df_wknn_scanvi_q2r.apply(lambda x: np.log1p(x), axis=0)
    .apply(lambda x: np.clip(x, np.percentile(x, 1), np.percentile(x, 99)))
    .apply(lambda x: (x - np.min(x)) / (np.max(x) - np.min(x)))
)

df_wknn_scanvi_q2r_norm.to_csv(
    f"{DATA_DIR}/0727_oa_wknn_scanvi_q2r_per_datasets_sm.tsv.gz",
    sep="\t",
    compression="gzip",
)


fig, axs = plt.subplots(5, 10, figsize=(30, 11))
for i in range(df_wknn_scanvi_q2r_norm.shape[1]):
    idx_row = i // 10
    idx_col = i % 10
    adata_oa.obs["cov_score"] = df_wknn_scanvi_q2r_norm.iloc[:, i]
    sc.pl.embedding(
        adata_oa,
        ax=axs[idx_row, idx_col],
        basis="X_umap_scpoli",
        color=["cov_score"],
        color_map="Greys",
        title=[adata_pasca.obs["condition"].cat.categories[i]],
        vmin=0,
        vmax=1,
        show=False,
        frameon=False,
    )

plt.savefig(
    "/home/fleckj/projects/atlas/plots/pasca_patscreen/pasca_oa_covscore_per_condition_umap.png"
)


#### Get max cov score ####
adata_oa.obs["log_num_wknn_scanvi_q2r_ds_max"] = df_wknn_scanvi_q2r_norm.max(1)
# What's detected in hBO?
adata_oa.obs["observed"] = (
    adata_oa.obs["log_num_wknn_scanvi_q2r_ds_max"] > 0.2
).astype("category")

avg_cov_cl_ref = np.array(
    [
        adata_oa.obs["log_num_wknn_scanvi_q2r_ds_max"][
            adata_oa.obs["leiden_scpoli_80"] == x
        ].mean()
        for x in np.sort(adata_oa.obs["leiden_scpoli_80"].unique())
    ]
)
adata_oa.obs["is_in_cov_cl"] = np.isin(
    adata_oa.obs["leiden_scpoli_80"], np.where(avg_cov_cl_ref > 0.2)[0]
)
adata_oa.obs["is_in_cov_cl"] = adata_oa.obs["is_in_cov_cl"].astype("category")

# Plot max coverage as umap
fig, ax = plt.subplots(1, 3, figsize=(15, 4))
sc.pl.embedding(
    adata_oa,
    ax=ax[0],
    basis="X_umap_scpoli",
    color=["log_num_wknn_scanvi_q2r_ds_max"],
    color_map="RdBu",
    title="Max coverage (m)",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)
sc.pl.embedding(
    adata_oa,
    ax=ax[1],
    basis="X_umap_scpoli",
    color=["observed"],
    title="Paired in hBO (m > 0.2)",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)
sc.pl.embedding(
    adata_oa,
    ax=ax[2],
    basis="X_umap_scpoli",
    color=["is_in_cov_cl"],
    title="Paired in hBO (avg-cl > 0.2)",
    frameon=False,
    size=0.2,
    sort_order=False,
    show=False,
)
plt.savefig(
    "/home/fleckj/projects/atlas/plots/pasca_patscreen/pasca_oa_max_covscore_umap.png"
)


#### Write stuff ####
# OA umap as data frame
ua_umap = adata_oa.obsm["X_umap_scpoli"]
df_ua_umap = pd.DataFrame(
    ua_umap,
    columns=["UMAP_1", "UMAP_2"],
    index=adata_oa.obs_names,
)
df_ua_umap.to_csv(
    f"../data/metadata/0727_oa_umap.tsv",
    sep="\t",
)
