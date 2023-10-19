import os
import jax
import jax.numpy as jnp

import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import cellrank as cr

from scanpy.tools._dpt import DPT

os.chdir("/home/fleckj/projects/atlas/scripts")


#### FUNCS ####
def get_symmetric_transition_matrix(transition_matrix):
    sym_mat = (transition_matrix + transition_matrix.T) / 2
    # normalise transition matrix
    row_sums = sym_mat.sum(axis=1).A1
    sym_mat.data = sym_mat.data / row_sums[sym_mat.nonzero()[0]]
    return sym_mat


#### PATHS ####
PROJECT_DIR = "/home/fleckj/projects/atlas/data/results/neural_ot/v5_marg_scores/"
adata = sc.read(PROJECT_DIR + "230510_phase3_nOT_marg_score_velocity_pull.h5ad")
PLOT_DIR = "/home/fleckj/projects/atlas/data/results/neural_ot/v5_marg_scores/plots/"
os.makedirs(PLOT_DIR, exist_ok=True)

#### Put into cellrank ####
adata.obsm["X_scpoli"] = adata.layers["scpoli"].copy()
sc.pp.neighbors(adata, use_rep="X_scpoli", n_neighbors=30, method="rapids")

adata.write(PROJECT_DIR + "230510_phase3_nOT_velocity_neighbors.h5ad")

vk = cr.kernels.VelocityKernel(adata, xkey="scpoli").compute_transition_matrix()
ck = cr.kernels.ConnectivityKernel(adata).compute_transition_matrix()

adata.uns["velocity_graph"] = vk.transition_matrix

scv.set_figure_params("scvelo", transparent=True, fontsize=20)
p = scv.pl.velocity_embedding_stream(
    adata,
    basis="umap_scpoli",
    size=0.5,
    color="binned_ages",
    smooth=None,
    show=False,
)
p.figure.savefig(PLOT_DIR + "230510_phase3_nOT_k30_vk_ck_stream.png", dpi=300)
p.figure.savefig(PLOT_DIR + "230510_phase3_nOT_k30_vk_ck_stream.pdf", dpi=300)

p = scv.pl.velocity_embedding_grid(
    adata,
    basis="umap_scpoli",
    color="binned_ages",
    smooth=0.5,
    show=False,
)
p.figure.savefig(PLOT_DIR + "230510_phase3_nOT_k30_vk_ck_grid.png", dpi=300)

adata.write(PROJECT_DIR + "230510_phase3_nOT_cellrank.h5ad")


#### compute the DPT ####
dpt = DPT(adata=adata, neighbors_key="neighbors")
dpt._transitions_sym = get_symmetric_transition_matrix(adata.uns["velocity_graph"])
dpt.compute_eigen(n_comps=15, random_state=0)

adata.obsm["X_diffmap"] = dpt.eigen_basis
adata.uns["diffmap_evals"] = dpt.eigen_values

## chose a root cell (maybe we should improve this?)
# adata.uns["iroot"] = np.flatnonzero(adata.obs["binned_ages"] == "15")[100]
# Choose cell randomly
adata.uns["iroot"] = np.flatnonzero(adata.obs["binned_ages"] == "15")[42424]
sc.tl.dpt(adata)

adata.obs["pseudotime_ranks"] = adata.obs["dpt_pseudotime"].rank()
adata.obs["DC1"] = -dpt.eigen_basis[:, 1]
adata.obs["DC1_ranks"] = adata.obs["DC1"].rank()
adata.obs["DC2"] = -dpt.eigen_basis[:, 2]
adata.obs["DC2_ranks"] = adata.obs["DC2"].rank()

p = sc.pl.scatter(
    adata,
    basis="umap_scpoli",
    color=[
        "pseudotime_ranks",
        "dpt_pseudotime",
        "DC1",
        "DC1_ranks",
        "DC2",
        "DC2_ranks",
    ],
    show=False,
)
p[0].figure.savefig(PLOT_DIR + "230510_phase3_nOT_vk1_ck0_pseudotime.png", dpi=300)


p = sc.pl.scatter(
    adata,
    basis="umap_scpoli",
    color="binned_ages",
    show=False,
)
p.figure.savefig(PLOT_DIR + "230510_phase3_nOT_ages.png", dpi=300)

adata.write(PROJECT_DIR + "230510_phase3_nOT_k30_cellrank_dpt.h5ad")

# Write out the pseudotime
adata.obs[
    [
        "pseudotime_ranks",
        "dpt_pseudotime",
        "DC1",
        "DC1_ranks",
        "DC2",
        "DC2_ranks",
    ]
].to_csv(PROJECT_DIR + "230510_phase3_nOT_k30_cellrank_dpt_pseudotime.tsv", sep="\t")


#### Subset only cortex and rescale pt ####
adata = sc.read("/home/fleckj/projects/atlas/scratch/v3/230510_phase3.h5ad")
manual_annot = pd.read_csv(
    "/home/fleckj/projects/atlas/scratch/v3/tab_manual_annot_global.csv",
    index_col=0,
)

pt_meta = pd.read_csv(
    PROJECT_DIR + "230510_phase3_nOT_k30_cellrank_dpt_pseudotime.tsv",
    sep="\t",
    index_col=0,
)

ctx_class = ["cortical NPC", "cortical ExN"]
ctx_cells = manual_annot.index[manual_annot["manual_annot_ct"].isin(ctx_class)]

# Remove batch numbers (-1) from cell names in adata
clean_names = adata.obs.index.str.replace(r"-\d+$", "", regex=True)
# Add manual annot to adata
adata.obs.index = clean_names
adata.obs = pd.concat([adata.obs, manual_annot], axis=1)
pt_meta.index = pt_meta.index.str.replace(r"-\d+$", "", regex=True)
adata.obs = pd.concat([adata.obs, pt_meta], axis=1)
# Select ctx
ctx_idx = adata.obs.index.isin(ctx_cells)
adata_ctx = adata[ctx_idx, :].copy()

# UMAP
sc.pp.neighbors(adata_ctx, use_rep="X_scpoli", n_neighbors=30, method="rapids")
sc.tl.umap(adata_ctx, min_dist=0.5, spread=1.0, random_state=0, method="rapids")

adata_ctx.obs["DC1_ranks"] = adata_ctx.obs["DC1"].rank()
adata_ctx.obs["DC2_ranks"] = adata_ctx.obs["DC2"].rank()

p = sc.pl.scatter(
    adata_ctx,
    basis="umap",
    color=[
        "manual_annot_ct",
        "DC1_ranks",
        "DC2_ranks",
        "DC1",
        "DC2",
        "snapseed_scpoli_level_2",
    ],
    show=False,
)
p[0].figure.savefig(PLOT_DIR + "230510_phase3_nOT_ctx_pt.png", dpi=300)

# Feature plots
p = sc.pl.scatter(
    adata_ctx,
    basis="umap",
    color=["NEUROD6", "SOX2", "MKI67", "TBR1", "EOMES"],
    show=False,
)
p[0].figure.savefig(PLOT_DIR + "230510_phase3_nOT_ctx_features.png", dpi=300)

# Write full and ctx meta and reductions
adata.write("/home/fleckj/projects/atlas/scratch/v3/230510_phase3_pt_annot.h5ad")
adata_ctx.write(
    "/home/fleckj/projects/atlas/scratch/v3/230510_phase3_ctx_pt_annot.h5ad"
)

# Write meta with embeddings
scpoli_umap = pd.DataFrame(adata.obsm["X_umap_scpoli"], index=adata.obs.index)
meta = pd.concat([adata.obs, scpoli_umap], axis=1)
meta.to_csv(
    "/home/fleckj/projects/atlas/scratch/v3/230510_phase3_pt_annot_meta.tsv",
    sep="\t",
)

scpoli_umap_ctx = pd.DataFrame(adata_ctx.obsm["X_umap"], index=adata_ctx.obs.index)
meta_ctx = pd.concat([adata_ctx.obs, scpoli_umap_ctx], axis=1)
meta_ctx.to_csv(
    "/home/fleckj/projects/atlas/scratch/v3/230510_phase3_ctx_pt_annot_meta.tsv",
    sep="\t",
)

# Subsample 100k cells
adata_ctx_subs_100k = sc.pp.subsample(
    adata_ctx,
    n_obs=100000,
    copy=True,
    random_state=0,
)
adata_ctx_subs_100k.write(
    "/home/fleckj/projects/atlas/scratch/v3/230510_phase3_ctx_pt_annot_subs100k.h5ad"
)


#### Subset vanteral telencephalon and rescale pt ####
ge_class = ["subcortical NPC", "subcortical InN"]
ge_cells = manual_annot.index[manual_annot["manual_annot_ct"].isin(ge_class)]

# Select GE
ge_idx = adata.obs.index.isin(ge_cells)
adata_ge = adata[ge_idx, :].copy()

adata_ge.obs["DC1"] = pt_meta["DC1"]
adata_ge.obs["DC2"] = pt_meta["DC2"]
adata_ge.obs["DC1_ranks"] = adata_ge.obs["DC1"].rank()
adata_ge.obs["DC2_ranks"] = adata_ge.obs["DC2"].rank()

p = sc.pl.scatter(
    adata_ge,
    basis="umap_scpoli",
    color=[
        "manual_annot_ct",
        "DC1_ranks",
        "DC2_ranks",
        "DC1",
        "DC2",
        "snapseed_scpoli_level_2",
    ],
    show=False,
)
p[0].figure.savefig(PLOT_DIR + "230510_phase3_nOT_ge_pt.png", dpi=300)

# Feature plots
p = sc.pl.scatter(
    adata_ge,
    basis="umap_scpoli",
    color=["TAC3", "NR2F2", "ST18", "DLX5", "ISL1", "SOX2", "PAX6", "NKX2-1"],
    show=False,
)
p[0].figure.savefig(PLOT_DIR + "230510_phase3_nOT_ge_features.png", dpi=300)


### Write stuff ####
adata_ge.write("/home/fleckj/projects/atlas/scratch/v3/230510_phase3_ge_pt_annot.h5ad")

# Write meta with embeddings
scpoli_umap = pd.DataFrame(adata.obsm["X_umap_scpoli"], index=adata.obs.index)
meta = pd.concat([adata.obs, scpoli_umap], axis=1)
meta.to_csv(
    "/home/fleckj/projects/atlas/scratch/v3/230510_phase3_pt_annot_meta.tsv",
    sep="\t",
)

# Subsample 100k cells
adata_ge_subs_100k = sc.pp.subsample(
    adata_ge,
    n_obs=100000,
    copy=True,
    random_state=0,
)
adata_ge_subs_100k.write(
    "/home/fleckj/projects/atlas/scratch/v3/230510_phase3_ge_pt_annot_subs100k.h5ad"
)
