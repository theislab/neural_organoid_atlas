import os
import jax
import jax.numpy as jnp

import numpy as np
import pandas as pd

# import single-cell packages
import scanpy as sc
from anndata import AnnData
from moscot.problems.time import TemporalNeuralProblem

from scanpy.tools._dpt import DPT

os.chdir("/home/fleckj/projects/atlas/scripts")


#### PATHS ####
MODEL_DIR = "/home/fleckj/projects/atlas/data/results/neural_ot/v5_marg_scores/TemporalNeuralProblem.pkl"

tnp = TemporalNeuralProblem.load(MODEL_DIR)

PROJECT_DIR = "/home/fleckj/projects/atlas/data/results/neural_ot/v5_marg_scores/"

ATLAS_DATA = "/home/fleckj/projects/atlas/scratch/v3/230510_phase3_hv2k_braun.h5ad"
adata = sc.read_h5ad(ATLAS_DATA)

## create bins for the ages
adata.obs["organoid_age_days"] = adata.obs["organoid_age_days"].astype("int")

## organize the timepoints into bins
bins = [0, 15, 30, 60, 90, 120, 150, 450]
labels = [15, 30, 60, 90, 120, 150, 450]

adata.obs["binned_ages"] = pd.cut(
    adata.obs["organoid_age_days"], bins=bins, labels=labels
)
print(adata.obs["binned_ages"].value_counts())

adata.obs["binned_ages"] = pd.to_numeric(adata.obs["binned_ages"])

latent_dim = adata.obsm["X_scpoli"].shape[1]
end_time = bins[-1]
velo_list = []
for src, tgt in list(tnp):
    solution = tnp[(src, tgt)].solution
    adata_src = adata[adata.obs["binned_ages"] == src, :]
    # Make fake anndata with velocity for source
    source = jnp.array(adata_src.obsm["X_scpoli"].copy())
    adata_velo = adata_src[:, :latent_dim].copy()
    adata_velo.layers["velocity"] = solution.pull(source) - source
    adata_velo.layers["velocity"] = np.asarray(adata_velo.layers["velocity"])
    adata_velo.layers["scpoli"] = adata_velo.obsm["X_scpoli"].copy()
    velo_list.append(adata_velo)
    # Check if target is end time
    if tgt != end_time:
        continue
    # Make fake anndata with velocity for target
    adata_tgt = adata[adata.obs["binned_ages"] == tgt, :]
    source = jnp.array(adata_tgt.obsm["X_scpoli"].copy())
    adata_velo = adata_tgt[:, :latent_dim].copy()
    adata_velo.layers["velocity"] = solution.pull(source) - source
    adata_velo.layers["velocity"] = np.asarray(adata_velo.layers["velocity"])
    adata_velo.layers["scpoli"] = adata_velo.obsm["X_scpoli"].copy()
    velo_list.append(adata_velo)

adata_velo_full = AnnData.concatenate(*velo_list, batch_key="time")
adata_velo_full.write(PROJECT_DIR + "230510_phase3_nOT_marg_score_velocity_pull.h5ad")
