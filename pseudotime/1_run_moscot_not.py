import warnings

warnings.filterwarnings("ignore")

import os
import sys

import jax
import jax.numpy as jnp

import numpy as np
import pandas as pd

import scvi

import scanpy as sc

from interfaces import vae_interface

from moscot.problems.time import TemporalNeuralProblem


# PATHS
PROJECT_DIR = os.path.join(
    "/home/fleckj/projects/atlas/data/results/neural_ot/v5_marg_scores/"
)
OUT_DIR = PROJECT_DIR

ATLAS_PATH = "/home/fleckj/projects/atlas/scratch/v3/230510_phase3.h5ad"

adata = sc.read(ATLAS_PATH)
print(adata)

adata.varm = dict()
adata.obs["batch"] = adata.obs["bio_sample"].astype(str).copy()

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

print("Setting up OT problem")
tnp = TemporalNeuralProblem(adata)
tnp = tnp.score_genes_for_marginals(
    gene_set_proliferation="human",
    gene_set_apoptosis="human",
    proliferation_key="prolif",
)
adata.obs["marginals"] = np.exp(4 * (adata.obs["prolif"] - adata.obs["apoptosis"]))
tnp = tnp.prepare(time_key="binned_ages", joint_attr="X_scpoli")

print("Solving OT problem")
tnp = tnp.solve(
    iterations=25000,
    compute_wasserstein_baseline=False,
    batch_size=1024,
    valid_freq=250,
    log_freq=10,
    patience=100,
    pretrain=True,
    train_size=1,
)

os.makedirs(OUT_DIR, exist_ok=True)
tnp.save(OUT_DIR, overwrite=True)
