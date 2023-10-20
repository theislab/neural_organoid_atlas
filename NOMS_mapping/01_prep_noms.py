import warnings

warnings.filterwarnings("ignore")

import os
import sys

import numpy as np
import pandas as pd

import scanpy as sc

REF_PATH = "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_common_hv2k.h5ad"

QUERY_PATH = "/home/fleckj/scratch/data/public_datasets/brain_organoids/AminPasca2023brx/230605_pasca_all_01.h5ad"

print("Loading reference & query data")
adata_ref = sc.read(REF_PATH)
print(adata_ref)

adata_q = sc.read(QUERY_PATH)
print(adata_q)

# Get feature intersection
common_feats = np.intersect1d(adata_ref.var_names, adata_q.var_names)
q_idx = adata_q.var.index.isin(adata_ref.var_names)
adata_q = adata_q[:, q_idx]
ref_idx = adata_ref.var.index.isin(adata_q.var_names)
adata_ref = adata_ref[:, ref_idx]

# Subset to common features
adata_q.varm = dict()
adata_q.obs["batch"] = adata_q.obs["orig.ident"].astype(str).copy()

adata_q.write_h5ad(
    "/home/fleckj/scratch/data/public_datasets/brain_organoids/AminPasca2023brx/230605_pasca_all_v1_common_hv2k.h5ad"
)
