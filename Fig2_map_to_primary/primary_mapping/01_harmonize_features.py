import os
import numpy as np
import scanpy as sc
from scipy import sparse

adata = sc.read(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/human_dev_GRCh38-3.0.0_compressed.h5ad"
)

adata.obs["region_class"] = (
    adata.obs["Region"].astype(str) + "_" + adata.obs["CellClass"].astype(str)
)
adata.obs["subregion_class"] = (
    adata.obs["Subregion"].astype(str) + "_" + adata.obs["CellClass"].astype(str)
)

adata.write(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2.h5ad"
)


#### Do QC ####
adata = sc.read(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2.h5ad"
)

sc.pp.filter_cells(adata, min_genes=300)


# Run some formatting
adata.var["Gene"] = adata.var["Gene"].astype(str)
adata.var = adata.var.set_index("Gene", drop=False)
adata.var.index.name = None
adata.var_names_make_unique()

adata.write(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2qc.h5ad"
)

#### Intersect features with atlas ####
atlas = sc.read(
    "/home/fleckj/projects/atlas/scratch/v3/230510_08_organoids_labelled.h5ad"
)

atlas.var["gene_symbol"] = atlas.var["gene_symbol"].astype(str)
atlas.var = atlas.var.set_index("gene_symbol", drop=False)
atlas.var.index.name = None
atlas.var_names_make_unique()

atlas.write("/home/fleckj/projects/atlas/scratch/v3/230510_phase3.h5ad")

common_genes = np.intersect1d(adata.var.index, atlas.var.index)
np.save(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/atlasv3_devbrainv2_feature_intersect.npy",
    common_genes,
)

common_idx = np.isin(adata.var.index, common_genes)
adata_harm = adata[:, common_idx].copy()

adata_harm.write(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_common.h5ad"
)


#### Select highvar genes ####
adata_harm = sc.read(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_common.h5ad"
)

sc.pp.log1p(adata_harm)
sc.pp.highly_variable_genes(adata_harm, n_top_genes=2000, batch_key="Donor")

adata_harm.write(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_common.h5ad"
)

# Subset to hvg
adata_harm = adata_harm[:, adata_harm.var.highly_variable].copy()
adata_harm.layers["counts"] = sparse.csr_matrix(adata_harm.layers["counts"])
adata_harm.write(
    "/home/fleckj/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_common_hv2k.h5ad"
)

atlas = sc.read("/home/fleckj/projects/atlas/scratch/v3/230510_phase3.h5ad")

# Also subset atlas to devbrain hvg to save memory
atlas = atlas[:, adata_harm.var.index]
atlas.write("/home/fleckj/projects/atlas/scratch/v3/230510_phase3_hv2k_braun.h5ad")
