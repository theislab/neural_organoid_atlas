import os
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
from matplotlib import pyplot as plt
from scib_metrics.benchmark import Benchmarker, BatchCorrection, BioConservation
from scib_metrics.nearest_neighbors import NeighborsOutput
import time
from rich import print
import faiss
from scipy.sparse import csr_matrix
import gc
import pickle


basepath = "/data/20230426_organoid_atlas_leander.dony/"
outpath = f"{basepath}output/metrics/{time.strftime('%y%m%d%H%M')}/"
os.makedirs(outpath)

adata = sc.read(os.path.join(basepath, "adatas", "230510_09_organoids_morphogens_manualannot_integrations_cleanedbenchmarking.h5ad"))

latents = [
    'X_pca_unintegrated',  
    'X_pca_rss', 
    'X_scpoli', 
    'X_benchmark_aggr_scpoli_level1', 
    'X_benchmark_aggr_scpoli_level123', 
    'X_benchmark_scpoli_level1', 
    'X_benchmark_scpoli_level12', 
    'X_benchmark_scpoli_level123', 
    'X_benchmark_scanvi_level1',
    'X_benchmark_scanvi_level12', 
    'X_benchmark_scanvi_level123', 
    'X_benchmark_scvi',
]

# Remove duplicate cells
sets = []
for l in latents:
    sets.append(set(pd.DataFrame(index=adata.obs.index, data=adata.obsm[l]).drop_duplicates(keep="first").index))
consensus_index = list(sets[0].intersection(*sets[1:]))
while len(consensus_index) != len(adata):
    print(f"Removed {adata.shape[0] - len(consensus_index)} duplicate cells ({(adata.shape[0] - len(consensus_index))/adata.shape[0]:.3f}%)")
    adata = adata[consensus_index].copy()
    sets = []
    for l in latents:
        sets.append(set(pd.DataFrame(index=adata.obs.index, data=adata.obsm[l]).drop_duplicates(keep="first").index))
    consensus_index = list(sets[0].intersection(*sets[1:]))

biocons = BioConservation(
    isolated_labels=True,
    nmi_ari_cluster_labels_leiden=True,
    nmi_ari_cluster_labels_kmeans=False,
    silhouette_label=True,
    clisi_knn=True
)
batchcorr = BatchCorrection(
    silhouette_batch=True,
    ilisi_knn=True,
    kbet_per_label=True,
    graph_connectivity=True,
    pcr_comparison=False
)

bm = Benchmarker(
    adata,
    batch_key="batch",
    label_key='snapseed_pca_unintegrated_level_123',
    embedding_obsm_keys=latents,
    pre_integrated_embedding_obsm_key="X_pca_unintegrated",
    bio_conservation_metrics=biocons, 
    batch_correction_metrics=batchcorr,
    n_jobs=-1,
)
del adata
gc.collect()

bm.prepare()

# Print latents which have a non-unique number of neighbors
for l in latents:
    if len(np.unique(np.unique(bm._emb_adatas[l].obsp["90_distances"].nonzero()[0], return_counts=True)[1])) > 1:
        x = [len(i) == 89 for i in np.split(bm._emb_adatas[l].obsp["90_distances"].indices, bm._emb_adatas[l].obsp["90_distances"].indptr[1:-1])]
        bm._emb_adatas[l] = bm._emb_adatas[l][x].copy()
        print(f"In embedding {l} cells did not have the same number of neighbors. Removed {sum([not i for i in x])} cells with less than expected number of neighbors")
print("Done checking for latents with non-unique number of neighbors")

bm.benchmark()

with open(f"{outpath}metrics_benchmark.pickle", 'wb') as handle:
    pickle.dump(bm, handle, protocol=pickle.HIGHEST_PROTOCOL)

bm.plot_results_table(show=False, save_dir=f"{outpath}scaled_results.pdf")
bm.plot_results_table(show=False, min_max_scale=False, save_dir=f"{outpath}unscaled_results.pdf")

df = bm.get_results(min_max_scale=False)
df.to_csv(f"{outpath}/metrics_unscaled.csv")

df_scaled = bm.get_results(min_max_scale=True)
df_scaled.to_csv(f"{outpath}metrics_scaled.csv")
