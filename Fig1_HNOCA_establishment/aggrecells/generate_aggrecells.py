import os
import glob
import re
import anndata as ad
import scanpy as sc
import biomart
import sys
import pandas as pd
import numpy as np
import math

from scipy import sparse
from sklearn.preprocessing import scale

sys.path.append('/links/groups/treutlein/USERS/zhisong_he/Tools/')
import scripts.py_util as util

# load data
os.chdir("/links/groups/treutlein/USERS/zhisong_he/Work/brain_organoid_atlas/")
adata = sc.read("data/phase3_final_0516/230510_08_organoids_labelled.h5ad")
adata.obs['publication_protocol'] = [ adata.obs['publication'][i] + ": " + adata.obs['assay_differentiation'][i] for i in range(adata.shape[0]) ]


# determine aggrecell ratios
## how-to: 1) as little aggregation as possible;
##         2) try to cap ncell per publication 50k;
##         3) try to increase avg ncount per cell over 3k
##         4) make sure avg ncell per sample over 300
##         5) cap ncell per protocol 120k;
##         
ncell_per_dataset = adata.obs['publication'].value_counts()
ncell_per_sample_per_dataset = adata.obs.groupby(['publication','bio_sample']).size()
nsample_per_dataset = ncell_per_sample_per_dataset[ncell_per_sample_per_dataset!=0].groupby('publication').size()
avg_cell_per_sample = ncell_per_dataset[nsample_per_dataset.index] / nsample_per_dataset
avg_count_per_cell = adata.obs.groupby(['publication'])['total_counts'].mean()

exp_min_aggrecell_per_sample = 300
exp_max_aggrecell_per_dataset = 35000
exp_min_count_per_cell = 3000

ratio_cell_per_sample = exp_min_aggrecell_per_sample/avg_cell_per_sample
ratio_cell_per_sample[ratio_cell_per_sample > 1] = 1
ratio_cell_per_dataset = exp_max_aggrecell_per_dataset/ncell_per_dataset[nsample_per_dataset.index]
ratio_cell_per_dataset[ratio_cell_per_dataset > 1] = 1
ratio_count_per_cell = avg_count_per_cell / exp_min_count_per_cell
ratio_count_per_cell[ratio_count_per_cell > 1] = 1
aggrecell_ratio = pd.DataFrame(np.where(ratio_cell_per_dataset < ratio_count_per_cell, ratio_cell_per_dataset, ratio_count_per_cell), index = ratio_cell_per_dataset.index)[0]
aggrecell_ratio[aggrecell_ratio < ratio_cell_per_sample] = ratio_cell_per_sample[aggrecell_ratio < ratio_cell_per_sample]
expected_ncell_per_dataset = aggrecell_ratio * ncell_per_dataset[aggrecell_ratio.index]

exp_max_aggrecell_per_protocol = 120000
ncell_per_dataset_per_protocol = adata.obs.groupby(['assay_differentiation','publication_protocol','publication']).size()
ncell_per_dataset_per_protocol = ncell_per_dataset_per_protocol[ncell_per_dataset_per_protocol!=0].reset_index().set_axis(['assay_differentiation','publication_protocol','publication','ncell'], axis=1)
ncell_per_dataset_per_protocol['exp_naggrecell'] = np.array(ncell_per_dataset_per_protocol['ncell']) * np.array(aggrecell_ratio[ncell_per_dataset_per_protocol['publication']])
expected_ncell_per_protocol = ncell_per_dataset_per_protocol.groupby('assay_differentiation').sum()
extra_ratio_protocol = exp_max_aggrecell_per_protocol / expected_ncell_per_protocol['exp_naggrecell']
extra_ratio_protocol[extra_ratio_protocol>1] = 1

aggrecell_ratio_final = pd.DataFrame(np.array(aggrecell_ratio[ncell_per_dataset_per_protocol['publication']]) * np.array(extra_ratio_protocol[ncell_per_dataset_per_protocol['assay_differentiation']]), index = ncell_per_dataset_per_protocol['publication_protocol'], columns = ['aggrecell_ratio'])
ncell_per_dataset_per_protocol['exp_naggrecell'] = np.array(ncell_per_dataset_per_protocol['ncell']) * np.array(aggrecell_ratio_final['aggrecell_ratio'])
ncell_per_dataset_per_protocol['ratio'] = np.array(aggrecell_ratio_final)

ncell_per_dataset_per_protocol.iloc[:,:5]


# split the data by publication
adata.uns['log1p']['base'] = None
adata.var_names_make_unique()
adatas_datasets = [ adata[adata.obs['publication_protocol'] == x,:].copy() for x in aggrecell_ratio_final.index ]

# build aggrecells for each batch sample
import warnings
warnings.filterwarnings("ignore")

adatas_aggr_datasets = list()
idx_aggr_datasets = list()
for adata_this in adatas_datasets:
    print('processing data set: ' + adata_this.obs.publication_protocol[0])
    ratio = aggrecell_ratio_final.loc[adata_this.obs.publication_protocol[0],'aggrecell_ratio']
    minsize = np.array([math.floor(1/ratio / 1.5), 3]).min()
    
    if ratio == 1:
        adata_this.X = adata_this.layers['counts_lengthnorm'].copy()
        adatas_aggr_datasets.append(adata_this)
        idx_cell = pd.DataFrame({'cell': adata_this.obs_names,
                                 'idx_aggrecell': range(adata_this.shape[0]),
                                 'assay_differentiation' : adata_this.obs['assay_differentiation'],
                                 'publication': adata_this.obs['publication'],
                                 'publication_protocol': adata_this.obs['publication_protocol'],
                                 'bio_sample': adata_this.obs['bio_sample'],
                                 'batch' : adata_this.obs['batch'],
                                 'aggrecell': adata_this.obs_names,
                                 'is_aggrecell': False})
        idx_aggr_datasets.append(idx_cell)
        continue
    
    samples = pd.unique(adata_this.obs['batch'])
    adatas_samples = [ adata_this[adata_this.obs['batch'] == x,:].copy() for x in samples]
    adatas_aggr_samples = list()
    idx_aggr_samples = list()
    for adata_this in adatas_samples:
        if adata_this.shape[0] < 100:
            adata_this.X = adata_this.layers['counts_lengthnorm'].copy()
            adatas_aggr_samples.append(adata_this)
            idx_cell = pd.DataFrame({'cell': adata_this.obs_names,
                                     'idx_aggrecell': range(adata_this.shape[0]),
                                     'assay_differentiation' : adata_this.obs['assay_differentiation'],
                                     'publication': adata_this.obs['publication'],
                                     'publication_protocol': adata_this.obs['publication_protocol'],
                                     'bio_sample': adata_this.obs['bio_sample'],
                                     'batch' : adata_this.obs['batch'],
                                     'aggrecell': adata_this.obs_names,
                                     'is_aggrecell': False})
            idx_aggr_samples.append(idx_cell)
            continue
        
        raw_counts = adata_this.layers['counts_lengthnorm'].copy()
        raw_var = adata_this.var
        sc.pp.highly_variable_genes(adata_this, flavor='seurat')
        adata_this = adata_this[:, adata_this.var.highly_variable]
        sc.pp.scale(adata_this, max_value=10)
        sc.tl.pca(adata_this, svd_solver='arpack')
        sc.pp.neighbors(adata_this, n_neighbors=10, n_pcs=20)
        sc.tl.leiden(adata_this)
        
        idx_aggre = util.generate_aggrecell_idx(adata_this, ratio = ratio, minsize = minsize, constraint_on_grp = "leiden")
        adata_this = ad.AnnData(X = raw_counts, obs = adata_this.obs, var = raw_var)
        adata_aggr = util.generate_aggrecells(adata_this, idx_aggre, prefix = adata_this.obs['batch'][0] + "_aggrecell_", X_is_count = True)
        adatas_aggr_samples.append(adata_aggr)
        idx_aggre['assay_differentiation'] = adata_this.obs['assay_differentiation'][0]
        idx_aggre['publication'] = adata_this.obs['publication'][0]
        idx_aggre['publication_protocol'] = adata_this.obs['publication_protocol'][0]
        idx_aggre['bio_sample'] = adata_this.obs['bio_sample'][:]
        idx_aggre['batch'] = adata_this.obs['batch'][0]
        idx_aggre['aggrecell'] = [ None if np.isnan(x) else adata_this.obs['batch'][0] + "_aggrecell_" + str(int(x)) for x in idx_aggre['idx_aggrecell'] ]
        idx_aggre['is_aggrecell'] = True
        idx_aggr_samples.append(idx_aggre)
    
    adata_aggr = ad.concat(adatas_aggr_samples)
    adatas_aggr_datasets.append(adata_aggr)
    idx_aggr = pd.concat(idx_aggr_samples)
    idx_aggr_datasets.append(idx_aggr)

adata_aggr = ad.concat(adatas_aggr_datasets)
idx_aggr = pd.concat(idx_aggr_datasets)

warnings.filterwarnings("default")


pd.concat([pd.DataFrame(adata.obs.publication.value_counts()),
           pd.DataFrame(adata_aggr.obs.publication.value_counts())], axis=1)


adata_aggr.var = adata.var.copy()
sc.pp.calculate_qc_metrics(adata_aggr, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
adata_aggr.obs['log_n_genes_by_counts'] = np.log(adata_aggr.obs.n_genes_by_counts)
adata_aggr.obs['log1p_pct_counts_mt'] = np.log1p(adata_aggr.obs.pct_counts_mt)

pd.concat([pd.DataFrame({'ncells' : adata.obs.publication.value_counts()}, index=adata.obs.publication.value_counts().index),
           pd.DataFrame({'naggrecells' : adata_aggr.obs.publication.value_counts()}, index=adata_aggr.obs.publication.value_counts().index),
           pd.DataFrame({'ngenes_per_cell' : adata.obs.groupby('publication')['n_genes_by_counts'].median()}, index=adata.obs.groupby('publication')['n_genes_by_counts'].median().index),
           pd.DataFrame({'ngenes_per_aggrecell' : adata_aggr.obs.groupby('publication')['n_genes_by_counts'].median()}, index=adata_aggr.obs.groupby('publication')['n_genes_by_counts'].median().index)], axis=1)


adata_aggr.write_h5ad("data/phase3_final_0516/phase3_aggrecells.h5ad")
idx_aggr.to_csv("data/phase3_final_0516/idx_aggrecells.tsv", sep="\t")


adata_aggr.layers['counts_lengthnorm'] = adata_aggr.X.copy()
sc.pp.normalize_total(adata_aggr)
sc.pp.log1p(adata_aggr)

adata_aggr.write_h5ad("data/phase3_final_0516/phase3_aggrecells.h5ad")
