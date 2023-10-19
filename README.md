# The Human Neural Organoid Atlas (HNOCA)

Reproducibility repository for the [Human Neural Organoid Atlas preprint](https://www.biorxiv.org/content/10.1101/2023.10.05.561097v1).

> He, Z.; Dony, L.; Fleck, J. S. et al. An integrated transcriptomic cell atlas of human neural organoids. bioRxiv 2023.10.05.561097; doi: 10.1101/2023.10.05.561097

![Abstract](https://raw.githubusercontent.com/theislab/neural_organoid_atlas/main/supplemental_files/abstract.jpg)

## Code availability

> **Note**
> :construction: We are currently collecting and uploading all notebooks used in this analysis.
> This note wil disappear as soon as all code has been made available via this page.

The "notebooks" folder in this repository contains the code used to build an analyse HNOCA.

The "integration_evaluation" directory contains the script for scVI, scPoli, and scANVI integration of HNOCA which were used in the integration benchmarking.

The "supplemental_files" folder in this repository contains:
* Snapseed input file containing our celltype hierarchy with associated marker genes
* Model parameters of our scPoli integration model
* scPoli sample embedding of HNOCA

## Data availability

An pre-release version comprised of the public portion of HNOCA is avalable for download [here](https://hmgubox2.helmholtz-muenchen.de/index.php/s/jWCCdr3SQP7Tprb/download/hnoca_pre-release_public_subset.h5ad). Download details: `hnoca_pre-release_public_subset.h5ad` (40.5 GB).

<details>
<summary>HNOCA pre-release metadata legend</summary>

### .obs keys legend
\* sfaira-controlled metadata, most controlled by the respective ontology. see here for more details: https://sfaira.readthedocs.io/en/latest/adding_datasets.html#dataset-or-observation-wise

\*\* original values provided in the datasets prior to ontology matching (only available for ontology-controlled sfaira metadata fields)

| .obs key | Description |
| --- | --- |
| `annot_level_1` | Level 1 of the final cell-type annotation hierarchy |
| `annot_level_2` | Level 1 of the final cell-type annotation hierarchy |
| `annot_level_3` | Level 1 of the final cell-type annotation hierarchy |
| `annot_level_4` | Level 1 of the final cell-type annotation hierarchy |
| `annot_region` | Final brain region annotation for neuronal cells |
| `assay_sc` | * |
| `assay_differentiation` | * |
| `assay_type_differentiation` | * |
| `bio_sample` | * |
| `cell_line` | * |
| `cell_type` | * |
| `development_stage` | * |
| `disease` | * |
| `ethnicity` | * |
| `gm` | * |
| `individual` | * |
| `organ` | * |
| `organism` | * |
| `sex` | * |
| `state_exact` | * |
| `sample_source` | * |
| `source_doi` | * |
| `tech_sample` | * |
| `treatment` | * |
| `assay_sc_original` | ** |
| `cell_line_original` | ** |
| `cell_type_original` | ** |
| `development_stage_original` | ** |
| `disease_original` | ** |
| `ethnicity_original` | ** |
| `organ_original` | ** |
| `organism_original` | ** |
| `sex_original` | ** |
| `id` | Unique sfaira sample ID |
| `suspension_type` | Single nuclei or single cells |
| `suspension_type_original` | Single nuclei or single cells |
| `obs_names_original` | Original obs annotation in the author-provided dataset |
| `organoid_age_days` | Organoid age in dayys at the time of experiment |
| `publication` | First Author, Year for the source publication of the dataset |
| `doi` | DOI for the source publication of the dataset |
| `batch` | Concatenation of the columns `id`, `bio_sample` and `tech_sample` |
| `n_genes_by_counts` | Number of genes with >0 counts |
| `log1p_n_genes_by_counts` | natural logarithm of number of genes with >0 counts +1 |
| `total_counts` | total UMI count of this cell |
| `log1p_total_counts` | Natural logarithm of the total UMI count of this cell +1 |
| `total_counts_mt` | Total mitochondrial UMI count |
| `log1p_total_counts_mt` | Natural logarithm of total mitochondrial UMI count +1 |
| `pct_counts_mt` | Percentage of UMI counts mapped to mitochindrial genes |
| `leiden_pca_unintegrated_1` | Leiden clustering (resolution=1) of unintegrated PCA-derived knn graph |
| `leiden_pca_unintegrated_80` | Leiden clustering (resolution=80) of unintegrated PCA-derived knn graph |
| `leiden_pca_rss_1` | Leiden clustering (resolution=1) of pre-integrated RSS-PCA-derived knn graph |
| `leiden_pca_rss_80` | Leiden clustering (resolution=80) of pre-integrated RSS-PCA-derived knn graph |
| `leiden_scpoli_1` | Leiden clustering (resolution=1) of integrated scPoli embedding-derived knn graph |
| `leiden_scpoli_80` | Leiden clustering (resolution=80) of integrated scPoli embedding-derived knn graph |
| `snapseed_pca_unintegrated_level_1` | Level 1 label of automatic cell type annotation (snapseed) from unintegrated PCA representation |
| `snapseed_pca_unintegrated_level_2` | Level 2 label of automatic cell type annotation (snapseed) from unintegrated PCA representation |
| `snapseed_pca_unintegrated_level_3` | Level 3 label of automatic cell type annotation (snapseed) from unintegrated PCA representation |
| `snapseed_pca_unintegrated_level_4` | Level 4 label of automatic cell type annotation (snapseed) from unintegrated PCA representation |
| `snapseed_pca_unintegrated_level_5` | Level 5 label of automatic cell type annotation (snapseed) from unintegrated PCA representation |
| `snapseed_pca_unintegrated_level_12` | Merged level 1,2 annotations from unintegrated PCA representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_pca_unintegrated_level_123` | Merged level 1,2,3 annotations from unintegrated PCA representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_pca_unintegrated_level_1234` | Merged level 1,2,3,4 annotations from unintegrated PCA representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_pca_unintegrated_level_12345` | Merged level 1,2,3,4,5 annotations from unintegrated PCA representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_pca_rss_level_1` | Level 1 label of automatic cell type annotation (snapseed) from pre-integrated RSS-PCA representation |
| `snapseed_pca_rss_level_2` | Level 2 label of automatic cell type annotation (snapseed) from pre-integrated RSS-PCA representation |
| `snapseed_pca_rss_level_3` | Level 3 label of automatic cell type annotation (snapseed) from pre-integrated RSS-PCA representation |
| `snapseed_pca_rss_level_4` | Level 4 label of automatic cell type annotation (snapseed) from pre-integrated RSS-PCA representation |
| `snapseed_pca_rss_level_5` | Level 5 label of automatic cell type annotation (snapseed) from pre-integrated RSS-PCA representation |
| `snapseed_pca_rss_level_12` | Merged level 1,2 annotations from pre-integrated RSS-PCA representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_pca_rss_level_123` | Merged level 1,2,3 annotations from pre-integrated RSS-PCA representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_pca_rss_level_1234` | Merged level 1,2,3,4 annotations from pre-integrated RSS-PCA representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_pca_rss_level_12345` | Merged level 1,2,3,4,5 annotations from pre-integrated RSS-PCA representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_scpoli_level_1` | Level 1 label of automatic cell type annotation (snapseed) from integrated scPoli representation |
| `snapseed_scpoli_level_2` | Level 2 label of automatic cell type annotation (snapseed) from integrated scPoli representation |
| `snapseed_scpoli_level_3` | Level 3 label of automatic cell type annotation (snapseed) from integrated scPoli representation |
| `snapseed_scpoli_level_4` | Level 4 label of automatic cell type annotation (snapseed) from integrated scPoli representation |
| `snapseed_scpoli_level_5` | Level 5 label of automatic cell type annotation (snapseed) from integrated scPoli representation |
| `snapseed_scpoli_level_12` | Merged level 1,2 annotations from integrated scPoli representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_scpoli_level_123` | Merged level 1,2,3 annotations from integrated scPoli representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_scpoli_level_1234` | Merged level 1,2,3,4 annotations from integrated scPoli representation (iterative replacement of NAs with higher-level labels) |
| `snapseed_scpoli_level_12345` | Merged level 1,2,3,4,5 annotations from integrated scPoli representation (iterative replacement of NAs with higher-level labels) |
| `ECM_raw` | Binary morphogen indicator for differentiation protocol: Extra-cellular matrix |
| `ROCK_inhibitor_raw` | Binary morphogen indicator for differentiation protocol: ROCK inhibitor |
| `BMP_activator_raw` | Binary morphogen indicator for differentiation protocol: BMP activator |
| `TGF_B_activator_raw` | Binary morphogen indicator for differentiation protocol: TGFβ activator |
| `TGF_B_inhibitor_raw` | Binary morphogen indicator for differentiation protocol: TGFβ inhibitor |
| `BMP_inhibitor_raw` | Binary morphogen indicator for differentiation protocol: BMP inhibitor |
| `WNT_activator_raw` | Binary morphogen indicator for differentiation protocol: WNT activator |
| `WNT_inhibitor_raw` | Binary morphogen indicator for differentiation protocol: WNT inhibitor |
| `EGF_raw` | Binary morphogen indicator for differentiation protocol: EGF |
| `FGF2_raw` | Binary morphogen indicator for differentiation protocol: FGF2 |
| `FGF8_raw` | Binary morphogen indicator for differentiation protocol: FGF8 |
| `SHH_agonist_raw` | Binary morphogen indicator for differentiation protocol: Sonic hedgehog agonist |
| `RA_raw` | Binary morphogen indicator for differentiation protocol: Retinoic acid |
| `MEK_ERK_inhibitor_raw` | Binary morphogen indicator for differentiation protocol: MEK/ERK inhibitor |
| `Notch_inhibitor_raw` | Binary morphogen indicator for differentiation protocol: Notch Inhibitor |

</details>

:construction: All data will be released here in due course.


