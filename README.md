# The Human Neural Organoid Atlas (HNOCA)

Reproducibility repository for the [Human Neural Organoid Atlas publication](https://doi.org/10.1038/s41586-024-08172-8).

> He, Z.; Dony, L.; Fleck, J. S. et al. An integrated transcriptomic cell atlas of human neural organoids. Nature, 2024; doi: 10.1038/s41586-024-08172-8

![Abstract](https://raw.githubusercontent.com/theislab/neural_organoid_atlas/main/supplemental_files/abstract.jpg)


> **Note**
> :construction: We are currently collecting and uploading all notebooks used in this analysis.
> This note wil disappear as soon as all code has been made available via this page.

The repository is organized based on the main figures in the manuscript the analysis is linked with:
* `Fig1_HNOCA_establishment/`: codes and Jupyter notebooks related to Fig1, about the establishment of HNOCA including data curation and integration, integration benchmarking, cell type annotation and real time-informed pseudotime analysis with moscot neural OT.
* `Fig2_map_to_primary/`: codes and Jupyter notebooks related to Fig2, about the data projection and label transfer of HNOCA to the human developing brain cell atlas ([Braun et al. 2023](https://www.science.org/doi/10.1126/science.adf1226)), as well as compositional analysis of morphogen usage.
* `Fig3_DE_to_primary/`: codes and Jupyter notebooks related to Fig3, about transcriptomic comparison (including DE and transcriptomic similarity) between HNOCA cells and their counterparts in primary atlases; also includes the codes to curate an integrated primary cortical developmental cell atlas, to dive into heterogeneity of telencephalic cells in HNOCA, and to look into glycolysis activities.
* `Fig4_Amin_mapping/`: scripts related to Fig4, about mapping the neural organoid morphogen screen ([Amin et al. 2023](https://www.biorxiv.org/content/10.1101/2023.05.31.541819v1)) data onto the HNCOA and primary data.
* `Fig5_disease_atlas/`: Jupyter notebook related to Fig5, about projecting the integrated disease-modeling brain organoid atlas to HNOCA for compositional and DE analysis.
* `Fig6_HNOCA-extended/`: codes and Juputer notebooks related to Fig6, about extending HNOCA by projecting new data.
* `supplemental_files/` contains:
  * Snapseed input file containing our celltype hierarchy with associated marker genes
  * Model parameters of our scPoli integration model
  * scPoli sample embedding of HNOCA
  * Jupyter notebook to export H5AD files for cellxgene.
