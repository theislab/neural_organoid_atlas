library(Seurat)
library(Matrix)
library(anndata)
options(Seurat.object.assay.version = "v3")

# convert data to h5ad format
data_to_h5ad.Seurat <- function(object,
                                assay = DefaultAssay(object),
                                savefile = NULL,
                                verbose = F)
{
  if (verbose)
    message("start to create the anndata object...")
  adata <- anndata::AnnData(X = t(object[[assay]]@data),
                            obs = object@meta.data,
                            var = object[[assay]]@meta.features,
                            layers = list(counts = t(object[[assay]]@counts)),
                            obsm = setNames(lapply(names(object@reductions), function(x) Embeddings(object, x)), names(object@reductions))
                            )
  if (verbose)
    message("done.")
  
  if (!is.null(savefile)){
    if (verbose)
      message(paste0("saving the anndata object to file: ",savefile,"..."))
    adata$write_h5ad(savefile)
    if (verbose)
      message("done.")
  }
  
  return(adata)
}

data_to_h5ad.default <- function(object,
                                 vars = NULL,
                                 obs = NULL,
                                 obsm = list(),
                                 layers = list(),
                                 savefile = NULL,
                                 verbose=F)
{
  if (verbose)
    message("start to create the anndata object...")
  if (is.null(vars))
    vars <- data.frame(id = colnames(object), row.names = colnames(object))
  if (is.null(obs))
    obs <- data.frame(id = rownames(object), row.names = rownames(object))
  adata <- anndata::AnnData(X = object,
                            obs = obs,
                            var = vars,
                            obsm = obsm,
                            layers = layers
                            )
  if (verbose)
    message("done.")
  
  if (!is.null(savefile)){
    if (verbose)
      message(paste0("saving the anndata object to file: ",savefile,"..."))
    adata$write_h5ad(savefile)
    if (verbose)
      message("done.")
  }
  
  return(adata)
}

data_to_h5ad <- function(object, ...) {
  UseMethod(generic = 'data_to_h5ad', object = object)
}


cells <- read.table('cell_ids.txt')[,1]
cells_samples <- tapply(cells, sapply(strsplit(cells, '-'), '[', 2), list)
cells_samples <- lapply(cells_samples, function(x) gsub('-[0-9]$','-1',x))
filtered_cells <- lapply(paste0('processed/', grep('^SAMN', list.files('./processed'), value=T),'/outs/filtered_feature_bc_matrix/barcodes.tsv.gz'), function(x){
    read.table(x)[,1]
})
names(filtered_cells) <- grep('^SAMN', list.files('./processed'), value=T)
num_overlap <- sapply(filtered_cells, function(x) sapply(cells_samples, function(y) length(intersect(x,y))))
idx_samples <- apply(num_overlap, 2, which.max)
sample_names <- setNames(names(idx_samples), idx_samples)

counts_samples <- lapply(paste0('processed/', grep('^SAMN', list.files('./processed'), value=T), '/outs/raw_feature_bc_matrix.h5'), function(x){
    Read10X_h5(x)
})
names(counts_samples) <- grep('^SAMN', list.files('./processed'), value=T)

counts_samples <- lapply(names(counts_samples), function(sample){
    counts <- counts_samples[[sample]]
    idx <- idx_samples[sample]
    colnames(counts) <- gsub('-1$', paste0('-',idx), colnames(counts))
    return(counts)
})
counts <- do.call(cbind, counts_samples)
counts_filtered <- counts[,cells]

meta_lib <- read.csv('./SraRunTable.txt', sep=',')
meta_sample <- unique(meta_lib[,c('BioSample','cell_line','AGE')])
meta <- data.frame(meta_sample, row.names=meta_sample[,'BioSample'])[sample_names[sapply(strsplit(colnames(counts_filtered),'-'), '[', 2)],]
rownames(meta) <- colnames(counts_filtered)

seurat <- CreateSeuratObject(counts=counts_filtered, meta.data=meta, project='ALICOs')
saveRDS(seurat, file='Giandomenico_2019.seurat.rds')
data_to_h5ad(seurat, savefile='Giandomenico_2019.h5ad')
