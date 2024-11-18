source('~/scripts/R/master.R')
source('~/scripts/R/colors.R')

library(Pando)

setwd('~/projects/atlas/')

np <- reticulate::import('numpy')


#### Read stuff ####
pasca_meta <- read_tsv('data/metadata/0711_pasca_meta.tsv.gz') %>% 
    rename('cell'=1) %>% 
    mutate(
        cond_group=factor(case_when(
            str_detect(condition, 'CHIR.*SAG|chir.*sag') ~ 'CHIR,SAG',
            str_detect(condition, 'SAG|sag') ~ 'SAG',
            str_detect(condition, 'CHIR|chir') ~ 'CHIR',
            str_detect(condition, 'BMP|bmp') ~ 'BMP',
            str_detect(condition, 'fgf|FGF') ~ 'FGF',
            .default = 'other'
        ), levels=c('CHIR,SAG','SAG','CHIR','BMP','FGF','other'))
    )    

oa_meta <- read_tsv('data/metadata/0627_atlas_meta.tsv.gz') %>% 
    rename('cell'=1) %>% 
    mutate(pp_annot=paste0(publication_protocol, '_', final_annot))

oa_umap <- read_tsv('data/metadata/0727_oa_umap.tsv') %>% 
    rename('cell'=1)

braun_meta <- read_tsv('data/metadata/0711_braun_meta.tsv.gz') %>% 
    rename('cell'=1)

scpoli_scarches_pasca <- anndata::read_h5ad('data/results/pasca_scarches/v4_artur/scpoli_scarches_pasca.h5ad')
scanvi_braun_pasca <- anndata::read_h5ad('data/results/pasca_scarches/v4_artur/scanvi_braun_pasca.h5ad')


#### Colors ####
region_colors_devmouse <- c(
    "Head" = "#424242",
    "Brain" = "#bdbdbd",
    "Forebrain" = "#ab1673",
    "Telencephalon" = "#7b1fa2",
    "Cortex" = "#b44a6e",
    "Hippocampus" = "#F69B97",
    "Striatum" = "#9575cd",
    "Subcortex" = "#ce93d8",
    "Hypothalamus" = "#8EC7EC",
    "Diencephalon" = "#29b6f6",
    "Thalamus" = "#29D5F6",
    "Midbrain" = "#00897b",
    "Midbrain ventral" = "#26a69a",
    "Midbrain dorsal" = "#1B8466",
    "Pons" = "#C6B161",
    "Cerebellum" = "#95CB6E",
    "Medulla" = "#f9a825",
    "Hindbrain" = "#95CB6E"
)
names(region_colors_devmouse) <- names(region_colors_devmouse) %>% str_remove_all('b')

unspec_regions <- c('Brain', 'Head', 'Forebrain', 'Head', 'Hindbrain', 'Telencephalon') %>% 
    str_remove_all('b')

zh_region_colors <- c(
    Forebrain='#A6CEE3', 
    Telencephalon='#5EA0C9', 
    Cortex='#287EB1', 
    Subcortex='#77B59A', 
    Striatum='#A1D67D', 
    Hippocampus='#5DB54B', 
    Diencephalon='#5B9E41', 
    Hypothalamus='#C59B7B', 
    Thalamus='#F47777', 
    Midbrain='#E73335', 
    "Midbrain dorsal"='#EB5037', "Midbrain ventral"='#F9A863', 
    Hindbrain='#FDA542', 
    Cerebellum='#FE8307', Pons='#E69663', Medulla='#CAB2D6',
    'non-neural'='grey'
)
names(zh_region_colors) <- names(zh_region_colors) %>% str_remove_all('b')

zh_class_colors <- c(
    'Radial glia'='#8966A9', 
    'Neuronal IPC'='#598DC7', 
    Neuroblast='#3EA9CE', 
    Neuron='#3EB6B7', 
    Glioblast='#7DBD65', 
    Oligo='#BEC326', 
    Placodes='#EBC92E', 
    Immune='#F7B82E', 
    Vascular='#F59C2A', 
    Erythrocyte='#EE7426', 
    'Neural crest'='#E6532B', 
    Fibroblast='#DC3838'
)

mol_cols <- c(
    'SAG'='#7EA184',
    'CHIR,SAG'='#A989C4',
    'CHIR'='#8FAEDE',
    'BMP'='#B4957A',
    'FGF'='#FF9D86',
    'other'='darkgrey'
)

cols_annot_level_2 <- setNames(
    c("#fee8c8","#fee0d2",
      "#9ecae1","#bcbddc","#fcbba1","#6baed6",
      "#4292c6","#9e9ac8","#fb6a4a",
      "#addd8e","#31a354","#005a32","#fdbb84",
      "#303030","#bdbdbd","#696969","#909090"),
    c('PSC','Neuroepithelium',
      'Dorsal Telencephalic NPC','Ventral Telencephalic NPC','Non-telencephalic NPC','Dorsal Telencephalic IP',
      'Dorsal Telencephalic Neuron','Ventral Telencephalic Neuron','Non-telencephalic Neuron',
      'Glioblast','Astrocyte','OPC','CP','Microglia','NC Derivatives','EC','MC')
)

bupu_zh <- colorRampPalette(c("#b2182b","#d6604d","#f4a582","#fddbc7","#f7f7f7","#d1e5f0","#92c5de","#4393c3","#2166ac"))(100)


#### Plot presence scores on braun ####
wknn_per_dataset <- read_tsv('/home/fleckj/projects/atlas/data/results/pasca_scarches/v4_OA_artur/0727_oa_wknn_scanvi_q2r_per_datasets_sm.tsv.gz') %>% 
    rename('cell'=1)

wknn_cond_mat <- wknn_per_dataset %>% 
    column_to_rownames('cell') %>% 
    as.matrix() %>% 
    Matrix(sparse=T)

pasca_protocol_agg <- Pando::aggregate_matrix(wknn_cond_mat, groups = oa_meta$publication_protocol) %>% t()
pasca_annot_agg <- Pando::aggregate_matrix(wknn_cond_mat, groups = oa_meta$final_annot) %>% t()
pasca_cluster_agg <- Pando::aggregate_matrix(wknn_cond_mat, groups = oa_meta$leiden_scpoli_1) %>% t()
pasca_protocol_annot_agg <- Pando::aggregate_matrix(wknn_cond_mat, groups = oa_meta$pp_annot) %>% t()

pasca_protocol_agg_scaled <- apply(pasca_protocol_agg, 1, scale01) %>% t()
pasca_annot_agg_scaled <- apply(pasca_annot_agg, 1, scale01) %>% t()
pasca_cluster_agg_scaled <- apply(pasca_cluster_agg, 1, scale01) %>% t()
pasca_protocol_annot_agg_scaled <- apply(pasca_protocol_annot_agg, 1, scale01) %>% t()

pasca_annot_agg %>% pheatmap::pheatmap()
pasca_protocol_agg %>% pheatmap::pheatmap()
pasca_protocol_annot_agg_scaled %>% pheatmap::pheatmap()
pasca_cluster_agg_scaled %>% pheatmap::pheatmap()



#### Plot for figure ####
#### Coembedding UMAP ####
scpoli_umap <- scpoli_scarches_pasca$obsm$X_umap
rownames(scpoli_umap) <- scpoli_scarches_pasca$obs_names
colnames(scpoli_umap) <- c('UMAP_1', 'UMAP_2')
scpoli_umap_df <- scpoli_umap %>% as_tibble(rownames='cell')
scpoli_umap_df$dataset <- ifelse(
    str_detect(scpoli_umap_df$cell, '^homo'), 'HNOCA', 'NOMS'
)


plot_df <- scpoli_umap_df %>% group_by(dataset) %>% sample_n(min(100000, n())) %>% arrange(dataset=='NOMS')
ggplot(plot_df, aes(UMAP_1, UMAP_2, color=dataset, size=dataset)) +
    geom_point() +
    scale_size_manual(values=c(1,0.5)) +
    scale_color_manual(values=c('lightgrey', '#ce93d8')) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/sfig_pasca_oa_coembedding_umap.png', width=10, height=10)   


plot_df <- scpoli_umap_df %>% filter(dataset=='HNOCA') %>% sample_n(1000000) %>% inner_join(oa_meta)
ggplot(plot_df, aes(UMAP_1, UMAP_2, color=annot_level_2)) +
    geom_point(size=0.5) +
    scale_color_manual(values=cols_annot_level_2) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/sfig_pasca_oa_coembedding_celltype_umap.png', width=10, height=10)   


scanvi_umap <- scanvi_braun_pasca$obsm$X_umap
rownames(scanvi_umap) <- scanvi_braun_pasca$obs_names
colnames(scanvi_umap) <- c('UMAP_1', 'UMAP_2')
scanvi_umap_df <- scanvi_umap %>% as_tibble(rownames='cell')
scanvi_umap_df$dataset <- ifelse(
    str_detect(scanvi_umap_df$cell, '^10X'), 'HDBCA', 'NOMS'
)


plot_df <- scanvi_umap_df %>% group_by(dataset) %>% sample_n(min(100000, n())) %>% arrange(dataset=='NOMS')
ggplot(plot_df, aes(UMAP_1, UMAP_2, color=dataset, size=dataset)) +
    geom_point() +
    scale_size_manual(values=c(1,0.5)) +
    scale_color_manual(values=c('lightgrey', '#ce93d8')) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/sfig_pasca_braun_coembedding_umap.png', width=10, height=10)   


scanvi_umap_df$cell <- str_remove(scanvi_umap_df$cell, '-1$')
plot_df <- scanvi_umap_df %>% filter(dataset!='NOMS') %>% sample_n(1000000) %>% inner_join(dplyr::select(braun_meta, !starts_with('UMAP')))
ggplot(plot_df, aes(UMAP_1, UMAP_2, color=Subregion)) +
    geom_point(size=0.5) +
    scale_color_manual(values=zh_region_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/sfig_pasca_braun_coembedding_subregion_umap.png', width=10, height=10)   




#### UMAP ####
oa_meta$max_cov_score <- wknn_cond_mat %>% rowMaxs()

plot_df <- inner_join(oa_meta, oa_umap)

ggplot(arrange(plot_df, max_cov_score), aes(UMAP_1, UMAP_2, color=max_cov_score)) +
    geom_point(size=0.02, shape=16) +
    scale_color_gradientn(colors=bupu_zh, limits=c(0,1)) +
    theme_void() + no_legend() +
    ggtitle('Pasca max coverage score')
ggsave('plots/paper/fig4/sfig_pasca_oa_max_covscore_umap.png', width=10, height=10)


#### Per condition group UMAP ####
score_mats <- map(set_names(unique(pasca_meta$cond_group)), function(x){
    conds_use <- pasca_meta %>% filter(cond_group==x) 
    score_mat <- wknn_cond_mat[,unique(conds_use$condition)] %>% rowMaxs()
})


plot_df <- inner_join(oa_meta, oa_umap)
plots <- map(names(score_mats), function(n){
    x <- score_mats[[n]]
    plot_df$tmp_cov_score <- x
    p_df <- sample_n(plot_df, 1000000)
    p <- ggplot(arrange(p_df, tmp_cov_score), aes(UMAP_1, UMAP_2, color=tmp_cov_score)) +
        geom_point(size=0.05, shape=16) +
        scale_color_gradientn(colors=bupu_zh, limits=c(0,1)) +
        theme_void() + no_legend() +
        ggtitle(n)
    return(p)
})
wrap_plots(plots)
ggsave('plots/paper/fig4/sfig_pasca_oa_max_covscore_per_cond_umap.png', width=15, height=12)



pasca_braun_scores <- anndata::read_h5ad('data/presence_scores/pasca_score_per_dataset.h5ad')
braun_scores_mat <- pasca_braun_scores$X %>% Matrix(sparse=T)
colnames(braun_scores_mat) <- braun_meta$cell
braun_scores_mat <- braun_scores_mat %>% t()

braun_score_mats <- map(set_names(unique(pasca_meta$cond_group)), function(x){
    conds_use <- pasca_meta %>% filter(cond_group==x) 
    score_mat <- braun_scores_mat[,unique(conds_use$condition)] %>% rowMaxs()
})

plot_df <- braun_meta
plots <- map(names(braun_score_mats), function(n){
    x <- braun_score_mats[[n]]
    plot_df$tmp_cov_score <- x
    p_df <- sample_n(plot_df, 1000000)
    p <- ggplot(arrange(p_df, tmp_cov_score), aes(UMAP_1, UMAP_2, color=tmp_cov_score)) +
        geom_point(size=0.05, shape=16) +
        scale_color_gradientn(colors=bupu_zh, limits=c(0,1)) +
        theme_void() + no_legend() +
        ggtitle(n)
    return(p)
})
wrap_plots(plots)
ggsave('plots/paper/fig4/sfig_pasca_braun_max_covscore_per_cond_umap.png', width=15, height=12)



#### Heatmap ####
cond_clust <- dist(pasca_protocol_agg_scaled) %>% hclust()
cond_order <- cond_clust %>% {.$labels[.$order]}

protocol_clust <- dist(t(pasca_protocol_agg_scaled)) %>% hclust()
protocol_order <- protocol_clust %>% {.$labels[.$order]}

pasca_cluster_meta <- pasca_meta %>% distinct(condition, cond_group) 
oa_cluster_meta <- oa_meta %>% distinct(publication_protocol, assay_type_differentiation, publication) 


plot_df <- pasca_protocol_agg_scaled %>% 
    as_tibble(rownames='cond') %>% 
    pivot_longer(!cond, names_to='publication_protocol', values_to='covscore') %>% 
    inner_join(oa_cluster_meta) %>% 
    left_join(pasca_cluster_meta, by=c('cond'='condition')) %>%
    mutate(
        cond=factor(cond, levels=cond_order),
        protocol=factor(publication_protocol, levels=protocol_order),
    ) 

ph <- ggplot(plot_df, aes(protocol, cond, fill=covscore)) +
    geom_tile() +
    scale_fill_gradientn(colors=grad(pals::brewer.greys, 1.5)) +
    # facet_grid(~subregion, space='free', scales='free') +
    article_text() +
    no_x_text() + no_y_text() + no_label() + no_margin() + no_legend() +
    theme(
        panel.spacing.x = unit(0.05, 'line'),
        panel.border = element_rect(linewidth=0.2)
    )

pm <- ggplot(distinct(plot_df, cond_group, cond), aes(cond, '1', fill=cond_group)) +
    geom_tile(stat='identity', width=0.8) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_manual(values=mol_cols, na.value='white') +
    theme_void() + no_legend() + coord_flip() + no_margin()

pc <- ggplot(distinct(plot_df, cond_group, cond), aes(cond, '1', label=cond)) +
    geom_text(hjust=1, size=5/ggplot2::.pt) +
    scale_y_discrete(expand=c(0,0)) +
    theme_void() + coord_flip() + no_margin()

pg <- ggplot(distinct(plot_df, assay_type_differentiation, protocol), aes(protocol, '1', fill=assay_type_differentiation)) +
    geom_tile(stat='identity', width=0.8) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_manual(values=c('black', 'grey'), na.value='white') +
    theme_void() + no_legend() + no_margin()

pt <- ggplot(distinct(plot_df, publication, protocol), aes(protocol, '1', label=publication)) +
    geom_text(angle=-90, hjust=0, size=5/ggplot2::.pt) +
    scale_y_discrete(expand=c(0,0)) +
    theme_void() + no_margin()

design <- c(
    patchwork::area(t=1, b=30, l=1, r=3),
    patchwork::area(t=1, b=30, l=4, r=5),
    patchwork::area(t=1, b=30, l=6, r=50),
    patchwork::area(t=31, b=32, l=6, r=50),
    patchwork::area(t=32, b=35, l=6, r=50)
)
(pc + pm + ph + pg + pt) + plot_layout(design=design, guides='collect')
ggsave('plots/paper/fig4/sfig_pasca_oa_coverage_heatmap.pdf', width=8, height=10, unit='cm')    





#### Heatmap with clusters / cell types ####
cond_clust <- dist(pasca_cluster_agg_scaled) %>% hclust()
cond_order <- cond_clust %>% {.$labels[.$order]}

clust_clust <- dist(t(pasca_cluster_agg_scaled)) %>% hclust()
clust_order <- clust_clust %>% {.$labels[.$order]}

pasca_cluster_meta <- pasca_meta %>% distinct(condition, cond_group) 
oa_cluster_meta <- oa_meta %>% distinct(leiden_scpoli_1, annot_level_2) %>% mutate(leiden_scpoli_1=as.character(leiden_scpoli_1))


plot_df <- pasca_cluster_agg_scaled %>% 
    as_tibble(rownames='cond') %>% 
    pivot_longer(!cond, names_to='leiden_scpoli_1', values_to='covscore') %>% 
    inner_join(oa_cluster_meta) %>%
    left_join(pasca_cluster_meta, by=c('cond'='condition')) %>%
    mutate(
        cond=factor(cond, levels=cond_order),
        cluster=factor(leiden_scpoli_1, levels=clust_order),
        annot=factor(annot_level_2, levels=names(cols_annot_level_2))
    ) 

ph <- ggplot(plot_df, aes(cluster, cond, fill=covscore)) +
    geom_tile() +
    scale_fill_gradientn(colors=grad(pals::brewer.greys, 1.5)) +
    facet_grid(~annot, space='free', scales='free') +
    article_text() +
    no_x_text() +
    no_y_text() + no_label() + no_margin() + no_legend() +
    theme(
        strip.text = element_blank(),
        panel.spacing.x = unit(0.05, 'line'),
        panel.border = element_rect(linewidth=0.2)
    )

pm <- ggplot(distinct(plot_df, cond_group, cond), aes(cond, '1', fill=cond_group)) +
    geom_tile(stat='identity', width=0.8) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_manual(values=mol_cols, na.value='white') +
    theme_void() + no_legend() + coord_flip() + no_margin()

pc <- ggplot(distinct(plot_df, cond_group, cond), aes(cond, '1', label=cond)) +
    geom_text(hjust=1, size=5/ggplot2::.pt) +
    scale_y_discrete(expand=c(0,0)) +
    theme_void() + coord_flip() + no_margin()

pg <- ggplot(distinct(plot_df, cluster, annot), aes(cluster, '1', fill=annot)) +
    geom_tile(stat='identity', width=1.1) +
    scale_y_discrete(expand=c(0,0)) +
    facet_grid(~annot, space='free', scales='free') +
    scale_fill_manual(values=cols_annot_level_2, na.value='white') +
    theme_void() + no_legend() + no_margin() +
    theme(
        strip.text = element_blank(),
        panel.spacing.x = unit(0.05, 'line')
    )

design <- c(
    # patchwork::area(t=1, b=30, l=1, r=3),
    patchwork::area(t=1, b=30, l=1, r=2),
    patchwork::area(t=1, b=30, l=3, r=50),
    patchwork::area(t=31, b=32, l=3, r=50)
)
(pm + ph + pg) + plot_layout(design=design, guides='collect')
ggsave('plots/paper/fig4/sfig_pasca_oa_cluster_coverage_heatmap.pdf', width=4, height=4, unit='cm')    
















