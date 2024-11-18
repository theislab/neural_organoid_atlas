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
            str_detect(condition, 'SAG') ~ 'SAG',
            str_detect(condition, 'CHIR') ~ 'CHIR',
            str_detect(condition, 'bmp') ~ 'BMP',
            str_detect(condition, 'fgf|FGF') ~ 'FGF',
            .default = 'other'
        ), levels=c('CHIR,SAG','SAG','CHIR','BMP','FGF','other'))
    )    
oa_meta <- read_tsv('data/metadata/0627_atlas_meta.tsv.gz') %>% 
    rename('cell'=1)
braun_meta <- read_tsv('data/metadata/0711_braun_meta.tsv.gz') %>% 
    rename('cell'=1)


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


#### Plot presence scores on braun ####

plot_df <- braun_meta %>% sample_n(1000000)

ggplot(plot_df, aes(UMAP_1, UMAP_2, color=log_num_wknn_scanvi_q2r_ds_max)) +
    geom_point(size=0.01, shape=16) +
    scale_color_gradientn(colors=bupu_zh, limits=c(0,1)) +
    theme_void() + no_legend() +
    ggtitle('Pasca max coverage score')
ggsave('plots/paper/fig4/pasca_braun_max_covscore_umap.png', width=10, height=8)


ggplot(plot_df, aes(UMAP_1, UMAP_2, color=pasca_cov_gains)) +
    geom_point(size=0.01, shape=16) +
    scale_color_gradientn(colors=gybu(1.5), limits=c(0,1)) +
    theme_void() + no_legend() +
    ggtitle('Pasca coverage gains over organoid atlas')
ggsave('plots/paper/fig4/pasca_braun_cov_gains_umap.png', width=10, height=8)



plot_df <- braun_meta %>% 
    group_by(Clusters) %>% 
    mutate(
        mean_gains=mean(pasca_cov_gains),
        gained_cluster=ifelse(mean_gains>0.3, as.character(Clusters), 'none')
    ) %>% 
    ungroup() %>% sample_n(1000000)

plot_df$gained_cluster[plot_df$gained_cluster=='none'] <- NA
ggplot(arrange(plot_df, !is.na(gained_cluster)), aes(UMAP_1, UMAP_2, color=gained_cluster, size=!is.na(gained_cluster))) +
    geom_point(shape=16) +
    theme_void() + no_legend() +
    scale_size_manual(values=c(0.01, 0.2)) +
    scale_color_manual(values=pals::brewer.blues(71)[30:71], na.value='grey')
ggsave('plots/paper/fig4/pasca_braun_gained_clusters_blues_umap.png', width=10, height=8)
    
ggplot(arrange(plot_df, !is.na(gained_cluster)), aes(UMAP_1, UMAP_2, color=!is.na(gained_cluster))) +
    geom_point(size=0.01, shape=16) +
    theme_void() + no_legend() +
    scale_color_manual(values=c('grey', pals::brewer.blues(5)[4]))
ggsave('plots/paper/fig4/pasca_braun_gained_clusters_blue_umap.png', width=10, height=8)

gains_meta <- braun_meta %>% 
    filter(!Subregion%in%unspec_regions, !Region%in%c('Forerain')) %>% 
    group_by(Clusters) %>% 
    mutate(
        consensus_region=mode(Region),
        consensus_subregion=case_when(
            str_detect(Subregion, 'Midrain') ~ 'Midrain',
            .default = mode(Subregion)
        ),
    ) %>% 
    summarise(
        mean_gains=mean(pasca_cov_gains), 
        sd_gain=sd(pasca_cov_gains),
        consensus_region=consensus_region[1], 
        consensus_subregion=consensus_subregion[1]
    ) %>% 
    arrange(desc(mean_gains)) %>% 
    mutate(
        Clusters=factor(Clusters, levels=unique(.$Clusters)),
        consensus_region=factor(consensus_region, levels=names(region_colors_devmouse))
    )


plot_df <- gains_meta
ggplot(plot_df, aes(Clusters, mean_gains, fill=consensus_region)) +
    geom_errorbar(aes(ymin=mean_gains-sd_gain, ymax=mean_gains+sd_gain), color='lightgrey', linewidth=0.1) +
    geom_point(shape=21, stroke=0.05, color='white', size=0.8) +
    geom_hline(yintercept=0.3, linewidth=0.3, color='darkgrey') +
    scale_fill_manual(values=zh_region_colors) +
    scale_y_continuous(limits=c(NA,1)) +
    scale_x_discrete(expand=c(0,10)) +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    
    labs(y='Gained coverage\n(in Pasca over organoid atlas)', x='Clusters')
ggsave('plots/paper/fig4/pasca_braun_gained_clusters_dots.pdf', width=8, height=4, units='cm')


ggplot(plot_df, aes(Clusters, mean_gains, fill=consensus_region)) +
    geom_errorbar(aes(ymin=mean_gains-sd_gain, ymax=mean_gains+sd_gain), color='lightgrey', linewidth=0.1) +
    geom_point(shape=21, stroke=0.05, color='white', size=0.8) +
    geom_hline(yintercept=0.3, linewidth=0.3, color='darkgrey') +
    scale_fill_manual(values=zh_region_colors) +
    scale_y_continuous(limits=c(NA,1)) +
    scale_x_discrete(expand=c(0,10)) +
    facet_grid(~consensus_region, space='free', scales='free_x') +
    scale_axis_rangeframe() + theme_rangeframe() +
    article_text() +
    no_x_text() + no_legend() +
    theme(
        panel.spacing.x = unit(0, 'line')
    ) +
    labs(y='Gained coverage\n(in Pasca over organoid atlas)', x='Clusters')
ggsave('plots/paper/fig4/pasca_braun_gained_clusters_split_dots.pdf', width=8, height=3, units='cm')



##### Assess auto annotations #####
pasca_region_colors <- c(
    'forebrain'=zh_region_colors[['Sucortex']], 'midbrain'=zh_region_colors[['Midrain']], 'hindbrain'=zh_region_colors[['Hindrain']], 
    'non-specific'='#555555', 'unknown'='grey'
)

ggplot(arrange(pasca_meta, region!='unknown'), aes(UMAP_1, UMAP_2, color=region)) +
    geom_point(size=0.3, shape=16) +
    scale_color_manual(values=pasca_region_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/pasca_region_annots_umap.png', width=10, height=8, units='cm')


ggplot(pasca_meta, aes(UMAP_1, UMAP_2, color=wknn_region)) +
    geom_point(size=0.2, shape=16) +
    scale_color_manual(values=zh_region_colors) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/pasca_region_wknn_umap.png', width=10, height=8, units='cm')

ggplot(arrange(pasca_meta, cond_group!='other'), aes(UMAP_1, UMAP_2, color=cond_group)) +
    geom_point(size=0.3, shape=16) +
    scale_color_manual(values=mol_cols) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/pasca_condition_umap.png', width=10, height=8, units='cm')



plot_df <- pasca_meta %>% 
    filter(!wknn_region %in% unspec_regions) %>% 
    mutate(
        wknn_region=case_when(
            str_detect(wknn_region, 'Midrain') ~ 'Midrain',
            .default = wknn_region
        )
    ) %>% 
    group_by(region) %>% 
    mutate(n_region=n()) %>% 
    group_by(wknn_region) %>% 
    mutate(n_wknn_region=n()) %>% 
    group_by(region, wknn_region) %>% 
    summarize(
        count = n(),
        n_region = n_region[1], 
        n_wknn_region = n_wknn_region[1],
        frac_of_annot = n() / n_region,
        frac_of_pred = n() / n_wknn_region
    ) %>% distinct() %>% 
    mutate(
        region=factor(region, levels=rev(names(pasca_region_colors))),
        wknn_region=factor(wknn_region, levels=names(region_colors_devmouse))
    )    
    

ggplot(plot_df, aes(wknn_region, region, fill=frac_of_pred)) +
    geom_tile() + 
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    scale_fill_gradientn(colors=pals::brewer.greys(100)) +
    article_text() +
    rotate_x_text(45) +
    theme(
        axis.ticks = element_line(size=0.3),
        axis.ticks.length=unit(0.05, "cm")
    ) +
    labs(y='Annotated region', x='Predicted region')
ggsave('plots/paper/fig4/pasca_annot_heatmap.pdf', width=6, height=2.7, units='cm')



plot_df <- pasca_meta %>% 
    filter(!wknn_region %in% unspec_regions) %>% 
    mutate(
        wknn_region=case_when(
            str_detect(wknn_region, 'Midrain') ~ 'Midrain',
            .default = wknn_region
        )
    ) %>% 
    mutate(
        region=factor(region, levels=rev(names(pasca_region_colors))),
        wknn_region=factor(wknn_region, levels=names(region_colors_devmouse))
    )   

comp_mat <- table(plot_df$condition, plot_df$wknn_region) %>% {./rowSums(.)}
cond_order <- comp_mat %>% dist() %>% hclust() %>% {.$labels[.$order]}

ggplot(plot_df, aes(factor(condition, levels=cond_order), fill=wknn_region)) +
    geom_bar(position='fill') +
    scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
    scale_fill_manual(values=zh_region_colors) +
    facet_grid(~cond_group, scales='free', space='free') +
    article_text() +
    no_x_text() +
    theme(
        panel.spacing = unit(0.05,'line'),
        panel.border = element_rect(linewidth=0.3)
    ) +
    labs(x='Condition', y='Fraction')
ggsave('plots/paper/fig4/pasca_pred_cond_comp_barplot.pdf', width=8, height=4, units='cm')



#### Embed the conditions by cluster-aggregated coverage ####
pasca_scores_adata <- anndata::read_h5ad('data/presence_scores/pasca_score_per_dataset.h5ad')
oa_scores_adata <- anndata::read_h5ad('data/presence_scores/oa_score_per_dataset.h5ad')

pasca_scores_agg <- Pando::aggregate_matrix(t(pasca_scores_adata$X), groups = braun_meta$Clusters) %>% t()
rownames(pasca_scores_agg) <- rownames(pasca_scores_adata$obs)

pasca_scores_agg %>% t() %>% as_tibble(rownames='Clusters') %>% write_tsv('data/Table_SX_NOMS_HDBCA_cov_per_cluster.tsv')

oa_scores_agg <- Pando::aggregate_matrix(t(oa_scores_adata$X), groups = braun_meta$Clusters) %>% t()
rownames(oa_scores_agg) <- rownames(oa_scores_adata$obs)

oa_scores_agg %>% t() %>% as_tibble(rownames='Clusters') %>% write_tsv('data/Table_SX_HNOCA_HDBCA_cov_per_cluster.tsv')

combined_scores <- rbind(pasca_scores_agg, oa_scores_agg)
# combined_scores %>% pheatmap::pheatmap()

scores_cor <- Pando::sparse_cor(t(pasca_scores_agg), t(oa_scores_agg))
# scores_cor %>% pheatmap::pheatmap()

score_pca <- prcomp_irlba(scale(combined_scores), n = 2)
rownames(score_pca$x) <- rownames(combined_scores)

plot_df <- score_pca$x %>% 
    as_tibble(rownames='sample')

ggplot(plot_df, aes(PC1, PC2, label=sample)) +
    geom_text()


#### Plot for figure ####
sample_clust <- dist(combined_scores) %>% hclust()
sample_order <- sample_clust %>% {.$labels[.$order]}

cluster_clust <- dist(t(combined_scores)) %>% hclust()
cluster_order <- cluster_clust %>% {.$labels[.$order]}

pasca_cluster_meta <- pasca_meta %>% distinct(condition, cond_group) 
oa_cluster_meta <- oa_meta %>% mutate(protocol=paste0(publication, ': ', assay_differentiation)) %>% distinct(protocol, assay_type_differentiation) 

cluster_meta <- braun_meta %>%
    filter(!Subregion%in%c(unspec_regions, 'Midrain', 'Diencephalon'), !Region%in%c('Forerain')) %>% 
    mutate(
        Subregion=case_when(
            str_detect(Subregion, 'Striatum|Sucortex') ~'Sucortex',
            str_detect(Subregion, 'Hippocampus') ~'Cortex',
            .default=Subregion
            ),
        CellClass=factor(case_when(
            str_detect(CellClass, 'Neurolast') ~'Neuroblast',
            str_detect(CellClass, 'Gliolast') ~'Glioblast',
            str_detect(CellClass, 'Firolast') ~'Fibroblast',
            .default=CellClass
        ), levels=names(zh_class_colors))
    ) %>% 
    group_by(Clusters) %>%
    summarize(
        subregion=case_when(
            !CellClass[1]%in%names(zh_class_colors)[1:4] ~ 'non-neural',
            .default=mode(Subregion)
        ),
        Clusters=as.character(Clusters)[1],
        CellClass=as.character(CellClass)[1]
    )

plot_df <- combined_scores %>% 
    as_tibble(rownames='sample') %>% 
    pivot_longer(!sample, names_to='Clusters', values_to='covscore') %>% 
    inner_join(cluster_meta) %>% 
    left_join(pasca_cluster_meta, by=c('sample'='condition')) %>% 
    left_join(oa_cluster_meta, by=c('sample'='protocol')) %>% 
    mutate(
        sample=factor(sample, levels=sample_order),
        cluster=factor(Clusters, levels=cluster_order),
        origin=ifelse(str_detect(sample, 'doi'), 'OA', 'Pasca'),
        subregion=factor(subregion, levels=names(zh_region_colors)),
        CellClass=factor(CellClass, levels=names(zh_class_colors))
    ) 

ph <- ggplot(plot_df, aes(cluster, sample, fill=covscore)) +
    geom_tile() +
    scale_fill_gradientn(colors=grad(pals::brewer.greys, 1.5)) +
    facet_grid(~subregion, space='free', scales='free') +
    article_text() +
    no_x_text() + no_y_text() + no_label() + no_margin() +
    theme(
        panel.spacing.x = unit(0.05, 'line'),
        panel.border = element_rect(linewidth=0.2)
    )

pr <- ggplot(distinct(plot_df, subregion, cluster), aes(cluster, 1, fill=subregion)) +
    geom_tile(stat='identity') +
    scale_y_continuous(expand=c(0,0)) +
    facet_grid(~subregion, space='free', scales='free') +
    scale_fill_manual(values=zh_region_colors) +
    theme_void() + no_legend() + no_margin() +
    theme(
        strip.text = element_blank(),
        panel.spacing.x = unit(0.05, 'line'),
    )

pc <- ggplot(distinct(plot_df, CellClass, cluster, subregion), aes(cluster, '1', fill=CellClass)) +
    geom_tile(stat='identity', width=0.8) +
    scale_y_discrete(expand=c(0,0)) +
    facet_grid(~subregion, space='free', scales='free') +
    scale_fill_manual(values=zh_class_colors) +
    theme_void() + no_margin() + 
    article_text() + no_x_text() + no_y_text() +
    theme(
        strip.text = element_blank(),
        panel.spacing.x = unit(0.05, 'line'),
    )
    
po <- ggplot(distinct(plot_df, origin, sample), aes(sample, '1', fill=origin)) +
    geom_tile(stat='identity', width=0.8) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_manual(values=c('#8EC7EC', '#ce93d8')) +
    theme_void() + no_legend() + coord_flip() + no_margin()

pm <- ggplot(distinct(plot_df, cond_group, sample), aes(sample, '1', fill=cond_group)) +
    geom_tile(stat='identity', width=0.8) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_manual(values=mol_cols, na.value='white') +
    theme_void() + coord_flip() + no_margin()

pg <- ggplot(distinct(plot_df, assay_type_differentiation, sample), aes(sample, '1', fill=assay_type_differentiation)) +
    geom_tile(stat='identity', width=0.8) +
    scale_y_discrete(expand=c(0,0)) +
    scale_fill_manual(values=c('black', 'grey'), na.value='white') +
    theme_void() + coord_flip() + no_margin()
pg

design <- c(
    patchwork::area(t=1, b=30, l=1, r=2),
    patchwork::area(t=1, b=30, l=3, r=4),
    patchwork::area(t=1, b=30, l=5, r=6),
    patchwork::area(t=1, b=30, l=7, r=50),
    patchwork::area(t=31, b=32, l=7, r=50),
    patchwork::area(t=33, b=34, l=7, r=50)
)
(po + pm + pg + ph + pr + pc) + plot_layout(design=design, guides='collect')
ggsave('plots/paper/fig4/pasca_braun_oa_coverage_heatmap.pdf', width=10, height=3.8, unit='cm')    
    

library(ggraph)
ggraph(sample_clust) + geom_edge_elbow()
ggsave('plots/paper/fig4/pasca_braun_oa_coverage_dend.pdf', width=10, height=3.8, unit='cm')    




#### Feature plots ####
pasca_ad <- anndata::read_h5ad('~/scratch/data/public_datasets/brain_organoids/AminPasca2023brx/230605_pasca_all_01.h5ad')
pasca_mat <- t(pasca_ad$X) %>% Matrix::Matrix(sparse=T)

pasca_ad$obs$nFeature_RNA <- NULL

pasca_srt <- CreateSeuratObject(
    pasca_mat,  meta.data=pasca_ad$obs
)


dim_plot(pasca_srt, group.by=c('seurat_clusters'), order=T, pt.size=0.7, label=T)

ggplot(pasca_meta, aes(UMAP_1, UMAP_2, color=seurat_clusters==6)) +
    geom_point(size=0.3, shape=16) +
    theme_void() + no_legend()

ggplot(pasca_meta, aes(UMAP_1, UMAP_2, color=seurat_clusters==36)) +
    geom_point(size=0.3, shape=16) +
    theme_void() + no_legend()

ggplot(pasca_meta, aes(UMAP_1, UMAP_2, color=seurat_clusters==5)) +
    geom_point(size=0.3, shape=16) +
    theme_void() + no_legend()



umap_mat <- pasca_ad$obsm['X_umap']$X_umap
rownames(umap_mat) <- rownames(pasca_ad$obs)
pasca_srt[['umap']] <- CreateDimReducObject(umap_mat, key='UMAP_')

gset <- c('ESRRB', 'CA8', 'TFAP2A', 'UBASH3B')
gset <- c('CA8', 'SKOR2')
pasca_srt$coex_score <- pasca_srt[['RNA']][gset,] %>% colMins()
feature_plot(pasca_srt, features=c('coex_score'), order=T, pt.size=0.7)
ggsave('plots/paper/fig4/pasca_medulla_purk_ca8_skor2_coex_umap.png', width=10, height=8)

plot_df <- pasca_srt@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(pasca_meta, by='cell') %>% 
    filter(coex_score>0)

ggplot(plot_df, aes('1', fill=cond_group)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=mol_cols) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/pasca_medulla_purk_ca8_skor2_cond_frac.pdf', width=3, height=4)




gset <- c('LHX6', 'ACKR3', 'MPPED1')
pasca_srt$coex_score <- pasca_srt[['RNA']][gset,] %>% colMins()
feature_plot(pasca_srt, features=c('coex_score'), order=T, pt.size=0.7)
ggsave('plots/paper/fig4/pasca_telen_lhx6_coex_umap.png', width=10, height=8)

plot_df <- pasca_srt@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(pasca_meta, by='cell') %>% 
    filter(coex_score>0)

ggplot(plot_df, aes('1', fill=cond_group)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=mol_cols) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/pasca_telen_lhx6_cond_frac.pdf', width=3, height=4)

ggplot(plot_df, aes('1', fill=condition.x)) +
    geom_bar(position='fill') +
    # scale_fill_manual(values=mol_cols) +
    theme_void()


gset <- c('LHX8','CA10')
pasca_srt$coex_score <- pasca_srt[['RNA']][gset,] %>% colMins()
feature_plot(pasca_srt, features=c('coex_score'), order=T, pt.size=0.7)
ggsave('plots/paper/fig4/pasca_dien_lhx8_coex_umap.png', width=10, height=8)

gset <- c('GATA3', 'CDH8')
pasca_srt$coex_score <- pasca_srt[['RNA']][gset,] %>% colMins()
feature_plot(pasca_srt, features=c('coex_score'), order=T, pt.size=0.7)
ggsave('plots/paper/fig4/pasca_pons_gata3_coex_umap.png', width=10, height=8)

gset <- c('TH', 'EN1')
pasca_srt$coex_score <- pasca_srt[['RNA']][gset,] %>% colMins()
feature_plot(pasca_srt, features=c('coex_score'), order=T, pt.size=0.7)
ggsave('plots/paper/fig4/pasca_midbrain_th_en1_coex_umap.png', width=10, height=8)

plot_df <- pasca_srt@meta.data %>% as_tibble(rownames='cell') %>% 
    inner_join(pasca_meta, by='cell') %>% 
    filter(coex_score>0)

ggplot(plot_df, aes('1', fill=cond_group)) +
    geom_bar(position='fill') +
    scale_fill_manual(values=mol_cols) +
    theme_void() + no_legend()
ggsave('plots/paper/fig4/pasca_midbrain_th_en1_cond_frac.pdf', width=3, height=4)




feature_plot(pasca_srt, features=c('OLIG1', 'MBP', 'SOX10'), order=T, pt.size=0.7)



plot_df <- braun_meta %>% sample_n(50000)

gains_meta %>% filter(consensus_region=='Telencephalon') %>% top_n(3,mean_gains)
plot_df$this <- plot_df$Clusters==364
ggplot(arrange(plot_df, this), aes(UMAP_1, UMAP_2, color=this)) +
    geom_point()

gains_meta %>% filter(consensus_region=='Cereellum') %>% top_n(3,mean_gains)
plot_df$this <- plot_df$Clusters==524
ggplot(arrange(plot_df, this), aes(UMAP_1, UMAP_2, color=this)) +
    geom_point()

gains_meta %>% filter(consensus_region=='Diencephalon') %>% top_n(3,mean_gains)
plot_df$this <- plot_df$Clusters==443
ggplot(arrange(plot_df, this), aes(UMAP_1, UMAP_2, color=this)) +
    geom_point()

gains_meta %>% filter(consensus_region=='Pons') %>% top_n(3,mean_gains)
plot_df$this <- plot_df$Clusters==450
ggplot(arrange(plot_df, this), aes(UMAP_1, UMAP_2, color=this)) +
    geom_point()

gains_meta %>% filter(consensus_region=='Midrain') %>% top_n(3,mean_gains)
plot_df$this <- plot_df$Clusters==385
ggplot(arrange(plot_df, this), aes(UMAP_1, UMAP_2, color=this)) +
    geom_point()



#### Feature plots in Braun ####
braun_ad <- anndata::read_h5ad('~/scratch/data/public_datasets/primary_brain/BraunLinnarsson2022/braun_2022_fetal_brain_v2_subs200k.h5ad')
braun_mat <- t(braun_ad$X) %>% Matrix::Matrix(sparse=T)

braun_srt <- CreateSeuratObject(
    braun_mat,  meta.data=braun_ad$obs
)

umap_mat <- braun_meta %>% dplyr::select(cell, UMAP_1, UMAP_2) %>% column_to_rownames('cell') %>% as.matrix()
umap_mat <- umap_mat[rownames(braun_ad$obs), ]
braun_srt[['umap']] <- CreateDimReducObject(umap_mat, key='UMAP_')

gset <- c('ESRRB', 'CA8', 'TFAP2A', 'UBASH3B')
gset <- c('CA8', 'SKOR2')
braun_srt$coex_score <- braun_srt[['RNA']][gset,] %>% colMins()
feature_plot(braun_srt, features=c('coex_score'), order=T, pt.size=0.7, raster=F)
ggsave('plots/paper/fig4/braun_medulla_purk_ca8_skor2_coex_umap.png', width=10, height=8)

gset <- c('LHX6', 'ACKR3', 'MPPED1')
braun_srt$coex_score <- braun_srt[['RNA']][gset,] %>% colMins()
feature_plot(braun_srt, features=c('coex_score'), order=T, pt.size=0.7, raster=F)
ggsave('plots/paper/fig4/braun_telen_lhx6_coex_umap.png', width=10, height=8)

gset <- c('LHX8','CA10')
braun_srt$coex_score <- braun_srt[['RNA']][gset,] %>% colMins()
feature_plot(braun_srt, features=c('coex_score'), order=T, pt.size=0.7, raster=F)
ggsave('plots/paper/fig4/braun_dien_lhx8_coex_umap.png', width=10, height=8)

gset <- c('GATA3', 'CDH8')
braun_srt$coex_score <- braun_srt[['RNA']][gset,] %>% colMins()
feature_plot(braun_srt, features=c('coex_score'), order=T, pt.size=0.7, raster=F)
ggsave('plots/paper/fig4/braun_pons_gata3_coex_umap.png', width=10, height=8)

gset <- c('TH', 'EN1')
braun_srt$coex_score <- braun_srt[['RNA']][gset,] %>% colMins()
feature_plot(braun_srt, features=c('coex_score'), order=T, pt.size=0.7, raster=F)
ggsave('plots/paper/fig4/braun_midbrain_th_en1_coex_umap.png', width=10, height=8)












