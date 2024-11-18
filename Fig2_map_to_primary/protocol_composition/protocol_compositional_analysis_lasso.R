source('~/scripts/R/master.R')
source('~/scripts/R/colors.R')

library(Pando)

select <- dplyr::select

setwd('~/projects/atlas/')

pub_map <- c(
    "Andersen et al. 2020" = "andersenjimena_003_d10_1016",
    "Bhaduri et al. 2020 (unguided)" = "Least Directed",
    "Bhaduri et al. 2020 (most guided)" = "Most Directed",
    "Bhaduri et al. 2020 (guided)" = "_Directed",
    "Birey et al. 2017 (hCS)" = "hCS",
    "Birey et al. 2017 (hSS)" = "hSS",
    "Bowles et al. 2021" = "bowleskathryn_003_d10_1016",
    "Fiorenzano et al. 2021" = "fiorenzanoalessandro_001_d10_1038",
    "Fleck et al. 2022" = "fleckjonassimon_001_d10_1038",
    "He el al. 2022" = "hezhisong_001_d10_1038",
    "Huang et al. 2021" = "unknown",
    "Kanton et al. 2019" = "kantonsabina_001_d10_1038",
    "Kelava et al. 2022" = "kelavaiva_001_d10_1038",
    "Khan et al. 2020" = "Khan, 2020",
    "Marton et al. 2019" = "Marton, 2019",
    "Miura et al. 2020" = "Miura, 2020",
    "Paulsen et al. 2022" = "Paulsen, 2022",
    "Pellegrini et al. 2020 (ChPO)" = "5626Tel",
    "Pellegrini et al. 2020 (CO)" = "ChP",
    "Qian et al. 2020" = "Qian, 2020",
    "Samarasinghe et al. 2021 (CO)" = "Samarasinghe",
    "Samarasinghe et al. 2021 (GEO)" = "Samarasinghe",
    "Sawada et al. 2020" = "Sawada, 2020",
    "Sloan et al. 2017" = "Sloan, 2017",
    "Trujillo et al. 2019" = "Trujillo, 2019",
    "Unpublished [Quadrato]" = "Quadrato, 2023",
    "Unpublished [Treutlein] (Jo)" = "Midbrain_Jo_",
    "Unpublished [Treutlein] (Qian)" = "Midbrain_Qian_",
    "Uzquiano et al. 2022 (Quadrato)" = "2022_09_010Quad",
    "Uzquiano et al. 2022 (Velasco)" = "cell_2022_09_010Org",
    "Velasco et al. 2019" = "Velasco, 2019",
    "Vértesy et al. 2022" = "Vértesy, 2022",
    "Xiang et al. 2019" = "Xiang, 2019",
    "Yoon et al. 2019" = "Yoon, 2019"
)



#### Load metadata ####
full_meta <- read_tsv('data/metadata/0627_atlas_meta.tsv.gz') %>% 
    rename('cell'=1) %>% 
    mutate(pp_annot=paste0(publication_protocol, '_', final_annot))

morph_mat <- read_tsv('data/morphogen_meta.tsv') 

morph_long <- morph_mat %>% 
    pivot_longer(`Birey et al. 2017 (hCS)`:`Paulsen et al. 2022`, names_to='protocol', values_to='has_mol') %>% 
    rename('molecule_group'=1, 'media'=2) %>% 
    mutate(has_mol=!is.na(has_mol), protocol=str_remove(protocol, '\\.\\.\\.\\d+'))

match_df <- map_dfr(unique(morph_long$protocol), function(x){
    det_str = pub_map[x]
    match_df <- meta_annot %>% 
        filter(str_detect(bio_sample, det_str) | str_detect(publication, det_str)) %>% 
        distinct(bio_sample, publication)
    match_df$protocol <- x
    return(match_df)
})

morph_meta <- morph_long %>% 
    distinct(molecule_group, protocol, has_mol) %>% 
    mutate(has_mol=as.numeric(has_mol)) %>% 
    group_by(protocol) %>% 
    filter(has_mol==max(has_mol)) %>% 
    pivot_wider(names_from='molecule_group', values_from=has_mol, values_fill=0)

region_meta <- full_meta %>% 
    mutate(region = case_when(
        final_region2 == 'Unspecific' ~ final_class,
        .default=final_region2
    )) %>% 
    inner_join(match_df) %>% 
    inner_join(morph_meta) %>% 
    filter(!region%in%c('PSC', 'neuroepithelium'))

ggplot(region_meta, aes(publication, fill=region)) +
    geom_bar(position='fill')

region_meta %>% write_tsv('data/results/protocol_composition/280715_region_meta.tsv')


#### PCA on region comp ####
protocol_meta <- region_meta %>% 
    distinct(bio_sample, protocol)

region_comp <- region_meta %>% 
    # filter(level_1=='neuron') %>% 
    dplyr::filter(!is.na(region)) %>% 
    dplyr::group_by(bio_sample) %>% 
    dplyr::mutate(prot_count=n()) %>% 
    dplyr::group_by(bio_sample, region) %>% 
    dplyr::summarize(
        prot_reg_count=n(),
        reg_frac=n()/prot_count
    ) %>% distinct() 

region_comp_mat <- region_comp %>%
    ungroup() %>% 
    select(bio_sample, region, reg_frac) %>% 
    distinct() %>% 
    pivot_wider(bio_sample, names_from=region, values_from=reg_frac, values_fill=0) %>% 
    column_to_rownames('bio_sample') %>% as.matrix()

region_comp_mat <- log1p(region_comp_mat*100)

protocol_pca <- prcomp_irlba(region_comp_mat, n = 2)$x
rownames(protocol_pca) <- rownames(region_comp_mat)

prot_pca_meta <- protocol_pca %>% 
    as_tibble(rownames='bio_sample') %>% 
    inner_join(protocol_meta) %>%
    inner_join(region_comp) %>%
    inner_join(morph_meta) %>%
    inner_join(as_tibble(region_comp_mat, rownames='bio_sample')) %>% 
    pivot_longer(`1-ECM`:`15-Notch inhibitor`, names_to='mol_group', values_to='has_mol') %>%
    group_by(region) %>% 
    mutate(reg_frac_scaled=zscale(reg_frac))

p1 <- ggplot(prot_pca_meta, aes(PC1, PC2, color=reg_frac_scaled)) +
    geom_point(size=1.5) +
    facet_wrap(~region) +
    scale_color_gradientn(colors=pals::magma(100))

p2 <- ggplot(arrange(prot_pca_meta, has_mol), aes(PC1, PC2, color=factor(has_mol))) +
    geom_point() +
    facet_wrap(~mol_group) +
    scale_color_manual(values=c('grey', 'black'))

p3 <- ggplot(prot_pca_meta, aes(PC1, PC2, color=protocol)) +
    geom_point()

p1 | p2
ggsave('plots/atlasv3/composition/bio_sample_region_composition_pca.png', bg='white', width=16, height=6)



ggplot(prot_pca_meta, aes(PC1, PC2, color=protocol=='Birey et al. 2017 (hSS)')) +
    geom_point()

ggplot(prot_pca_meta, aes(PC1, PC2, label=protocol, color=protocol)) +
    geom_text()



#### Train LM on region comp ####
library(glmnet)
library(glmnetUtils)

## For all regions
region_comp_mat <- region_comp %>%
    inner_join(protocol_meta) %>% 
    distinct(bio_sample, region, reg_frac) %>% 
    pivot_wider(bio_sample, names_from=region, values_from=reg_frac, values_fill=0) %>% 
    column_to_rownames('bio_sample') %>% as.matrix()

region_comp_mat <- log1p(region_comp_mat*100)
colnames(region_comp_mat) <- colnames(region_comp_mat) %>% str_replace(' ', '_')

morph_meta <- morph_long %>% 
    inner_join(protocol_meta) %>% 
    distinct(molecule_group, bio_sample, has_mol) %>% 
    mutate(has_mol=as.numeric(has_mol)) %>% 
    group_by(bio_sample) %>% 
    filter(has_mol==max(has_mol)) %>% 
    pivot_wider(names_from='molecule_group', values_from=has_mol, values_fill=0)

model_all <- map_dfr(1:ncol(region_comp_mat), function(i){
    
    y <- region_comp_mat[,i,drop=FALSE]
    
    x <- column_to_rownames(morph_meta, 'bio_sample')
    mol_names <- colnames(x) %>% 
        str_replace_all('[\\s-/]', '_') %>% str_remove('^\\d+_')
    colnames(x) <- mol_names
    
    prot_intersect <- intersect(rownames(x), rownames(y))
    
    mol_int <- levels(interaction(mol_names, mol_names, sep = ':'))
    mol_add <- paste(mol_names, collapse=' + ')
    formula_str <- reformulate(
        paste0(
            mol_add
            # paste(mol_int, collapse=' + ')
        ),
        response = colnames(region_comp_mat)[i],
        intercept = FALSE
    )
    
    model_mat <- as.data.frame(cbind(y[prot_intersect, ], x[prot_intersect, ]))
    colnames(model_mat)[1] <- colnames(region_comp_mat)[i]
    
    mframe <- model.frame(formula=formula_str, model_mat)
    
    model_fit <- cv.glmnet(formula=formula_str, data=mframe, alpha=1)
    coefs <- model_fit$glmnet.fit$beta[,model_fit$lambda==model_fit$lambda.min]
    coefs %>% 
        enframe('term', 'coef') %>% 
        mutate(region_coarse=colnames(region_comp_mat)[i]) %>% 
        return()
})


coef_mat <- model_all %>% 
    dplyr::group_by(term) %>% 
    dplyr::filter(!is.na(coef), abs(sum(coef))>0.1) %>% 
    pivot_wider(names_from=region_coarse, values_from=coef) %>% 
    column_to_rownames('term') %>% as.matrix()

mol_order <- hclust(dist(coef_mat)) %>% {.$labels[.$order]}

model_all %>% write_tsv('data/results/protocol_composition/230728_glm_a1_lmin_molecule_coefs.tsv')

plot_df <- model_all %>% 
    dplyr::group_by(term) %>% 
    dplyr::filter(!is.na(coef), abs(sum(coef))>0.1) %>% 
    mutate(term=factor(term, levels=mol_order))

p1 <- ggplot(plot_df, aes(term, coef)) +
    geom_hline(yintercept = 0) +
    geom_bar(stat='identity', color='black', fill='darkgrey') +
    facet_grid(region_coarse~.) +
    rotate_x_text(90) +
    theme(
        panel.grid.major.x = element_line(color='#EEEEEE', linewidth=0.5)
    ) +
    labs(x='Molecule group', y='Coefficient')
p1
ggsave('plots/atlasv3/composition/protocol_region_composition_lm_coefs_strict.png', bg='white', width=4, height=20)




### With interactions
model_all <- map_dfr(1:ncol(region_comp_mat), function(i){
    
    y <- region_comp_mat[,i,drop=FALSE]
    
    x <- column_to_rownames(morph_meta, 'bio_sample')
    mol_names <- colnames(x) %>% 
        str_replace_all('[\\s-/]', '_') %>% str_remove('^\\d+_')
    colnames(x) <- mol_names
    
    prot_intersect <- intersect(rownames(x), rownames(y))
    
    mol_int <- levels(interaction(mol_names, mol_names, sep = ':'))
    mol_add <- paste(mol_names, collapse=' + ')
    formula_str <- reformulate(
        paste0(
            # mol_add
            paste(mol_int, collapse=' + ')
        ),
        response = colnames(region_comp_mat)[i],
        intercept = FALSE
    )
    
    model_mat <- as.data.frame(cbind(y[prot_intersect, ], x[prot_intersect, ]))
    colnames(model_mat)[1] <- colnames(region_comp_mat)[i]
    
    mframe <- model.frame(formula=formula_str, model_mat)
    
    model_fit <- cv.glmnet(formula=formula_str, data=mframe, alpha=1)
    coefs <- model_fit$glmnet.fit$beta[,model_fit$lambda==model_fit$lambda.min]
    coefs %>% 
        enframe('term', 'coef') %>% 
        mutate(region_coarse=colnames(region_comp_mat)[i]) %>% 
        return()
})


coef_mat <- model_all %>% 
    dplyr::group_by(term) %>% 
    dplyr::filter(!is.na(coef), abs(sum(coef))>0.1) %>% 
    pivot_wider(names_from=region_coarse, values_from=coef) %>% 
    column_to_rownames('term') %>% as.matrix()

mol_order <- hclust(dist(coef_mat)) %>% {.$labels[.$order]}

plot_df <- model_all %>% 
    dplyr::group_by(term) %>% 
    dplyr::filter(!is.na(coef), abs(sum(coef))>0.1) %>% 
    mutate(term=factor(term, levels=mol_order))

p1 <- ggplot(plot_df, aes(term, coef)) +
    geom_hline(yintercept = 0) +
    geom_bar(stat='identity', color='black', fill='darkgrey') +
    facet_grid(region_coarse~.) +
    rotate_x_text(90) +
    theme(
        panel.grid.major.x = element_line(color='#EEEEEE', linewidth=0.5)
    ) +
    labs(x='Molecule group', y='Coefficient')
p1
ggsave('plots/atlasv3/composition/protocol_region_composition_lm_coefs_interaction.png', bg='white', width=8, height=8)




### PCA on mol coefs

coef_mat <- model_all %>% 
    dplyr::group_by(term) %>% 
    dplyr::filter(!is.na(coef), abs(sum(coef))>0.1) %>% 
    pivot_wider(names_from=region_coarse, values_from=coef) %>% 
    column_to_rownames('term') %>% as.matrix()
    
mol_pca <- prcomp_irlba(scale(coef_mat), n = 2)$x
rownames(mol_pca) <- rownames(coef_mat)

prot_pca_meta <- mol_pca %>% 
    as_tibble(rownames='term') 

ggplot(prot_pca_meta, aes(PC1, PC2, label=term)) +
    geom_text(size=3)
ggsave('plots/atlasv3/composition/protocol_region_composition_lm_coef_pca.png', bg='white', width=5, height=5)





