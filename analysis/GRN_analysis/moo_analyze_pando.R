source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/perturbator/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

library(Signac)
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M, 
        lib.loc = '/local1/USERS/tomasgomes/multiome_analysis/')
library(Pando)
library(doParallel)

registerDoParallel(40)

setwd('~/projects/axolotl/')


#### Read Stuff and format ####
mome_atac <- read_rds('data/muo_glut_v1.3.2pando_ATAC_links10m_srt.rds')
rna_expr <- t(Seurat::GetAssayData(mome_atac, assay='RNA', slot='data'))

new_trajectories <- read_csv('ref_glut_dat.csv') %>% 
    mutate(cell=str_replace_all(newcellnames, '_', '-')) %>% 
    filter(cell%in%colnames(mome_atac))
    
fate_df <- new_trajectories %>% 
    select(newcellnames, fate, new_pseudotime) %>% distinct() %>% 
    mutate(
        newcellnames=str_replace_all(newcellnames, '_', '-')
    ) %>% 
    pivot_wider(names_from = 'fate', values_from = 'new_pseudotime') %>% 
    column_to_rownames('newcellnames')

group_df <- new_trajectories %>% 
    select(newcellnames, group) %>% distinct() %>% 
    mutate(group=paste0('group_', group)) %>% 
    mutate(
        val=group,
        newcellnames=str_replace_all(newcellnames, '_', '-')
    ) %>% 
    pivot_wider(names_from = 'group', values_from = 'val') %>% 
    column_to_rownames('newcellnames')

mome_atac <- AddMetaData(mome_atac, fate_df)
mome_atac <- AddMetaData(mome_atac, group_df)

dim_plot(mome_atac, reduction='rnaumap', group.by=c('group_hippocampus', 'group_lc'))

mome_atac %>% write_rds('data/muo_glut_v1.4annot_links10m_srt.rds')

mome_atac <- find_modules(
    mome_atac,
    network = 'glm_network',
    p_thresh = 0.05,
    nvar_thresh = 2
)

glm_modules <- NetworkModules(mome_atac)
glm_net <- glm_modules@meta

glm_coefs <- coef(mome_atac, network='glm_network')
glmnet_coefs <- coef(mome_atac, network='glmnet_network')
bagridge_coefs <- coef(mome_atac, network='bagging_ridge_network')

glm_net <- glm_coefs %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
    filter(padj<0.05) %>% 
    group_by(tf, target) %>% 
    filter(pval==min(pval)) %>% 
    summarize(
        estimate=mean(estimate),
        pval=min(pval),
        padj=min(padj),
        regions=region[1]
    )
glmnet_net <- glmnet_coefs %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
    filter(estimate>0) %>% 
    group_by(tf, target) %>% 
    filter(estimate==max(abs(estimate))) %>% 
    summarize(
        estimate=mean(estimate),
        regions=region[1]
    )
bagridge_net <- bagridge_coefs %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
    filter(padj<0.05) %>% 
    group_by(tf, target) %>% 
    filter(pval==min(pval)) %>% 
    summarize(
        estimate=mean(estimate),
        pval=min(pval),
        padj=min(padj),
        regions=region[1]
    )


net_use <- bagridge_net


#### Global motif matches ####
#### Global motif matches ####
# Initiate GRN object and select candidate regions
pfm <- read_rds('~/resources/DB/JASPAR/JASPAR2020_vertebrates_motifs.rds')
motif_df <- read_tsv('data/JASPAR2020_vertebrates_motif2tf.tsv')

mome_atac2 = initiate_grn(
    mome_atac, 
    peak_assay = 'ATAC', 
    rna_assay = 'RNA', 
    exclude_exons = F
)

# Scan candidate regions for TF binding motifs
mome_atac2 = find_motifs(
    mome_atac2, 
    pfm = pfm, 
    motif_tfs = motif_df,
    genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M
)

global_motifs <- mome_atac2@grn@regions@motifs@data
global_motifs %>% write_rds('data/grn/global_motif_matches.rds')

motif_counts <- colSums2(global_motifs)
names(motif_counts) <- colnames(global_motifs)
motif_counts <- motif_counts[motif_df$motif]
names(motif_counts) <- motif_df$tf

motif_count_df <- motif_counts %>% 
    enframe('tf', 'motif_count') %>% 
    arrange(desc(motif_count)) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf)))

tf_count <- net_use %>% 
    group_by(tf) %>% 
    mutate(gene_count=n()) %>% 
    distinct(tf, gene_count) %>% 
    inner_join(motif_count_df)

p1 <- ggplot(tf_count, aes(gene_count, motif_count, label=tf)) +
    geom_text()
p2 <- ggplot(tf_count, aes(log10(gene_count), log10(motif_count), label=tf)) +
    geom_text()

p1 / p2


#### Get avg expression for genes ####
region_summary <- Pando::aggregate_matrix(rna_expr[, union(net_use$tf, net_use$target)], groups=mome_atac$pred_regions_all)
subclass_summary <- Pando::aggregate_matrix(rna_expr[, union(net_use$tf, net_use$target)], groups=mome_atac$subclasses)

region_summary_df <- region_summary %>% t() %>% 
    as_tibble(rownames='gene') 

subclass_summary_df <- subclass_summary %>% t() %>% 
    as_tibble(rownames='gene') 

gene_scores <- inner_join(region_summary_df, subclass_summary_df)


### Get coex and umap ####
gene_cor <- Pando::sparse_cor(rna_expr[, union(net_use$tf, net_use$target)])

reg_mat <- net_use %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
    # filter(target%in%.$tf) %>%
    distinct(target, tf, estimate) %>%
    pivot_wider(names_from=tf, values_from=estimate, values_fill=0) %>% 
    column_to_rownames('target') %>% as.matrix() %>% Matrix::Matrix(sparse=T)
reg_factor_mat <- abs(reg_mat) + 1

weighted_coex_mat <- gene_cor[rownames(reg_factor_mat), colnames(reg_factor_mat)] * sqrt(reg_factor_mat)
weighted_coex_mat <- as.matrix(gene_cor[rownames(reg_factor_mat), colnames(reg_factor_mat)])
weight_coex_umap <- uwot::umap(weighted_coex_mat, n_neighbors=5)
rownames(weight_coex_umap) <- rownames(weighted_coex_mat)
colnames(weight_coex_umap) <- c('UMAP1', 'UMAP2')


#### Plot network ####
weight_coex_meta <- weight_coex_umap %>% 
    as_tibble(rownames='gene') %>% 
    left_join(gene_scores)


glm_graph <- as_tbl_graph(net_use) %>% 
    activate(edges) %>% 
    mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>% 
    activate(nodes) %>% 
    mutate(
        central_pr=centrality_pagerank(weights = estimate),
        central_betw=centrality_betweenness(),
        central_eig=centrality_eigen(),
        central_deg=centrality_degree(),
        outdegree=centrality_degree(mode='out'),
        indegree=centrality_degree(mode='in')
    ) %>% 
    inner_join(weight_coex_meta, by=c('name'='gene')) %>% 
    activate(edges) %>%
    filter(padj<1e-2)


p1 <- ggraph(glm_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) + 
    geom_node_point(aes(fill=dorsal, size=central_eig), shape=21, stroke=0.2) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void()  + ggtitle('dorsal') + no_legend()

p2 <- ggraph(glm_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) + 
    geom_node_point(aes(fill=lateral, size=central_eig), shape=21, stroke=0.2) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void() + ggtitle('lateral') + no_legend()

p3 <- ggraph(glm_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) + 
    geom_node_point(aes(fill=medial, size=central_eig), shape=21, stroke=0.2) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void() + ggtitle('medial') + no_legend()

p4 <- ggraph(glm_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) + 
    geom_node_point(aes(fill=Ependymal, size=central_eig), shape=21, stroke=0.2) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void() + ggtitle('Ependymal') + no_legend()

p5 <- ggraph(glm_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(alpha=-log10(padj)), color='darkgray', width=0.2) + 
    geom_node_point(aes(fill=NPC, size=central_eig), shape=21, stroke=0.2) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void() + ggtitle('NPC') + no_legend()


(p1 | p2 | p3) / (p4 | p5)



ggraph(glm_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(alpha=-log10(padj)), color='grey', width=0.1) + 
    geom_node_point(shape=21, stroke=0.2, color='black', fill='darkgrey') +
    geom_node_text(aes(label=name), size=2, repel=T, max.overlaps=100) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void()

ggsave('plots/glut_grn_umap.png', width=20, height=20)


glm_graph %>% write_rds('data/grn/bagridge_graph_1m_1.rds')
# glm_graph <- read_rds('data/grn/bagridge_graph_1m_1.rds')



#### Gene DE between regions ####
mome_region_de <- de(mome_atac, groups='region', assay='RNA')
perc_expr <- Pando::aggregate_matrix(t(mome_atac[['RNA']]@counts>0), groups=mome_atac$region)
perc_expr_df <- perc_expr %>% 
    as_tibble(rownames='region') %>% 
    pivot_longer(!region, names_to='gene', values_to='perc_expr') %>% 
    pivot_wider(names_from = region, values_from = perc_expr, names_prefix = 'perc_expr_')



#### Prep object for DA ####
mome_atac@active.assay <- 'ATAC'
mome_atac_test <- DietSeurat(mome_atac, assays = c('ATAC'))
mome_atac_peaks_bin <- as(mome_atac@assays$ATAC@counts>0, 'dgCMatrix')
mome_atac_test[['ATAC_bin']] <- CreateAssayObject(mome_atac_peaks_bin)
peaks_test <- unique(mome_atac@assays$ATAC@links$peak)


#### DA between regions ####
regions <- set_names(unique(mome_atac$region))
region_da <- map_dfr(regions, function(reg){
    print(reg)
    mome_atac_test$test_var <- mome_atac_test$region == reg
    test_df <- lr_de(
        object = mome_atac_test,
        features_use = peaks_test,
        test_var = 'test_var',
        covariates = 'nFeature_ATAC',
        family = 'binomial',
        assay = 'ATAC_bin'
    )
    detection_rates <- Pando::aggregate_matrix(t(mome_atac_test[['ATAC_bin']]@counts[peaks_test, ]), mome_atac_test$test_var)
    detrate_df <- tibble(
        feature = colnames(detection_rates),
        detect_self = as.numeric(detection_rates['TRUE', ]),
        detect_other = as.numeric(detection_rates['FALSE', ])
    )
    return(inner_join(test_df, detrate_df))
}, .id='group')

mome_atac <- aggregate_assay(mome_atac, group_name='cellclusters', assay='ATAC_bin')
cluster_detection <- mome_atac@assays$ATAC_bin@misc$summary$cellclusters
mome_atac <- aggregate_assay(mome_atac, group_name='region', assay='ATAC_bin')
region_detection <- mome_atac@assays$ATAC_bin@misc$summary$region

detection_rate <- tibble(
    feature = colnames(cluster_detection),
    max_cluster_detection = colMaxs(cluster_detection),
    max_region_detection = colMaxs(region_detection)
)

region_da <- region_da %>% 
    mutate(
        detect_ratio=detect_self/detect_other,
        padj=p.adjust(pval, method='fdr')
    ) %>% 
    inner_join(detection_rate)

region_da %>% write_tsv('data/diff_expression/ATAC_linked_peaks_region_da.tsv')
# region_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_region_da.tsv')

region_da_sig <- region_da %>% 
    filter(padj<0.01 & log2(detect_ratio)>1)

ggplot(region_da, aes(log2(detect_ratio), -log10(pval), alpha=feature%in%region_da_sig$feature)) +
    geom_point() +
    facet_grid(~group)


reg_peaks <- net_use$regions %>% str_split(';') %>% unlist() 
reg_peaks_ranges <- reg_peaks %>% StringToGRanges()

region_da_sig <- region_da %>% 
    filter(padj<0.01 & detect_ratio > 2, feature%in%reg_peaks)

medial_peaks <- region_da_sig %>% 
    filter(group=='medial')
dorsal_peaks <- region_da_sig %>% 
    filter(group=='dorsal')
lateral_peaks <- region_da_sig %>% 
    filter(group=='lateral')


glm_graph <- glm_graph %N>% 
    inner_join(perc_expr_df, by=c('name'='gene'))

#### Region-specific networks ####

medial_graph <- glm_graph %E>% 
    filter(regions%in%medial_peaks$feature) %N>% 
    filter(perc_expr_medial>0.1) %>% 
    filter(!node_is_isolated())

p1 <- ggraph(medial_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(filter=regions%in%medial_peaks$feature), color='darkgray', width=0.2) + 
    geom_node_point(aes(fill=medial, size=central_eig), shape=21, stroke=0.2) +
    geom_node_text(aes(label=name), size=2, repel=T) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void() + no_legend() + ggtitle('medial')

p1


dorsal_graph <- glm_graph %E>% 
    filter(regions%in%dorsal_peaks$feature) %N>% 
    filter(perc_expr_dorsal>0.1) %>% 
    filter(!node_is_isolated())

p2 <- ggraph(dorsal_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(filter=regions%in%dorsal_peaks$feature), color='darkgray', width=0.2) + 
    geom_node_point(aes(fill=dorsal, size=central_eig), shape=21, stroke=0.2) +
    geom_node_text(aes(label=name), size=2, repel=T) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void() + no_legend() + ggtitle('dorsal')

p2

lateral_graph <- glm_graph %E>% 
    filter(regions%in%lateral_peaks$feature) %N>% 
    filter(perc_expr_lateral>0.1) %>% 
    filter(!node_is_isolated())

p3 <- ggraph(lateral_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(filter=regions%in%lateral_peaks$feature), color='darkgray', width=0.2) + 
    geom_node_point(aes(fill=lateral, size=central_eig), shape=21, stroke=0.2) +
    geom_node_text(aes(label=name), size=2, repel=T) +
    scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    scale_fill_viridis(option='magma', direction = -1) +
    theme_void() + no_legend() + ggtitle('lateral')

p1 | p2 | p3
ggsave('plots/glut_region_grn.png', width=14, height=8)


#### DA between fates ####
fates <- set_names(colnames(fate_df))

fate_da <- map_dfr(fates, function(fat){
    print(fat)
    mome_atac_test$test_var <- !is.na(mome_atac_test@meta.data[fat])
    test_df <- lr_de(
        object = mome_atac_test,
        features_use = peaks_test,
        test_var = 'test_var',
        covariates = 'nFeature_ATAC',
        family = 'binomial',
        assay = 'ATAC_bin'
    )
    detection_rates <- Pando::aggregate_matrix(t(mome_atac_test[['ATAC_bin']]@counts[peaks_test, ]), mome_atac_test$test_var)
    detrate_df <- tibble(
        feature = colnames(detection_rates),
        detect_self = as.numeric(detection_rates['TRUE', ]),
        detect_other = as.numeric(detection_rates['FALSE', ])
    )
    return(inner_join(test_df, detrate_df))
}, .id='group')

fate_da <- fate_da %>% 
    mutate(
        detect_ratio=detect_self/detect_other,
        padj=p.adjust(pval, method='fdr')
    )

fate_da %>% write_tsv('data/diff_expression/ATAC_linked_peaks_fate_da.tsv')
# fate_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_fate_da.tsv')

fate_da_sig <- fate_da %>% 
    filter(padj<0.01 & log2(detect_ratio)>1)

ggplot(fate_da, aes(log2(detect_ratio), -log10(pval), alpha=feature%in%region_da_sig$feature)) +
    geom_point() +
    facet_wrap(~group, scales='free_y')

fate_da_sig <- fate_da %>% 
    filter(padj<0.01 & detect_ratio > 3, feature%in%reg_peaks)

fates_peaks <- fate_da_sig %>% 
    group_by(group) %>% 
    group_split()



#### Fate specific networks ####
fate_grnplots <- map(fates_peaks, function(peaks){
    plot_graph <- glm_graph %E>% 
        filter(regions%in%peaks$feature) %N>% 
        filter(!node_is_isolated())
    p <- ggraph(plot_graph, x=UMAP1, y=UMAP2) + 
        geom_edge_diagonal(color='darkgray', width=0.2) + 
        geom_node_point(shape=21, stroke=0.2, fill='magenta') +
        geom_node_text(aes(label=name), size=2, repel=T) +
        scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
        scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
        scale_fill_viridis(option='magma', direction = -1) +
        theme_void() + no_legend() + ggtitle(peaks$groups[1])
    return(p)
})

wrap_plots(fate_grnplots)

ggsave('plots/glut_fate_grn.png', width=14, height=8)


#### DA between fates ####
groups <- set_names(colnames(group_df))

group_da <- map_dfr(groups, function(fat){
    print(fat)
    mome_atac_test$test_var <- !is.na(mome_atac_test@meta.data[fat])
    test_df <- lr_de(
        object = mome_atac_test,
        features_use = peaks_test,
        test_var = 'test_var',
        covariates = 'nFeature_ATAC',
        family = 'binomial',
        assay = 'ATAC_bin'
    )
    detection_rates <- Pando::aggregate_matrix(t(mome_atac_test[['ATAC_bin']]@counts[peaks_test, ]), mome_atac_test$test_var)
    detrate_df <- tibble(
        feature = colnames(detection_rates),
        detect_self = as.numeric(detection_rates['TRUE', ]),
        detect_other = as.numeric(detection_rates['FALSE', ])
    )
    return(inner_join(test_df, detrate_df))
}, .id='group')


group_da <- group_da %>% 
    mutate(
        detect_ratio=detect_self/detect_other,
        padj=p.adjust(pval, method='fdr')
    )

group_da %>% write_tsv('data/diff_expression/ATAC_linked_peaks_group_da.tsv')
# group_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_group_da.tsv')

group_da_sig <- group_da %>% 
    filter(padj<0.01 & log2(detect_ratio)>1)

ggplot(group_da, aes(log2(detect_ratio), -log10(pval), alpha=feature%in%region_da_sig$feature)) +
    geom_point() +
    facet_wrap(~group, scales='free_y')

group_da_sig <- group_da %>% 
    filter(padj<0.01 & detect_ratio > 3, feature%in%reg_peaks)

groups_peaks <- group_da_sig %>% 
    group_by(group) %>% 
    group_split()


#### Group specific networks ####
group_grnplots <- map(groups_peaks, function(peaks){
    plot_graph <- glm_graph %E>% 
        filter(regions%in%peaks$feature) %N>% 
        filter(!node_is_isolated())
    p <- ggraph(plot_graph, x=UMAP1, y=UMAP2) + 
        geom_edge_diagonal(color='darkgray', width=0.2) + 
        geom_node_point(shape=21, stroke=0.2, fill='magenta') +
        geom_node_text(aes(label=name), size=2, repel=T) +
        scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
        scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
        scale_fill_viridis(option='magma', direction = -1) +
        theme_void() + no_legend() + ggtitle(peaks$group[1])
    return(p)
})

wrap_plots(group_grnplots)

ggsave('plots/glut_group_grn.png', width=14, height=8)




#### Other graph layout ####
glm_graph <- glm_graph %N>% 
    mutate(
        xpos=ifelse(name%in%net_use$tf, 0, 1),
        yranks=rank(lateral)
    )

ggraph(glm_graph, x=xpos, y=yranks) +
    geom_node_point() +
    geom_edge_link(color='darkgray', width=0.2)

group_grnplots <- map(groups_peaks, function(peaks){
    plot_graph <- glm_graph %E>% 
        filter(regions%in%peaks$feature) %N>% 
        filter(!node_is_isolated())
    p <- ggraph(plot_graph, x=xpos, y=yranks) + 
        geom_edge_diagonal(color='darkgray', width=0.2) + 
        geom_node_point(shape=21, stroke=0.2, fill='magenta') +
        geom_node_text(aes(label=name), size=2) +
        scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
        scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
        scale_fill_viridis(option='magma', direction = -1) +
        theme_void() + no_legend() + ggtitle(peaks$group[1])
    return(p)
})

wrap_plots(group_grnplots)
ggsave('plots/glut_grn_groups_alluv.pdf', width=10, height=20)



fate_grnplots <- map(fates_peaks, function(peaks){
    plot_graph <- glm_graph %E>% 
        filter(regions%in%peaks$feature) %N>% 
        filter(!node_is_isolated())
    p <- ggraph(plot_graph, x=xpos, y=yranks) + 
        geom_edge_diagonal(color='darkgray', width=0.2) + 
        geom_node_point(shape=21, stroke=0.2, fill='magenta') +
        geom_node_text(aes(label=name), size=2) +
        scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
        scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
        scale_fill_viridis(option='magma', direction = -1) +
        theme_void() + no_legend() + ggtitle(peaks$group[1])
    return(p)
})

wrap_plots(fate_grnplots)
ggsave('plots/glut_grn_fates_alluv.pdf', width=10, height=20)


region_peaks <- region_da_sig %>% 
    group_by(group) %>% 
    group_split()

region_grnplots <- map(region_peaks, function(peaks){
    plot_graph <- glm_graph %E>% 
        filter(regions%in%peaks$feature) %N>% 
        filter(!node_is_isolated())
    p <- ggraph(plot_graph, x=xpos, y=yranks) + 
        geom_edge_diagonal(color='darkgray', width=0.2) + 
        geom_node_point(shape=21, stroke=0.2, fill='magenta') +
        geom_node_text(aes(label=name), size=2) +
        scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
        scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
        scale_fill_viridis(option='magma', direction = -1) +
        theme_void() + no_legend() + ggtitle(peaks$group[1])
    return(p)
})

wrap_plots(region_grnplots)
ggsave('plots/glut_grn_regions_alluv.pdf', width=6, height=10)



#### Get stats for linked / GRN peaks ####
links_df <- mome_atac@assays$ATAC@links %>% 
    as_tibble() %>% 
    mutate(dist=width(mome_atac@assays$ATAC@links)) %>% 
    inner_join(net_use, by=c('peak'='regions')) %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
    mutate(
        class = case_when(
            dist > 5000 ~ 'distal',
            T ~ 'promotor',
        )
    )

ggplot(links_df, aes(dist)) +
    geom_histogram(color='black', fill='grey') +
    labs(x='Distance to gene', y='Count')

ggsave('plots/glut_peaks_dist_hist.png', width=5, height=4)

ggplot(links_df, aes(class)) +
    geom_bar(color='black', fill='grey') +
    scale_x_discrete(labels=c('distal (>5kb)', 'promotor (<5kb)'))
    
ggsave('plots/glut_peaks_class_bar.png', width=5, height=4)

plot_df <- links_df %>% 
    group_by(target) %>% 
    mutate(npeaks=n()) %>% 
    ungroup() %>% 
    # filter(npeaks>1) %>%
    arrange(desc(npeaks)) %>% 
    mutate(target=factor(target, levels=unique(.$target)))

ggplot(plot_df, aes(target, npeaks)) +
    geom_bar(stat='summary') +
    scale_y_continuous(expand=c(0,0)) +
    no_x_text() +
    labs(x='Target gene', y='# linked regulatory peaks')

ggsave('plots/glut_npeaks_bar.png', width=12, height=5)


#### Peak meta ####
peak_meta <- links_df %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) 
    
peaks_test <- unique(mome_atac@assays$ATAC@links$peak)

#### Aggregate peak accessibility over clusters ####
lpeak_acc <- t(GetAssayData(mome_atac, assay='ATAC', slot='counts')[peaks_test, ])
lpeak_bin <- as(lpeak_acc>0, 'dgCMatrix')
lpeak_clusters <- Pando::aggregate_matrix(
    lpeak_bin[new_trajectories$cell, ],
    groups = new_trajectories$cellclusters
)
# Summarize acc over linked genes 
lpeak_target_clusters <- Pando::aggregate_matrix(
    t(lpeak_clusters)[peak_meta$peak, ],
    groups = peak_meta$target
)

lpeak_clusters_df <- lpeak_clusters %>% 
    as_tibble(rownames='cellclusters') %>% 
    pivot_longer(!cellclusters, names_to='feature', values_to='acc') %>% 
    group_by(feature) %>% 
    filter(max(acc)>0.1)

lpeak_targets_df <- lpeak_target_clusters %>% 
    as_tibble(rownames='target') %>% 
    pivot_longer(!target, names_to='cellclusters', values_to='acc') %>% 
    group_by(target) %>% 
    filter(max(acc)>0.15) %>% 
    mutate(acc_scaled=scale01(acc))

cluster_meta <- new_trajectories %>% 
    group_by(cellclusters, fate) %>% 
    summarize(
        new_pseudotime=mean(new_pseudotime),
        normalised_pseudotime=mean(normalised_pseudotime),
        latent_time=mean(latent_time)
    )


#### Plot accessibility of linked peaks over pseudotime ####
cluster_acc_df <- lpeak_targets_df %>% 
    inner_join(cluster_meta) %>% 
    arrange(new_pseudotime) %>% 
    mutate(cellclusters=factor(cellclusters, levels=unique(.$cellclusters)))

plots <- map(unique(cluster_acc_df$fate), function(fat){
    plot_df <- cluster_acc_df %>% 
        filter(fate==fat)
    
    clust_mat <- plot_df %>% select(target, cellclusters, acc) %>% 
        pivot_wider(names_from = 'cellclusters', values_from = 'acc') %>% 
        column_to_rownames('target') %>% as.matrix()
    
    target_clust <- dist(clust_mat) %>% hclust(method='ward.D2')
    target_order <- target_clust$labels[target_clust$order]
    
    p <- ggplot(plot_df, aes(cellclusters, factor(target, levels=target_order), fill=acc_scaled)) +
        geom_tile() +
        facet_wrap(~fate, scales='free') +
        scale_fill_gradientn(colors=pals::brewer.greys(100)) +
        ggtitle(plot_df$fate[1]) +
        labs(x='pt-ordered cellclusters', y='genes', fill='Scaled\naccessibility') +
        rotate_x_text(30)
    
    return(p)
})

wrap_plots(plots)
ggsave('plots/glut_fate_cluster_target_acc_heatmap.pdf', width=20, height=45)

plot_df <- cluster_acc_df %>% 
    filter(
        fate=='glut_SUBSET_8',
        target%in%c('BCL11B', 'NOTCH2', 'KIF13A', 'SOX2')
    ) 

ggplot(plot_df, aes(cellclusters, acc, color=target, group=1)) +
    geom_point(stat='summary', fun=sum) +
    stat_summary(fun=sum, geom='line') +
    facet_wrap(~target, scales='free') +
    rotate_x_text(40) +
    ggtitle('Fate: glut_SUBSET_8') +
    labs(y='Mean detection probability of regulatory peaks')
ggsave('plots/glut_fate_glut_SUBSET_8_clusterselesct_acc_line.pdf', width=10, height=5)


plot_df <- cluster_acc_df %>% 
    filter(
        fate=='glut_SUBSET_0',
        target%in%c('FZD2', 'MCTP1', 'GAD1', 'BCL11B', 'PCSK1N', 'CNIH2')
    ) 

ggplot(plot_df, aes(cellclusters, acc, color=target, group=1)) +
    geom_point(stat='summary', fun=sum) +
    stat_summary(fun=sum, geom='line') +
    facet_wrap(~target, scales='free') +
    rotate_x_text(40) +
    ggtitle('Fate: glut_SUBSET_0') +
    labs(y='Mean detection probability of regulatory peaks')


plot_df <- cluster_acc_df %>% 
    filter(
        fate=='glut_SUBSET_0',
        target%in%c('FZD2', 'MCTP1', 'GAD1', 'BCL11B', 'PCSK1N', 'CNIH2')
    ) 

ggplot(plot_df, aes(cellclusters, acc, color=target, group=1)) +
    geom_point(stat='summary', fun=sum) +
    stat_summary(fun=sum, geom='line') +
    facet_wrap(~target, scales='free') +
    rotate_x_text(40) +
    ggtitle('Fate: glut_SUBSET_0') +
    labs(y='Mean detection probability of regulatory peaks')









#### Plot regions ####
DefaultAssay(mome_atac) <- 'ATAC'

mome_atac_cov <- DietSeurat(mome_atac, assays=c('RNA', 'ATAC'))


#### EOMES ####
region <- FindRegion(mome_atac, 'chr14q-371777101-371789106') %>% 
    resize(width = 2, fix = 'end') %>%
    Extend(upstream=50000, downstream=50000)

CoveragePlot(
    mome_atac_cov, 
    region = region, 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

mome_atac@assays$ATAC@links %>% 
    as_tibble() %>% 
    filter(peak%in%reg_peaks)


mome_atac@assays$ATAC@links[unique(queryHits(findOverlaps(mome_atac@assays$ATAC@links, mome_atac@grn@regions@ranges))),]
mome_atac@assays$ATAC@links[unique(queryHits(findOverlaps(mome_atac@assays$ATAC@links, mome_atac@assays$ATAC@ranges))),]

links_str <- GRangesToString(mome_atac@assays$ATAC@links)
grn_str <- GRangesToString(mome_atac@grn@regions@ranges)
peaks_str <- rownames(mome_atac@assays$MACS)

links_str%in%grn_str
links_str%in%peaks_str













