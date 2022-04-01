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
registerDoParallel(48) 

setwd('~/projects/axolotl/')



#### Read stuff ####
rna_all <- read_rds('/links/groups/treutlein/USERS/tomasgomes/projects/pallium_evo/data/expression/axolotl_reclust/all_nuclei_clustered_highlevel_anno.RDS')

newcellnames = colnames(rna_all)
newcellnames = gsub("-1_1", "_a1_1", newcellnames)
newcellnames = gsub("-1_2", "_a1_2", newcellnames)
newcellnames = gsub("-1_3", "_a3_1", newcellnames)
newcellnames = gsub("-1_4", "_a3_2", newcellnames)

newcellnames = gsub("-1_1_5", "_D1", newcellnames)
newcellnames = gsub("-1_2_5", "_D2", newcellnames)
newcellnames = gsub("-1_3_5", "_L1", newcellnames)
newcellnames = gsub("-1_4_5", "_L2", newcellnames)
newcellnames = gsub("-1_5_5", "_M1", newcellnames)
newcellnames = gsub("-1_6_5", "_M2", newcellnames)

rna_all <- RenameCells(rna_all, new.names = newcellnames)
rna_expr <- t(Seurat::GetAssayData(rna_all, assay='RNA', slot='data'))

mome_atac <- read_rds('data/muo_glut_v1.4annot_links10m_srt.rds')
# rna_expr <- t(Seurat::GetAssayData(mome_atac, assay='RNA', slot='data'))
new_trajectories <- read_csv('ref_glut_dat.csv') %>% 
    mutate(cell=str_replace_all(newcellnames, '_', '-')) %>% 
    filter(cell%in%colnames(mome_atac))

fate_df <- new_trajectories %>% 
    select(newcellnames, fate, new_pseudotime) %>% distinct() %>% 
    pivot_wider(names_from = 'fate', values_from = 'new_pseudotime') %>% 
    column_to_rownames('newcellnames')

group_df <- new_trajectories %>% 
    select(newcellnames, group) %>% distinct() %>% 
    mutate(group=paste0('group_', group)) %>% 
    mutate(
        val=group,
    ) %>% 
    pivot_wider(names_from = 'group', values_from = 'val') %>% 
    column_to_rownames('newcellnames')

rna_all <- AddMetaData(rna_all, fate_df)
rna_all <- AddMetaData(rna_all, group_df)

tfs <- read_tsv('~/resources/DB/animal_tfdb/tf_human.tsv')

glm_coefs <- coef(mome_atac, network='glm_network')
bagridge_coefs <- coef(mome_atac, network='bagging_ridge_network')

grn_net <- glm_coefs %>% 
    filter(padj<0.05) %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
    group_by(tf, target) %>% 
    filter(pval==min(pval)) %>% 
    summarize(
        estimate=mean(estimate),
        pval=min(pval),
        padj=min(padj),
        regions=region[1]
    )

grn_net %>% write_tsv('data/grn/glm_modules.tsv')

fate_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_fate_da.tsv')
region_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_region_da.tsv')
group_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_group_da.tsv')



#### Group average expression ####
traj_groups <- new_trajectories %>% 
    group_by(group) %>% 
    group_split()
group_names <- traj_groups %>% map_chr(function(x) x$group[1])

group_expr <- Pando::map_par(traj_groups, function(x){
    colMeans2(rna_expr[x$newcellnames, ])
}, parallel=T)
group_expr <- do.call(rbind, group_expr) %>% Matrix::Matrix(sparse=T)
rownames(group_expr) <- group_names
colnames(group_expr) <- colnames(rna_expr)

group_expr_df <- group_expr %>% 
    as_tibble(rownames='group') %>% 
    pivot_longer(!group, values_to = 'expr')

group_prcexpr <- Pando::map_par(traj_groups, function(x){
    colMeans2(rna_expr[x$newcellnames, ]>0)
}, parallel=T)
group_prcexpr <- do.call(rbind, group_prcexpr) %>% Matrix::Matrix(sparse=T)
rownames(group_prcexpr) <- group_names
colnames(group_prcexpr) <- colnames(rna_expr)

group_prcexpr_df <- group_prcexpr %>% 
    as_tibble(rownames='group') %>% 
    pivot_longer(!group, values_to = 'prcex')


#### Net with only TFs ####
#### Get avg expression for genes ####
region_summary <- Pando::aggregate_matrix(rna_expr[, union(grn_net$tf, grn_net$target)], groups=mome_atac$pred_regions_all)
subclass_summary <- Pando::aggregate_matrix(rna_expr[, union(grn_net$tf, grn_net$target)], groups=mome_atac$subclasses)

region_summary_df <- region_summary %>% t() %>% 
    as_tibble(rownames='gene') 

subclass_summary_df <- subclass_summary %>% t() %>% 
    as_tibble(rownames='gene') 

gene_scores <- inner_join(region_summary_df, subclass_summary_df)


### Get coex and umap ####
gene_cor <- Pando::sparse_cor(rna_expr[, union(grn_net$tf, grn_net$target)])

reg_mat <- grn_net %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
    filter(target%in%tfs$symbol) %>%
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


tf_graph <- as_tbl_graph(grn_net) %>% 
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
    filter(padj<0.05) %N>% 
    filter(!node_is_isolated())



ggraph(tf_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(alpha=-log10(padj), color=factor(sign(estimate))), width=0.5) + 
    geom_node_point(aes(size=outdegree), shape=21, color='black', fill='grey') +
    geom_node_text(aes(label=name), size=5/ggplot2::.pt, repel=T) +
    scale_edge_color_manual(values=c('#f5b7b1', '#7dcea0')) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    theme_void() 
ggsave('plots/tf_grn_umap.png', width=6, height=6)
ggsave('plots/tf_grn_umap.pdf', width=6, height=6)


ggraph(tf_graph, x=UMAP1, y=UMAP2) + 
    geom_edge_diagonal(aes(alpha=-log10(padj), color=factor(sign(estimate))), width=0.5) + 
    geom_node_point(aes(size=outdegree), shape=21, color='black', fill='grey') +
    # geom_node_text(aes(label=name), size=5/ggplot2::.pt, repel=T) +
    scale_edge_color_manual(values=c('#f5b7b1', '#7dcea0')) +
    scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
    theme_void() 
ggsave('plots/tf_grn_unlabelled_umap.png', width=6, height=6)
ggsave('plots/tf_grn_unlabelled_umap.pdf', width=6, height=6)


group_da_sig <- group_da %>% 
    filter(padj<0.05)

groups_peaks <- group_da_sig %>% 
    group_by(group) %>% 
    group_split()

group_grnplots <- map(groups_peaks, function(peaks){
    group_name <- peaks$group[1] %>% str_remove('group_')
    gex <- group_expr_df %>% filter(group==group_name)
    gprc <- group_prcexpr_df %>% filter(group==group_name)
    plot_graph <- tf_graph %E>% 
        filter(regions%in%peaks$feature) %N>% 
        inner_join(gex) %>% 
        inner_join(gprc) %>% 
        filter(prcex>0.1) %>% 
        arrange(expr)
    p <- ggraph(plot_graph, x=UMAP1, y=UMAP2) + 
        geom_edge_diagonal(color='darkgray', width=0.2) + 
        geom_node_point(aes(fill=expr, size=prcex), shape=21, stroke=0.2) +
        geom_node_text(aes(label=name), size=5/ggplot2::.pt, repel=T) +
        scale_size_continuous(range=c(0.5, 3)) +
        scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
        scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
        scale_fill_viridis(option='magma') +
        theme_void() + ggtitle(peaks$group[1])
    return(p)
})

wrap_plots(group_grnplots) + plot_layout(guides='collect')
ggsave('plots/tf_group_grn_umap.png', width=8, height=4)
ggsave('plots/tf_group_grn_umap.pdf', width=8, height=4)

group_grnplots <- map(groups_peaks, function(peaks){
    group_name <- peaks$group[1] %>% str_remove('group_')
    gex <- group_expr_df %>% filter(group==group_name)
    gprc <- group_prcexpr_df %>% filter(group==group_name)
    plot_graph <- tf_graph %E>% 
        filter(regions%in%peaks$feature) %N>% 
        inner_join(gex) %>% 
        inner_join(gprc) %>% 
        filter(prcex>0.1)
    p <- ggraph(plot_graph, x=UMAP1, y=UMAP2) + 
        geom_edge_diagonal(color='darkgray', width=0.2) + 
        geom_node_point(aes(fill=expr, size=prcex), shape=21, stroke=0.2) +
        scale_size_continuous(range=c(0.5, 3)) +
        scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
        scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
        scale_fill_viridis(option='magma') +
        theme_void() + ggtitle(peaks$group[1])
    return(p)
})

wrap_plots(group_grnplots) + plot_layout(guides='collect')
ggsave('plots/tf_group_grn_unlabelled_umap.png', width=8, height=4)
ggsave('plots/tf_group_grn_unlabelled_umap.pdf', width=8, height=4)






region_da_sig <- region_da %>% 
    filter(padj<0.05)

regions_peaks <- region_da_sig %>% 
    group_by(group) %>% 
    group_split()

region_grnplots <- map(regions_peaks, function(peaks){
    group_name <- peaks$group[1] 
    plot_graph <- tf_graph %E>% 
        filter(regions%in%peaks$feature) %N>% 
        filter(!node_is_isolated())
    p <- ggraph(plot_graph, x=UMAP1, y=UMAP2) + 
        geom_edge_diagonal(color='darkgray', width=0.2) + 
        geom_node_point(aes_string(fill=group_name, size=group_name), shape=21, stroke=0.2) +
        geom_node_text(aes(label=name), size=5/ggplot2::.pt, repel=T) +
        scale_size_continuous(range=c(0.5, 8)) +
        scale_edge_color_gradientn(colors=rev(rdbu()), limits=c(-6,6)) +
        scale_edge_alpha_continuous(range=c(0.01,0.8), limits=c(2,20)) +
        scale_fill_viridis(option='magma') +
        theme_void() + ggtitle(peaks$group[1])
    return(p)
})

wrap_plots(region_grnplots) + plot_layout(guides='collect')
ggsave('plots/tf_group_grn_umap.png', width=8, height=4)
ggsave('plots/tf_group_grn_umap.pdf', width=8, height=4)

