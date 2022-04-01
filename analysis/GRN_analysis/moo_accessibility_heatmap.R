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

glm_coefs <- coef(mome_atac, network='glm_network')
grn_net <- read_tsv('data/grn/glm_modules.tsv')

trajectory_meta <- read_csv('ref_glut_dat.csv') %>% 
    mutate(cell=str_replace_all(newcellnames, '_', '-')) %>% 
    filter(cell%in%colnames(mome_atac))


fate_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_fate_da.tsv')
region_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_region_da.tsv')
group_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_group_da.tsv')


#### Peak meta ####
peak_meta <- Links(mome_atac) %>% 
    as_tibble() %>% 
    filter(!str_detect(gene, '^(AMEX|LOC)')) 
peaks_test <- unique(mome_atac@assays$ATAC@links$peak)


#### Get peak stats ####
peak_counts <- GetAssayData(mome_atac, assay='ATAC', slot='counts')
# peak_counts <- peak_counts[peaks_test, ]
idx_keep <- sparseMatrixStats::rowSums2(peak_counts) > 0
peak_counts <- peak_counts[idx_keep, , drop = FALSE]

peak_ranges <- StringToGRanges(rownames(peak_counts))

chromvar_obj <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = peak_counts),
    rowRanges = peak_ranges
)

chromvar_obj <- chromVAR::addGCBias(
    object = chromvar_obj,
    genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M
)

bg <- chromVAR::getBackgroundPeaks(
    object = chromvar_obj
)



#### Get enrichments for gene target peaks ####
target_peaks <- glm_coefs %>% 
    filter(padj<0.05) %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
    group_by(target) %>% 
    group_split() %>% 
    {names(.) <- map(., ~.x$target[1]);.} %>% 
    map(~unique(.x$region)) 

module_matrix <- target_peaks %>%
    map(function(x) {
        out <- rep(1, length(x))
        names(out) <- x
        return(out)
    }) %>% bind_rows() %>% as.matrix()
module_matrix[is.na(module_matrix)] <- 0
module_matrix <- t(Matrix(module_matrix, sparse=T))
colnames(module_matrix) <- names(target_peaks)

module_matrix_use <- Matrix(data=0, nrow=nrow(peak_counts), ncol=ncol(module_matrix))
colnames(module_matrix_use) <- colnames(module_matrix)
rownames(module_matrix_use) <- rownames(peak_counts)
module_matrix_use[rownames(module_matrix), ] <- module_matrix

dev <- chromVAR::computeDeviations(
    object = chromvar_obj,
    annotations = module_matrix_use,
    background_peaks = bg
)
chromvar_z <- SummarizedExperiment::assays(dev)[[2]]
rownames(chromvar_z) <- colnames(module_matrix)

chromvar_z %>% write_rds('data/grn/target_peaks_chromvar.rds')


#### Get enrichments for tf regulatory peaks ####
target_peaks <- glm_coefs %>% 
    filter(padj<0.05) %>% 
    filter(!str_detect(target, '^(AMEX|LOC)')) %>% 
    group_by(tf) %>% 
    group_split() %>% 
    {names(.) <- map(., ~.x$tf[1]);.} %>% 
    map(~unique(.x$region)) 

module_matrix <- target_peaks %>%
    map(function(x) {
        out <- rep(1, length(x))
        names(out) <- x
        return(out)
    }) %>% bind_rows() %>% as.matrix()
module_matrix[is.na(module_matrix)] <- 0
module_matrix <- t(Matrix(module_matrix, sparse=T))
colnames(module_matrix) <- names(target_peaks)

module_matrix_use <- Matrix(data=0, nrow=nrow(peak_counts), ncol=ncol(module_matrix))
colnames(module_matrix_use) <- colnames(module_matrix)
rownames(module_matrix_use) <- rownames(peak_counts)
module_matrix_use[rownames(module_matrix), ] <- module_matrix

dev <- chromVAR::computeDeviations(
    object = chromvar_obj,
    annotations = module_matrix_use,
    background_peaks = bg
)
chromvar_z <- SummarizedExperiment::assays(dev)[[2]]
rownames(chromvar_z) <- colnames(module_matrix)

chromvar_z %>% write_rds('data/grn/tf_peaks_chromvar.rds')


#### Get enrichments for gene linked peaks ####
target_peaks <- peak_meta %>% 
    filter(!str_detect(gene, '^(AMEX|LOC)')) %>% 
    group_by(gene) %>% 
    group_split() %>% 
    {names(.) <- map(., ~.x$gene[1]);.} %>% 
    map(~unique(.x$peak)) 

module_matrix <- target_peaks %>%
    map(function(x) {
        out <- rep(1, length(x))
        names(out) <- x
        return(out)
    }) %>% bind_rows() %>% as.matrix()
module_matrix[is.na(module_matrix)] <- 0
module_matrix <- t(Matrix(module_matrix, sparse=T))
colnames(module_matrix) <- names(target_peaks)

module_matrix_use <- Matrix(data=0, nrow=nrow(peak_counts), ncol=ncol(module_matrix))
colnames(module_matrix_use) <- colnames(module_matrix)
rownames(module_matrix_use) <- rownames(peak_counts)
module_matrix_use[rownames(module_matrix), ] <- module_matrix

dev <- chromVAR::computeDeviations(
    object = chromvar_obj,
    annotations = module_matrix_use,
    background_peaks = bg
)
chromvar_z <- SummarizedExperiment::assays(dev)[[2]]
rownames(chromvar_z) <- colnames(module_matrix)

chromvar_z %>% write_rds('data/grn/linked_peaks_chromvar.rds')



#### Get pseudotime bins ####
trajectory_list <- trajectory_meta %>% 
    filter(cell%in%colnames(peak_counts)) %>% 
    group_by(group) %>% 
    group_split() %>% 
    map(function(x){
        x %>% 
            mutate(pt_ranks=rank(new_pseudotime)) %>% 
            mutate(pt_norm=pt_ranks/max(pt_ranks)) %>% 
            mutate(pt_bins=cut(pt_norm, breaks=20, labels=F)) %>% 
            return()
    })


#### Aggregate peak accessibility over bins ####
linked_peaks <- peak_meta %>% 
    filter(!str_detect(gene, '^(AMEX|LOC)')) 

lpeak_acc <- t(GetAssayData(mome_atac, assay='ATAC', slot='counts')[peaks_test, ])
lpeak_bin <- as(lpeak_acc>0, 'dgCMatrix')
lpeak_clusters <- map(
    trajectory_list,
    function(traj){
        Pando::aggregate_matrix(
            lpeak_bin[traj$cell, ],
            groups = traj$pt_bins
        )
    }
)
names(lpeak_clusters) <- map_chr(trajectory_list, ~.x$group[1])


lpeak_clusters_df <- lpeak_clusters %>% 
    map_dfr(function(x){
        x %>% as_tibble(rownames='pt_bins') %>% 
            pivot_longer(!pt_bins, names_to='feature', values_to='acc') %>% 
            return()
    }, .id='group')


da_peaks <- group_da %>% 
    filter(padj<0.05, detect_ratio>2, coef>0, detect_self>0.02) %>% 
    group_by(group) %>% 
    top_n(30, detect_ratio) %>% 
    rename('da_group'=group) %>% 
    inner_join(lpeak_clusters_df) %>% 
    inner_join(linked_peaks, by=c('feature'='peak')) %>% 
    group_by(feature) %>% 
    mutate(
        acc01=scale01(acc),
        acc_scale=zscale(acc),
        pt_bins=factor(pt_bins, levels=sort(unique(as.numeric(pt_bins))))
    )

da_peaks_all <- lpeak_clusters %>% map(~t(.x)[unique(da_peaks$feature), ])
da_peaks_all <- do.call(cbind, da_peaks_all)

da_peaks_all_gene <- da_peaks_all
rownames(da_peaks_all_gene) <- linked_peaks$gene[match(rownames(da_peaks_all), linked_peaks$peak)]

peaks_order <- da_peaks_all %>% dist() %>% hclust(method='ward.D2') %>% {.$labels[.$order]}
gene_order <- da_peaks_all_gene %>% dist() %>% hclust(method='ward.D2') %>% {.$labels[.$order]}


plot_df <- da_peaks %>% 
    mutate(
        gene=factor(gene, levels=unique(gene_order)),
        feature=factor(feature, levels=unique(peaks_order))
    ) %>% 
    filter(!is.na(gene))

ggplot(plot_df, aes(pt_bins, gene, fill=acc_scale)) +
    geom_tile() +
    facet_grid(da_group~group, scales='free', space='free') +
    scale_fill_gradientn(colors=pals::inferno(100)) +
    no_x_text()

ggsave('plots/glut_trajectories_da_peaks_lgenes_heatmap.pdf', width=15, height=15)
ggsave('plots/glut_trajectories_da_peaks_lgenes_heatmap.png', width=15, height=15)


ggplot(plot_df, aes(pt_bins, feature, fill=acc_scale)) +
    geom_tile() +
    facet_grid(da_group~group, scales='free', space='free') +
    scale_fill_gradientn(colors=pals::inferno(100))





#### Summarize to target peaks ####
target_peaks <- glm_coefs %>% 
    filter(padj<0.05) %>% 
    filter(!str_detect(target, '^(AMEX|LOC)'))

linked_peaks <- peak_meta %>% 
    filter(!str_detect(gene, '^(AMEX|LOC)')) 

target_clusters <- map(lpeak_clusters, function(x){
    t(Pando::aggregate_matrix(t(x)[target_peaks$region, ], target_peaks$target))
})

tf_clusters <- map(lpeak_clusters, function(x){
    t(Pando::aggregate_matrix(t(x)[target_peaks$region, ], target_peaks$tf))
})

linked_clusters <- map(lpeak_clusters, function(x){
    t(Pando::aggregate_matrix(t(x)[linked_peaks$peak, ], linked_peaks$gene))
})


target_order <- target_clusters[[1]] %>% t() %>% dist() %>% hclust() %>% {.$labels[.$order]}

lpeak_clusters_df <- target_clusters[[1]] %>% 
    as_tibble(rownames='pt_bins') %>% 
    pivot_longer(!pt_bins, names_to='feature', values_to='acc') %>% 
    group_by(feature) %>% 
    # filter(max(acc)>0.1) %>% 
    mutate(
        acc01=scale01(acc),
        acc_scale=zscale(acc),
        pt_bins=factor(pt_bins, levels=sort(as.numeric(pt_bins))),
        feature=factor(feature, levels=target_order)
    )

ggplot(lpeak_clusters_df, aes(pt_bins, feature, fill=acc01)) +
    geom_tile() +
    scale_fill_gradientn(colors=pals::magma(100))


tf_order <- tf_clusters[[1]] %>% t() %>% dist() %>% hclust() %>% {.$labels[.$order]}

lpeak_clusters_df <- tf_clusters[[1]] %>% 
    as_tibble(rownames='pt_bins') %>% 
    pivot_longer(!pt_bins, names_to='feature', values_to='acc') %>% 
    group_by(feature) %>% 
    # filter(max(acc)>0.1) %>% 
    mutate(
        acc01=scale01(acc),
        acc_scale=zscale(acc),
        pt_bins=factor(pt_bins, levels=sort(as.numeric(pt_bins))),
        feature=factor(feature, levels=tf_order)
    )

ggplot(lpeak_clusters_df, aes(pt_bins, feature, fill=acc)) +
    geom_tile() +
    scale_fill_gradientn(colors=pals::magma(100))


gene_order <- linked_clusters[[4]] %>% t() %>% dist() %>% hclust() %>% {.$labels[.$order]}

lpeak_clusters_df <- linked_clusters[[4]] %>% 
    as_tibble(rownames='pt_bins') %>% 
    pivot_longer(!pt_bins, names_to='feature', values_to='acc') %>% 
    group_by(feature) %>% 
    filter(max(acc)>0.1) %>%
    mutate(
        acc01=scale01(acc),
        acc_scale=zscale(acc),
        pt_bins=factor(pt_bins, levels=sort(as.numeric(pt_bins))),
        feature=factor(feature, levels=gene_order)
    )

ggplot(lpeak_clusters_df, aes(pt_bins, feature, fill=acc)) +
    geom_tile() +
    scale_fill_gradientn(colors=pals::viridis(100))









