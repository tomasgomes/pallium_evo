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
mome_atac <- read_rds('data/muo_glut_v1.4annot_links10m_srt.rds')
rna_expr <- t(Seurat::GetAssayData(mome_atac, assay='RNA', slot='data'))

glm_coefs <- coef(mome_atac, network='glm_network')
grn_net <- read_tsv('data/grn/glm_modules.tsv')

trajectory_meta <- read_csv('ref_glut_dat.csv') %>% 
    mutate(cell=str_replace_all(newcellnames, '_', '-')) %>% 
    filter(cell%in%colnames(mome_atac))

#### Prep DA ####
mome_atac@active.assay <- 'ATAC'
mome_atac_test <- DietSeurat(mome_atac, assays = c('ATAC'))
mome_atac_peaks_bin <- as(mome_atac@assays$ATAC@counts>0, 'dgCMatrix')
mome_atac_test[['ATAC_bin']] <- CreateAssayObject(mome_atac_peaks_bin)
peaks_test <- unique(mome_atac@assays$ATAC@links$peak)

mome_terminal <- mome_atac_test[, mome_atac_test$cellclusters %in% unique(trajectory_meta$fate)]
mome_neurons <- mome_atac_test[, mome_atac_test$classes=='neuronal']

peak_meta <- Links(mome_atac) %>% 
    as_tibble() %>% 
    filter(!str_detect(gene, '^(AMEX|LOC)')) %>% 
    filter(!str_detect(gene, '\\.\\.'))

#### DA between terminal states ####
terminal_states <- purrr::set_names(unique(trajectory_meta$fate))
terminal_da <- map_dfr(terminal_states, function(x){
    print(x)
    mome_terminal$test_var <- mome_terminal$cellclusters == x
    test_df <- lr_de(
        object = mome_terminal,
        features_use = peaks_test,
        test_var = 'test_var',
        covariates = 'nFeature_ATAC',
        family = 'binomial',
        assay = 'ATAC_bin'
    )
    detection_rates <- Pando::aggregate_matrix(
        t(mome_terminal[['ATAC_bin']]@counts[peaks_test, ]), 
        mome_terminal$test_var
    )
    detrate_df <- tibble(
        feature = colnames(detection_rates),
        detect_self = as.numeric(detection_rates['TRUE', ]),
        detect_other = as.numeric(detection_rates['FALSE', ])
    )
    return(inner_join(test_df, detrate_df))
}, .id='group')

mome_terminal <- aggregate_assay(mome_terminal, group_name='cellclusters', assay='ATAC_bin')
cluster_detection <- mome_terminal@assays$ATAC_bin@misc$summary$cellclusters
mome_terminal <- aggregate_assay(mome_terminal, group_name='region', assay='ATAC_bin')
region_detection <- mome_terminal@assays$ATAC_bin@misc$summary$region

detection_rate <- tibble(
    feature = colnames(cluster_detection),
    max_cluster_detection = colMaxs(cluster_detection),
    max_region_detection = colMaxs(region_detection)
)

terminal_annot_da <- terminal_da %>% 
    mutate(
        detect_ratio=detect_self/detect_other,
        padj=p.adjust(pval, method='fdr')
    ) %>% 
    inner_join(detection_rate) %>% 
    inner_join(peak_meta, by=c('feature'='peak')) %>% 
    mutate(
        distal = width>2000
    )

terminal_annot_da %>% write_tsv('data/diff_expression/ATAC_linked_peaks_terminal_da.tsv')


#### DA between neurons ####
neuron_da <- map_dfr(terminal_states, function(x){
    print(x)
    mome_neurons$test_var <- !is.na(mome_neurons@meta.data[x])
    test_df <- lr_de(
        object = mome_neurons,
        features_use = peaks_test,
        test_var = 'test_var',
        covariates = 'nFeature_ATAC',
        family = 'binomial',
        assay = 'ATAC_bin'
    )
    detection_rates <- Pando::aggregate_matrix(
        t(mome_neurons[['ATAC_bin']]@counts[peaks_test, ]), 
        mome_neurons$test_var
    )
    detrate_df <- tibble(
        feature = colnames(detection_rates),
        detect_self = as.numeric(detection_rates['TRUE', ]),
        detect_other = as.numeric(detection_rates['FALSE', ])
    )
    return(inner_join(test_df, detrate_df))
}, .id='group')

mome_neurons <- aggregate_assay(mome_neurons, group_name='cellclusters', assay='ATAC_bin')
cluster_detection <- mome_neurons@assays$ATAC_bin@misc$summary$cellclusters
mome_neurons <- aggregate_assay(mome_neurons, group_name='region', assay='ATAC_bin')
region_detection <- mome_neurons@assays$ATAC_bin@misc$summary$region

detection_rate <- tibble(
    feature = colnames(cluster_detection),
    max_cluster_detection = colMaxs(cluster_detection),
    max_region_detection = colMaxs(region_detection)
)

neuron_annot_da <- neuron_da %>% 
    mutate(
        detect_ratio=detect_self/detect_other,
        padj=p.adjust(pval, method='fdr')
    ) %>% 
    inner_join(detection_rate) %>%
    inner_join(peak_meta, by=c('feature'='peak')) %>% 
    mutate(
        distal = width>2000
    )

neuron_annot_da %>% write_tsv('data/diff_expression/ATAC_linked_peaks_neuron_da.tsv')


#### Get pseudotime bins ####
trajectory_list <- trajectory_meta %>% 
    filter(cell%in%colnames(peak_counts)) %>% 
    group_by(fate) %>% 
    group_split() %>% 
    map(function(x){
        x %>% 
            mutate(pt_ranks=rank(new_pseudotime)) %>% 
            mutate(pt_norm=pt_ranks/max(pt_ranks)) %>% 
            mutate(pt_bins=cut(pt_norm, breaks=20, labels=F)) %>% 
            return()
    })

fate_meta <- trajectory_meta %>% 
    select('da_group'=fate, 'traj_group'=group) %>% 
    distinct()


#### Heatmap with DA peaks ####
lpeak_acc <- t(GetAssayData(mome_atac, assay='ATAC', slot='counts')[peaks_test, ])
lpeak_bin <- as(lpeak_acc>0, 'dgCMatrix')
lpeak_clusters <- map(
    trajectory_list,
    function(traj){
        Pando::aggregate_matrix(
            lpeak_bin[traj$cell, ],
            groups = mome_atac$subclasses[traj$cell]
        )
    }
)
names(lpeak_clusters) <- map_chr(trajectory_list, ~.x$fate[1])

lpeak_clusters_df <- lpeak_clusters %>% 
    map_dfr(function(x){
        x %>% as_tibble(rownames='class') %>% 
            pivot_longer(!class, names_to='feature', values_to='acc') %>% 
            return()
    }, .id='group')


neuron_order <- c('glut_SUBSET_1', 'glut_SUBSET_11', 'glut_SUBSET_0', 'glut_SUBSET_7',
                  'glut_SUBSET_6', 'glut_SUBSET_8', 'glut_SUBSET_10', 'glut_SUBSET_22', 
                  'glut_SUBSET_2', 'glut_SUBSET_9')

da_peaks <- terminal_annot_da %>% 
    filter(padj<0.05, detect_ratio>1.5, coef>0, detect_self>0.02) %>% 
    group_by(group, distal) %>% 
    top_n(30, detect_ratio) %>% 
    rename('da_group'=group) %>% 
    inner_join(lpeak_clusters_df) %>% 
    inner_join(fate_meta) %>% 
    group_by(feature) %>% 
    mutate(
        acc01=scale01(acc),
        acc_scale=zscale(acc),
        class=factor(class, levels=c('Ependymal', 'NPC', 'Glutamatergic')),
        group=factor(group, levels=neuron_order),
        da_group=factor(da_group, levels=neuron_order),
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

ggplot(plot_df, aes(class, gene, fill=acc_scale)) +
    geom_tile() +
    facet_grid(da_group~group, scales='free', space='free') +
    scale_fill_gradientn(colors=pals::inferno(100)) +
    rotate_x_text(40) +
    ggtitle('All linked DA peaks grouped by class')

ggsave('plots/heatmaps/term_da_class_heatmap.pdf', width=8, height=25)
ggsave('plots/heatmaps/term_da_class_heatmap.png', width=8, height=25)

ggplot(filter(plot_df, distal), aes(class, gene, fill=acc_scale)) +
    geom_tile() +
    facet_grid(da_group~group, scales='free', space='free') +
    scale_fill_gradientn(colors=pals::inferno(100)) +
    rotate_x_text(40) +
    ggtitle('All distal linked DA peaks grouped by class')

ggsave('plots/heatmaps/term_da_class_distal_heatmap.pdf', width=8, height=25)
ggsave('plots/heatmaps/term_da_class_distal_heatmap.png', width=8, height=25)

ggplot(filter(plot_df, !distal), aes(class, gene, fill=acc_scale)) +
    geom_tile() +
    facet_grid(da_group~group, scales='free', space='free') +
    scale_fill_gradientn(colors=pals::inferno(100)) +
    rotate_x_text(40) +
    ggtitle('All proximal linked DA peaks grouped by class')

ggsave('plots/heatmaps/term_da_class_proximal_heatmap.pdf', width=8, height=10)
ggsave('plots/heatmaps/term_da_class_proximal_heatmap.png', width=8, height=10)


ggplot(filter(plot_df, !distal & class=='Glutamatergic'), aes(group, gene, fill=acc_scale)) +
    geom_tile() +
    facet_grid(da_group~., scales='free', space='free') +
    scale_fill_gradientn(colors=pals::inferno(100)) +
    rotate_x_text(40) +
    ggtitle('All proximal linked DA peaks in neurons')

ggsave('plots/heatmaps/term_da_neuron_proximal_heatmap.pdf', width=6, height=10)
ggsave('plots/heatmaps/term_da_neuron_proximal_heatmap.png', width=6, height=10)


da_peaks %>% write_rds('data/diff_expression/ATAC_linked_peaks_terminal_da_annot.tsv')

