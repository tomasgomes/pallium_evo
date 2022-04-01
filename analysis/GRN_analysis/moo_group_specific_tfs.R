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
mome_atac <- read_rds('data/muo_glut_v1.4annot_links10m_srt.rds')
new_trajectories <- read_csv('ref_glut_dat.csv') %>% 
    mutate(cell=str_replace_all(newcellnames, '_', '-')) %>% 
    filter(cell%in%colnames(mome_atac))

tfs <- read_tsv('~/resources/DB/animal_tfdb/tf_human.tsv')

glm_coefs <- coef(mome_atac, network='glm_network')
bagridge_coefs <- coef(mome_atac, network='bagging_ridge_network')
grn_net <- read_tsv('data/grn/glm_modules.tsv')

grn_graph <- as_tbl_graph(grn_net)

fate_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_fate_da.tsv')
region_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_region_da.tsv')
group_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_group_da.tsv')


#### Get group-specific GRNs ####
group_da_sig <- group_da %>% 
    filter(padj<0.05, coef>0, detect_ratio>1.5) 

group_da_peaks <- group_da_sig %>% 
    group_by(group) %>% group_split()

group_spec_peaks <- group_da_sig %>% 
    group_by(feature) %>% 
    filter(detect_ratio==max(detect_ratio)) %>% 
    distinct(group, feature, detect_ratio)

names(group_da_peaks) <- map_chr(group_da_peaks, function(x) x$group[1])
da_peaks_all <- group_da_peaks %>% map(~.x$feature) %>% purrr::reduce(union)

group_nets <- map(group_da_peaks, function(peaks){
    grn_net %>% 
        filter(regions%in%peaks$feature) %>% 
        group_by(tf) %>% 
        mutate(ngenes=n(), group=peaks$group[1]) %>% 
        group_by(target) %>% 
        mutate(ntfs=n()) %>% 
        return()
})

group_nets %>% write_rds('data/grn/group_grns.rds')



#### Get TFs with most group-specific connections ####

grn_net <- grn_net %>% 
    group_by(tf) %>% 
    mutate(
        n_common = sum(regions%in%da_peaks_all),
        n_group_111 = sum(regions%in%group_da_peaks$group_111$feature),
        n_group_8620_ep = sum(regions%in%group_da_peaks$group_8620_ep$feature),
        n_group_eomes = sum(regions%in%group_da_peaks$group_eomes$feature),
        n_group_hippocampus = sum(regions%in%group_da_peaks$group_hippocampus$feature),
        n_group_lc = sum(regions%in%group_da_peaks$group_lc$feature)
    ) %>% 
    mutate(
        frac_common = n_common / length(unique(regions)),
        frac_group_111 = n_group_111 / length(unique(regions)),
        frac_group_8620_ep = n_group_8620_ep / length(unique(regions)),
        frac_group_eomes = n_group_eomes / length(unique(regions)),
        frac_group_hippocampus = n_group_hippocampus / length(unique(regions)),
        frac_group_lc = n_group_lc / length(unique(regions))
    )


grn_net_annot <- grn_net %>% 
    left_join(group_spec_peaks, by=c('regions'='feature')) %>% 
    group_by(tf) %>% 
    mutate(count=n()) %>% 
    ungroup()

plot_df <- grn_net_annot %>% 
    arrange(desc(count)) %>% 
    mutate(
        tf=factor(tf, levels=unique(.$tf)),
        group=ifelse(is.na(group), 'commmon', group)
    ) %>% 
    filter(count>10)

p1 <- ggplot(plot_df, aes(tf, fill=group)) +
    geom_bar() +
    rotate_x_text(40)


plot_df <- grn_net_annot %>% 
    select(tf, n_common:frac_group_lc) %>% distinct() %>% 
    arrange(desc(n_group_111)) %>% 
    filter(row_number()<=25) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf))) %>% 
    filter(n_group_111>5)

p2 <- ggplot(plot_df, aes(tf, n_group_111, fill=frac_group_111)) +
    geom_bar(stat='identity') +
    rotate_x_text(40) +
    scale_fill_gradientn(colors=pals::brewer.rdpu(100)) +
    labs(fill='frac_specific') + ggtitle('group_111')


plot_df <- grn_net_annot %>% 
    select(tf, n_common:frac_group_lc) %>% distinct() %>% 
    arrange(desc(n_group_8620_ep)) %>% 
    filter(row_number()<=25) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf))) %>% 
    filter(n_group_8620_ep>5)

p3 <- ggplot(plot_df, aes(tf, n_group_8620_ep, fill=frac_group_8620_ep)) +
    geom_bar(stat='identity') +
    rotate_x_text(40) +
    scale_fill_gradientn(colors=pals::brewer.rdpu(100)) +
    labs(fill='frac_specific') + ggtitle('group_8620_ep')


plot_df <- grn_net_annot %>% 
    select(tf, n_common:frac_group_lc) %>% distinct() %>% 
    arrange(desc(n_group_eomes)) %>% 
    filter(row_number()<=25) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf))) %>% 
    filter(n_group_eomes>5)

p4 <- ggplot(plot_df, aes(tf, n_group_eomes, fill=frac_group_eomes)) +
    geom_bar(stat='identity') +
    rotate_x_text(40) +
    scale_fill_gradientn(colors=pals::brewer.rdpu(100)) +
    labs(fill='frac_specific') + ggtitle('group_eomes')


plot_df <- grn_net_annot %>% 
    select(tf, n_common:frac_group_lc) %>% distinct() %>% 
    arrange(desc(n_group_hippocampus)) %>% 
    filter(row_number()<=25) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf))) %>% 
    filter(n_group_hippocampus>5)

p6 <- ggplot(plot_df, aes(tf, n_group_hippocampus, fill=frac_group_hippocampus)) +
    geom_bar(stat='identity') +
    rotate_x_text(40) +
    scale_fill_gradientn(colors=pals::brewer.rdpu(100)) +
    labs(fill='frac_specific') + ggtitle('group_hippocampus')


plot_df <- grn_net_annot %>% 
    select(tf, n_common:frac_group_lc) %>% distinct() %>% 
    arrange(desc(n_group_lc)) %>% 
    filter(row_number()<=25) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf))) %>% 
    filter(n_group_lc>5)

p7 <- ggplot(plot_df, aes(tf, n_group_lc, fill=frac_group_lc)) +
    geom_bar(stat='identity') +
    rotate_x_text(40) +
    scale_fill_gradientn(colors=pals::brewer.rdpu(100)) +
    labs(fill='frac_specific') + ggtitle('group_lc')

p1 / (p2 + p3 + p4 + p6 + p7)
ggsave('plots/glut_trajectory_tfs_bar.pdf', width=25, height=8, bg='white')
ggsave('plots/glut_trajectory_tfs_bar.png', width=25, height=8, bg='white')



(p2 + p3 + p4 + p6 + p7)
ggsave('plots/glut_trajectory_top_tfs_bar.pdf', width=25, height=4, bg='white')
ggsave('plots/glut_trajectory_top_tfs_bar.png', width=25, height=4, bg='white')
