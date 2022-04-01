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

grn_graph <- as_tbl_graph(grn_net)

fate_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_fate_da.tsv')
region_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_region_da.tsv')
group_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_group_da.tsv')

group_genes <- read_rds('heatmap_gen/genes_fates.rds')

diff_genes <- read_rds('diff_traj_genes.rds')
names(diff_genes)[2] <- '8620_ep'

fate_genes <- list(
    glut_SUBSET_1 = fate_genes$`111`$glut_SUBSET_1,
    glut_SUBSET_11 = fate_genes$`111`$glut_SUBSET_11,
    glut_SUBSET_6 = fate_genes$`8620_ep`$glut_SUBSET_6,
    glut_SUBSET_8 = fate_genes$`8620_ep`$glut_SUBSET_8,
    glut_SUBSET_10 = fate_genes$eomes$glut_SUBSET_10,
    glut_SUBSET_22 = fate_genes$eomes$glut_SUBSET_22,
    glut_SUBSET_0 = fate_genes$hippocampus$glut_SUBSET_0,
    glut_SUBSET_7 = fate_genes$hippocampus$glut_SUBSET_7,
    glut_SUBSET_2 = fate_genes$lc$glut_SUBSET_2,
    glut_SUBSET_9 = fate_genes$lc$glut_SUBSET_9
)

#### Stats on the full network ####
net_stats <- grn_net %>% 
    group_by(tf) %>% 
    summarize(ngenes=n()) %>% 
    arrange(desc(ngenes)) %>% 
    mutate(tf=str_to_title(tf)) %>% 
    mutate(tf=factor(tf, levels=unique(.$tf))) %>% 
    top_n(20, ngenes)

ggplot(net_stats, aes(tf, ngenes)) +
    geom_bar(stat='identity', color='black', fill='grey', size=0.2) +
    rotate_x_text(40) +
    article_text() + theme_rangeframe() + scale_axis_rangeframe() +
    scale_y_continuous(expand=c(0,0))
ggsave('plots/glut_grn_tf_top20_bar.pdf', width=56, height=29, unit='mm')


#### Get fate-specific networks ####


fate_da_peaks <- fate_da %>% 
    filter(padj<0.05, coef>0, detect_ratio>1.5) %>% 
    group_by(group) %>% group_split()

group_da_peaks <- group_da %>% 
    filter(padj<0.05, coef>0, detect_ratio>1.5) %>% 
    group_by(group) %>% group_split()

group_nets <- map(group_da_peaks, function(peaks){
    grn_net %>% 
        filter(regions%in%peaks$feature) %>% 
        group_by(tf) %>% 
        mutate(ngenes=n(), group=peaks$group[1]) %>% 
        group_by(target) %>% 
        mutate(ntfs=n()) %>% 
        return()
})

plots <- group_nets %>% map(function(x){
    group_genes <- diff_genes[[str_remove(x$group[1], 'group_')]]
    plot_df <- x %>% 
        ungroup() %>% 
        distinct(tf, ngenes) %>% 
        arrange(desc(ngenes)) %>% 
        filter(row_number()<=20) %>% 
        mutate(tf=factor(tf, levels=unique(.$tf)))
    p <- ggplot(plot_df, aes(tf, ngenes, fill=tf%in%group_genes)) +
        geom_bar(stat='identity', fill='darkgrey', color='black', size=0.3) +
        scale_fill_manual(values=c('darkgrey', 'black')) +
        rotate_x_text(40) + no_legend() +
        scale_y_continuous(expand=c(0,0)) +
        theme(text = element_text(family='Helvetica')) +
        article_text() +
        theme_rangeframe() + scale_axis_rangeframe() +
        ggtitle(x$group[1])
    return(p)
})

wrap_plots(plots)

ggsave('plots/glut_trajectory_grn_tf_bar.pdf', width=10, height=4)
ggsave('plots/glut_trajectory_grn_tf_bar.png', width=35, height=12, bg = 'white')



fate_nets <- map(fate_da_peaks, function(peaks){
grn_net %>% 
        filter(regions%in%peaks$feature) %>% 
        group_by(tf) %>% 
        mutate(ngenes=n(), group=peaks$group[1]) %>% 
        group_by(target) %>% 
        mutate(ntfs=n()) %>% 
        return()
})

plots <- fate_nets %>% map(function(x){
    group_genes <- fate_genes[[x$group[1]]]
    plot_df <- x %>% 
        distinct(tf, ngenes) %>% 
        arrange(desc(ngenes)) %>% mutate(tf=factor(tf, levels=unique(.$tf)))
    p <- ggplot(plot_df, aes(tf, ngenes, fill=tf%in%group_genes)) +
        geom_bar(stat='identity') +
        scale_fill_manual(values=c('darkgrey', 'black')) +
        rotate_x_text(40) + no_legend() +
        ggtitle(x$group[1])
    return(p)
})

wrap_plots(plots)
ggsave('plots/glut_fate_grn_tf_bar.png', width=48, height=20, bg = 'white')







#### Center network around specific genes ####
all_da_peaks <- map(group_da_peaks, function(x) x$feature) %>% purrr::reduce(union)
all_genes <- grn_graph %N>% as_tibble() %>% pull(name)

genes_plot <- c('NFIA', 'NFIX', 'NFIC', 'NFIB', 'GLI2', 'GLI3', 'EOMES', 'MEIS2', 'STAT5B')

gene_graphs <- map(genes_plot, function(gene){
    
    spaths <- all_shortest_paths(grn_graph, gene, all_genes, mode='out')$res
    spath_list <- Pando::map_par(spaths, function(p){
        edg <- names(p)
        edg_graph <- grn_graph %N>%
            filter(name%in%edg) %>%
            convert(to_shortest_path, from=which(.N()$name==edg[1]), to=which(.N()$name==edg[length(edg)])) %E>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            as_tibble()
        edg_dir <- edg_graph %>% pull(estimate) %>% sign() %>% prod()
        edg_p <- edg_graph %>% pull(padj) %>% {-log10(.)} %>% mean()
        return(
            list(
                path = tibble(
                    start_node = edg[1],
                    end_node = edg[length(edg)],
                    dir = edg_dir,
                    path = paste(edg, collapse=';'),
                    path_regions = paste(edg_graph$regions, collapse=';'),
                    order = length(edg)-1,
                    mean_padj = edg_p
                ),
                graph = mutate(edg_graph, path=paste(edg, collapse=';'), end_node=edg[length(edg)], comb_dir=edg_dir)
            )
        )
    }, parallel = T)
    
    spath_dir <- map_dfr(spath_list, function(x) x$path) %>%
        mutate(
            path_genes=str_split(path, ';'),
            path_regions=str_split(path_regions, ';')
        )
    spath_graph <- map_dfr(spath_list, function(x) x$graph)
    
    grn_pruned <- spath_dir %>%
        select(start_node, end_node, everything()) %>%
        group_by(end_node) %>% filter(order<=2) %>% 
        filter(order==1 | mean_padj==max(mean_padj) | any(path_regions%in%all_da_peaks))
    
    spath_graph_pruned <- spath_graph %>%
        filter(path%in%grn_pruned$path) %>%
        select(from_node, to_node, end_node, comb_dir) %>% distinct()
    
    grn_graph_pruned <- grn_graph %E>%
        mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
        as_tibble() %>% distinct() %>%
        inner_join(spath_graph_pruned) %>%
        select(from_node, to_node, everything(), -from, -to) %>% arrange(comb_dir, regions%in%all_da_peaks) %>% as_tbl_graph()
    
    return(grn_graph_pruned)
})

gene_graphs <- set_names(gene_graphs, genes_plot)
map(names(gene_graphs), function(n){
    g <- gene_graphs[[n]]
    p <- ggraph(g, layout='tree', circular=T) +
        geom_edge_diagonal(aes(color=sign(estimate)), width=0.3) +
        geom_node_label(aes(label=name, filter=name%in%tfs$symbol), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fill='black', color='white') +
        geom_node_label(aes(label=name, filter=!name%in%tfs$symbol), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05) +
        scale_edge_color_gradientn(colors=c('#f5b7b1', '#7dcea0')) +
        scale_x_continuous(expand=c(0.1, 0)) +
        scale_y_continuous(expand=c(0.1, 0)) +
        article_text() + theme_void() + no_legend() 
    
    print(p)
    ggsave(paste0('plots/trees/', n, '_grn_ctree.png'), width=8, height=8, bg='white')
    ggsave(paste0('plots/trees/', n, '_grn_ctree.pdf'), width=8, height=8, bg='white')
})


gene_graphs <- set_names(gene_graphs, genes_plot)
map(names(gene_graphs), function(n){
    g <- gene_graphs[[n]]
    p <- ggraph(g, layout='tree', circular=T) +
        geom_edge_diagonal(aes(color=sign(estimate)), width=0.3) +
        geom_node_point(aes(filter=!name%in%tfs$symbol), size=1, shape=21, fill='grey', stroke=0.2) +
        geom_node_label(aes(label=str_to_title(name), filter=name%in%tfs$symbol), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fontface='italic') +
        scale_edge_color_gradientn(colors=c('#f5b7b1', '#7dcea0')) +
        scale_x_continuous(expand=c(0.1, 0)) +
        scale_y_continuous(expand=c(0.1, 0)) +
        theme(
            text=element_text(face='italic', size=5)
        ) +
        article_text() + theme_void() + no_legend() 
    print(p)
    ggsave(paste0('plots/trees/', n, '_grn_2layer_ctree.png'), width=8, height=8, bg='white')
    ggsave(paste0('plots/trees/', n, '_grn_2layer_ctree.pdf'), width=8, height=8, bg='white')
})




##### Group-specific trees ####

get_gene_network <- function(
    graph, gene, max_order=2, unique_paths=T
){
    gnet <- as_tbl_graph(graph)
    all_genes <- gnet %N>% as_tibble() %>% pull(name)
    spaths <- all_shortest_paths(gnet, gene, all_genes, mode='out')$res
    spath_list <- Pando::map_par(spaths, function(p){
        edg <- names(p)
        edg_graph <- gnet %N>%
            filter(name%in%edg) %>%
            convert(to_shortest_path, from=which(.N()$name==edg[1]), to=which(.N()$name==edg[length(edg)])) %E>%
            mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
            as_tibble()
        edg_dir <- edg_graph %>% pull(estimate) %>% sign() %>% prod()
        edg_p <- edg_graph %>% pull(padj) %>% {-log10(.)} %>% mean()
        return(
            list(
                path = tibble(
                    start_node = edg[1],
                    end_node = edg[length(edg)],
                    dir = edg_dir,
                    path = paste(edg, collapse=';'),
                    order = length(edg)-1,
                    mean_padj = edg_p
                ),
                graph = mutate(edg_graph, path=paste(edg, collapse=';'), end_node=edg[length(edg)], comb_dir=edg_dir)
            )
        )
    }, parallel = T)
    
    spath_dir <- map_dfr(spath_list, function(x) x$path) %>%
        mutate(path_genes=str_split(path, ';'))
    spath_graph <- map_dfr(spath_list, function(x) x$graph)
    
    
    grn_pruned <- spath_dir %>%
        select(start_node, end_node, everything()) %>%
        group_by(end_node) %>% filter(order<=max_order) 
    
    if (unique_paths){
        grn_pruned <- filter(grn_pruned, order==1 | mean_padj==max(mean_padj))
    }
    
    spath_graph_pruned <- spath_graph %>%
        filter(path%in%grn_pruned$path) %>%
        select(from_node, to_node, end_node, comb_dir) %>% distinct()
    
    gnet <- gnet %E>%
        mutate(from_node=.N()$name[from], to_node=.N()$name[to]) %>%
        as_tibble() %>% distinct() %>%
        inner_join(spath_graph_pruned) %>%
        select(from_node, to_node, everything(), -from, -to) %>% arrange(comb_dir) %>% as_tbl_graph() 
    
    return(gnet)
}


graph_use <- get_gene_network(group_nets[[2]], 'FOXP1')

ggraph(graph_use, layout='tree', circular=T) +
    geom_edge_diagonal(aes(color=sign(estimate)), width=0.4) +
    # geom_node_point(size=1, shape=21, fill='grey', stroke=0.2) +
    # geom_node_label(aes(label=str_to_title(name), filter=name%in%tfs$symbol), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fill='black', color='white', fontface='italic') +
    # geom_node_label(aes(label=str_to_title(name), filter=name%in%tfs$symbol), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fontface='italic') +            # geom_node_label(aes(label=str_to_title(name), filter=name%in%tfs$symbol), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fontface='italic') +
    geom_node_label(aes(label=str_to_title(name)), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fontface='italic') +
    scale_edge_color_gradientn(colors=c('#f5b7b1', '#7dcea0')) +
    scale_x_continuous(expand=c(0.1, 0)) +
    scale_y_continuous(expand=c(0.1, 0)) +
    article_text() + theme_void() + no_legend() 
ggsave('plots/trees/RFX4_lc.png', width=5, height=5)
ggsave('plots/trees/RFX4_lc.pdf', width=5, height=5)


ggraph(graph_use, layout='tree', circular=T) +
    geom_edge_diagonal(aes(color=sign(estimate)), width=0.4) +
    geom_node_point(size=1, shape=21, fill='grey', stroke=0.2) +
    # geom_node_label(aes(label=str_to_title(name), filter=name%in%tfs$symbol), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fill='black', color='white', fontface='italic') +
    geom_node_label(aes(label=str_to_title(name), filter=name%in%tfs$symbol), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fontface='italic') +            # geom_node_label(aes(label=str_to_title(name), filter=name%in%tfs$symbol), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fontface='italic') +
    # geom_node_label(aes(label=str_to_title(name)), size=5/ggplot2::.pt, label.padding=unit(0.05, 'cm'), label.size=0.05, fontface='italic') +
    scale_edge_color_gradientn(colors=c('#f5b7b1', '#7dcea0'), values = c(-1,1)) +
    scale_x_continuous(expand=c(0.1, 0)) +
    scale_y_continuous(expand=c(0.1, 0)) +
    article_text() + theme_void() + no_legend() 
ggsave('plots/trees/RFX4_lc_nolab.png', width=5, height=5)
ggsave('plots/trees/RFX4_lc_nolab.pdf', width=5, height=5)




