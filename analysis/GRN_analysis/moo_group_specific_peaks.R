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

new_trajectories_glut <- read_csv('ref_glut_dat.csv') %>% 
    mutate(cell=str_replace_all(newcellnames, '_', '-')) %>% 
    filter(cell%in%colnames(mome_atac)) %>% 
    filter(str_detect(cellclusters, 'glut_')) %>% 
    distinct(newcellnames, group, cell) %>% 
    column_to_rownames('cell')



#### Get group-specific peaks ####

linked_peaks <- Links(mome_atac) %>% 
    as_tibble() %>% mutate(dist=width(Links(mome_atac)))

group_da_sig <- group_da %>% 
    filter(padj<0.05, coef>0, detect_ratio>1.5) 

group_spec_peaks <- group_da_sig %>% 
    group_by(feature) %>% 
    filter(detect_ratio==max(detect_ratio)) %>% 
    distinct(group, feature, detect_self, detect_ratio, padj) %>% 
    inner_join(linked_peaks, by=c('feature'='peak'))


#### Get peak ranges for all ####

mome_groups <- list(
    group_111 = subset(mome_atac, group_111=='group_111'),
    group_8620_ep = subset(mome_atac, group_8620_ep=='group_8620_ep'),
    group_eomes = subset(mome_atac, group_eomes=='group_eomes'),
    group_lc = subset(mome_atac, group_lc=='group_lc'),
    group_hippocampus = subset(mome_atac, group_hippocampus=='group_hippocampus')
)

mome_glut_cells <- colnames(mome_atac)[str_detect(mome_atac$cellclusters, 'glut_') & colnames(mome_atac)%in%rownames(new_trajectories_glut)]
mome_glut <- subset(mome_atac, cells=mome_glut_cells)
mome_glut <- AddMetaData(mome_glut, new_trajectories_glut)

peaks_use <- group_spec_peaks %>% 
    filter(!str_detect(gene, '^(AMEX|LOC)')) %>% 
    filter(dist<10000, detect_ratio>2)

da_peak_regions <- peaks_use$feature %>% 
    StringToGRanges() 
    # resize(width=2, fix='center') %>% 
    # Extend(upstream = 1000, downstream = 1000)

plots <- map(1:length(da_peak_regions), function(n){
    p <- CoveragePlot(
            mome_glut, 
            region=da_peak_regions[n], 
            annotation=F, 
            peaks=F, 
            links=F, 
            group.by='group',
            window=500
    ) 
    return(p + ggtitle(da_peak_regions[n]$gene))
})

wrap_plots(plots, ncol=1) & 
    theme_void() & no_legend() & theme(strip.text = element_blank()) &
    scale_fill_manual(values=rep('black', 5))


# Good looking peaks
# c(4, 5, 6, 8, 12, 17, 20, 26, 29, 31, 49, 53)


da_peak_regions <- peaks_use$feature %>% 
    StringToGRanges() %>% 
    resize(width=2, fix='center') %>%
    Extend(upstream = 6000, downstream = 6000)

plots <- map(1:length(da_peak_regions[c(4, 5, 6, 8, 12, 17, 20, 26, 29, 31, 49, 53)]), function(n){
    p <- CoveragePlot(
        mome_glut, 
        region=da_peak_regions[n], 
        annotation=F, 
        peaks=F, 
        links=F, 
        group.by='group',
        window=500
    ) 
    return(p + ggtitle(da_peak_regions[n]$gene))
})
    
    
wrap_plots(plots) & 
    theme_void() & no_legend() & theme(strip.text = element_blank()) &
    scale_fill_manual(values=rep('black', 5))

peaks_use_here <- peaks_use[c(4, 5, 6, 8, 12, 17, 20, 26, 29, 31, 49, 53), ]
da_regs_use <- da_peak_regions[c(4, 5, 6, 8, 12, 17, 20, 26, 29, 31, 49, 53)]


#### ST18 ####
region_use <- da_regs_use[1] %>% 
    resize(width=2, fix='center') %>%
    Extend(upstream = 3000, downstream = 6000)

p1 <- CoveragePlot(
    mome_glut, 
    region=region_use, 
    annotation=F, 
    peaks=F, 
    links=F, 
    group.by='group',
    window=500
) 

p2 <- AnnotationPlot(mome_glut, region_use) + ggtitle(peaks_use_here[1,]$gene)

p1 / p2 + plot_layout(heights=c(5,1))
ggsave('plots/tracks/ST18_TSS.png', width=8, height=8)
ggsave('plots/tracks/ST18_TSS.pdf', width=8, height=8)


#### ST18 ####
region_use <- da_regs_use[2] %>% 
    resize(width=2, fix='center') %>%
    Extend(upstream = 6000, downstream = 6000)

p1 <- CoveragePlot(
    mome_glut, 
    region=region_use, 
    annotation=F, 
    peaks=F, 
    links=F, 
    group.by='group',
    window=500
) 

p2 <- AnnotationPlot(mome_glut, region_use) + ggtitle(peaks_use_here[2,]$gene)

p1 / p2 + plot_layout(heights=c(5,1))
ggsave('plots/tracks/ZFPM2_TSS.png', width=8, height=8)
ggsave('plots/tracks/ZFPM2_TSS.pdf', width=8, height=8)


#### ST18 ####
region_use <- da_regs_use[3] %>% 
    resize(width=2, fix='center') %>%
    Extend(upstream = 6000, downstream = 6000)

p1 <- CoveragePlot(
    mome_glut, 
    region=region_use, 
    annotation=F, 
    peaks=F, 
    links=F, 
    group.by='group',
    window=500
) 

p2 <- AnnotationPlot(mome_glut, region_use) + ggtitle(peaks_use_here[3,]$gene)

p1 / p2 + plot_layout(heights=c(5,1))
ggsave('plots/tracks/DRD3_TSS.png', width=8, height=8)
ggsave('plots/tracks/DRD3_TSS.pdf', width=8, height=8)


#### ST18 ####
region_use <- da_regs_use[4] %>% 
    resize(width=2, fix='center') %>%
    Extend(upstream = 6000, downstream = 3000)

p1 <- CoveragePlot(
    mome_glut, 
    region=region_use, 
    annotation=F, 
    peaks=F, 
    links=F, 
    group.by='group',
    window=500
) 

p2 <- AnnotationPlot(mome_glut, region_use) + ggtitle(peaks_use_here[4,]$gene)

p1 / p2 + plot_layout(heights=c(5,1))
ggsave('plots/tracks/PDYN_TSS.png', width=8, height=8)
ggsave('plots/tracks/PDYN_TSS.pdf', width=8, height=8)


#### ST18 ####
region_use <- da_regs_use[12] %>% 
    resize(width=2, fix='center') %>%
    Extend(upstream = 2000, downstream = 4000)

p1 <- CoveragePlot(
    mome_glut, 
    region=region_use, 
    annotation=F, 
    peaks=F, 
    links=F, 
    group.by='group',
    window=500
) 

p2 <- AnnotationPlot(mome_glut, region_use) + ggtitle(peaks_use_here[12,]$gene)

p1 / p2 + plot_layout(heights=c(5,1))
ggsave('plots/tracks/NMB_TSS.png', width=8, height=8)
ggsave('plots/tracks/NMB_TSS.pdf', width=8, height=8)







