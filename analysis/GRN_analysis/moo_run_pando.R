source('~/scripts/single_cell/atac.R')
source('~/scripts/single_cell/wrapper.R')
source('~/scripts/single_cell/markers.R')
source('~/scripts/single_cell/celltype.R')
source('~/scripts/single_cell/de.R')
source('~/scripts/single_cell/graphs.R')
source('~/scripts/grn/models.R')

library(Signac)
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M, 
        lib.loc = '/local1/USERS/tomasgomes/multiome_analysis/')
library(doParallel)
library(Pando)

doParallel::registerDoParallel(40)

setwd('~/projects/axolotl/')

mome_atac <- read_rds('data/muo_glut_v1.2.2pando_ATAC_links10m_srt.rds')

mome_atac = infer_grn(
    mome_atac, 
    upstream = 10e+6, 
    downstream = 10e+6,
    only_tss = T,
    parallel = T,
    peak_to_gene_method = 'Signac',
    aggregate_peaks_col = 'cellclusters',
    verbose = 2,
    method = 'glm',
    alpha = 0.5
)

mome_atac %>% write_rds('data/muo_glut_v1.3.2pando_ATAC_links10m_srt.rds')

mome_atac = infer_grn(
    mome_atac, 
    upstream = 10e+6, 
    downstream = 10e+6,
    only_tss = T,
    parallel = T,
    peak_to_gene_method = 'Signac',
    aggregate_peaks_col = 'cellclusters',
    verbose = 2,
    method = 'bagging_ridge',
    alpha = 1,
    p_method = 'wilcox'
)

mome_atac %>% write_rds('data/muo_glut_v1.3.2pando_ATAC_links10m_srt.rds')

mome_atac = infer_grn(
    mome_atac, 
    upstream = 10e+6, 
    downstream = 10e+6,
    only_tss = T,
    parallel = T,
    peak_to_gene_method = 'Signac',
    aggregate_peaks_col = 'cellclusters',
    verbose = 2,
    method = 'glmnet',
    alpha = 0.5
)

mome_atac %>% write_rds('data/muo_glut_v1.3.2pando_ATAC_links10m_srt.rds')

