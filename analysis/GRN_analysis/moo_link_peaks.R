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
library(Pando)
library(ggplot2)
library(igraph)

setwd('~/projects/axolotl/')

mome_atac <- read_rds('data/muo_glut_v1_srt.rds') 
pfm <- read_rds('~/resources/DB/JASPAR/JASPAR2020_vertebrates_motifs.rds')
motif_df <- read_tsv('data/JASPAR2020_vertebrates_motif2tf.tsv')


#### Link peaks to genes ####
mome_atac <- FindVariableFeatures(mome_atac, assay='RNA', nfeatures=10000)

mome_atac <- RegionStats(
    mome_atac,
    genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M,
    assay = 'ATAC'
)

mome_atac <- RegionStats(
    mome_atac,
    genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M,
    assay = 'MACS'
)

mome_atac <- LinkPeaks(
    mome_atac,
    expression.assay = 'RNA',
    distance = 10e+6,
    peak.assay = 'ATAC',
    genes.use = VariableFeatures(mome_atac, assay='RNA')
)

mome_atac <- LinkPeaks(
    mome_atac,
    expression.assay = 'RNA',
    distance = 10e+6,
    peak.assay = 'MACS',
    genes.use = VariableFeatures(mome_atac, assay='RNA')
)

mome_atac %>% write_rds('data/muo_glut_v1.2links10m_srt.rds')


#### Select peaks ####
ATAC_links <- StringToGRanges(mome_atac@assays$ATAC@links$peak)

# Initiate GRN object and select candidate regions
mome_atac_plus = initiate_grn(
    mome_atac, 
    peak_assay = 'ATAC', 
    rna_assay = 'RNA', 
    regions = ATAC_links,
    exclude_exons = F
)

# Scan candidate regions for TF binding motifs
mome_atac_plus = find_motifs(
    mome_atac_plus, 
    pfm = pfm, 
    motif_tfs = motif_df,
    genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M
)

mome_atac_plus %>% write_rds('data/muo_glut_v1.2.2pando_ATAC_links10m_srt.rds')





