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
dir = '/links/groups/treutlein/USERS/tomasgomes/projects/pallium_evo/'


# read data (CellRanger peaks and MACS3 peaks, here only using the former)
mome_atac_SCT = read_rds('/links/groups/treutlein/USERS/tomasgomes/projects/pallium_evo/data/processed/multiome/multiome_integATAC_SCT.RDS')
mome_macs_SCT = readRDS('/links/groups/treutlein/USERS/tomasgomes/projects/pallium_evo/data/processed/multiome/multiome_integMACS_SCT.RDS')

# annotations (not directly used here)
## cell type annotations
meta = read.csv(paste0(dir,'data/annotations/axolotl_all_umeta.csv'), 
                header = T, row.names = 1)
### these cell names need reformatting
newcellnames = rownames(meta)
newcellnames = gsub('-1_1', '-a1-1', newcellnames)
newcellnames = gsub('-1_2', '-a1-2', newcellnames)
newcellnames = gsub('-1_3', '-a3-1', newcellnames)
newcellnames = gsub('-1_4', '-a3-2', newcellnames)
rownames(meta) = newcellnames
mome_atac_SCT = AddMetaData(mome_atac_SCT, metadata = meta)


## predicted regions (top is used)
meta_reg = read.csv(paste0(dir,'data/processed/multiome/WP_region_predictions.csv'), header = T, row.names = 1)
mome_atac_SCT = AddMetaData(mome_atac_SCT, metadata = meta_reg)


# plotting most important labels
dim_plot(mome_atac_SCT, reduction = 'umap', group.by = 'subclasses', label = T)
dim_plot(mome_atac_SCT, reduction = 'umap', group.by = 'cellclusters', label = T) +
    no_legend()
dim_plot(mome_atac_SCT, reduction = 'umap', group.by = c('pred_regions_top', 'subclasses'), label = T)


#### Combine objects ####
common_cells <- intersect(colnames(mome_atac_SCT), colnames(mome_macs_SCT))
mome_atac_SCT <- subset(mome_atac_SCT, cells=common_cells)
mome_macs_SCT <- subset(mome_macs_SCT, cells=common_cells)

mome_atac_SCT[['MACS']] <- mome_macs_SCT[['ATAC']]


#### Compare MACS and cellranger by UMAP ####
DefaultAssay(mome_atac_SCT) <- 'ATAC'
mome_atac_SCT <- mome_atac_SCT %>% 
    RunTFIDF(assay='ATAC') %>% 
    FindTopFeatures(min.cutoff='q50') %>% 
    RunSVD(reduction.name='alsi')

mome_atac_SCT <- mome_atac_SCT %>% 
    RunUMAP(dims=2:20, reduction='alsi', reduction.name='atacumap')


DefaultAssay(mome_atac_SCT) <- 'MACS'
mome_atac_SCT <- mome_atac_SCT %>% 
    RunTFIDF(assay='MACS') %>% 
    FindTopFeatures(min.cutoff='q50') %>% 
    RunSVD(reduction.name='mlsi')

mome_atac_SCT <- mome_atac_SCT %>% 
    RunUMAP(dims=2:20, reduction='mlsi', reduction.name='macsumap')

mome_atac_SCT <- mome_atac_SCT %>% 
    RunUMAP(dims=1:20, reduction='pca', reduction.name='rnaumap')


dim_plot(mome_atac_SCT, reduction = 'rnaumap', group.by = 'subclasses', label = T)
dim_plot(mome_atac_SCT, reduction = 'macsumap', group.by = 'subclasses', label = T)
dim_plot(mome_atac_SCT, reduction = 'atacumap', group.by = 'subclasses', label = T)


mome_atac_SCT %>% write_rds('data/muo_combined_annot_v1_srt.rds')
mome_atac_SCT <- read_rds('data/muo_combined_annot_v1_srt.rds')



#### Summarize highres clusters ####
mome_atac_SCT[['ATAC_bin']] <- CreateAssayObject(mome_atac_SCT[['ATAC']]@counts>0)
mome_atac_SCT[['MACS_bin']] <- CreateAssayObject(mome_atac_SCT[['MACS']]@counts>0)

mome_atac_SCT <- Pando::aggregate_assay(mome_atac_SCT, assay = 'ATAC_bin', group_name = 'cellclusters', slot = 'counts')
mome_atac_SCT <- Pando::aggregate_assay(mome_atac_SCT, assay = 'MACS_bin', group_name = 'cellclusters', slot = 'counts')



#### Get new annotations ####
new_annot <- read_csv('data/ref_glut_dat.csv') %>% 
    mutate(cell=str_replace_all(newcellnames, '_', '-')) %>% 
    filter(cell%in%colnames(mome_atac_SCT))

region_annot <- new_annot %>% 
    distinct(cell, region) %>% column_to_rownames('cell')

mome_atac_SCT <- subset(mome_atac_SCT, cells=unique(new_annot$cell))
mome_atac_SCT <- AddMetaData(mome_atac_SCT, region_annot)


DefaultAssay(mome_atac_SCT) <- 'ATAC'
mome_atac_SCT <- mome_atac_SCT %>% 
    RunTFIDF(assay='ATAC') %>% 
    FindTopFeatures(min.cutoff='q50') %>% 
    RunSVD(reduction.name='alsi')

mome_atac_SCT <- mome_atac_SCT %>% 
    RunUMAP(dims=2:20, reduction='alsi', reduction.name='atacumap')

DefaultAssay(mome_atac_SCT) <- 'MACS'
mome_atac_SCT <- mome_atac_SCT %>% 
    RunTFIDF(assay='MACS') %>% 
    FindTopFeatures(min.cutoff='q50') %>% 
    RunSVD(reduction.name='mlsi')

mome_atac_SCT <- mome_atac_SCT %>% 
    RunUMAP(dims=2:20, reduction='mlsi', reduction.name='macsumap')

mome_atac_SCT <- mome_atac_SCT %>% 
    RunUMAP(dims=1:20, reduction='pca', reduction.name='rnaumap')


dim_plot(mome_atac_SCT, reduction = 'rnaumap', group.by = c('region', 'subclasses'), label = T)
    
mome_atac_SCT %>% write_rds('data/muo_glut_v1_srt.rds')




#### Link peaks to genes ####
mome_atac_SCT <- FindVariableFeatures(mome_atac_SCT, assay='RNA', nfeatures=10000)

mome_atac_SCT <- RegionStats(
    mome_atac_SCT,
    genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M,
    assay = 'ATAC'
)

mome_atac_SCT <- RegionStats(
    mome_atac_SCT,
    genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M,
    assay = 'MACS'
)

mome_atac_SCT <- LinkPeaks(
    mome_atac_SCT,
    expression.assay = 'RNA',
    distance = 1e+6,
    peak.assay = 'ATAC',
    genes.use = VariableFeatures(mome_atac_SCT, assay='RNA')
)

mome_atac_SCT <- LinkPeaks(
    mome_atac_SCT,
    expression.assay = 'RNA',
    distance = 1e+6,
    peak.assay = 'MACS',
    genes.use = VariableFeatures(mome_atac_SCT, assay='RNA')
)

mome_atac_SCT %>% write_rds('data/muo_glut_v1.1links_srt.rds')
mome_atac_SCT <- read_rds('data/muo_glut_v1.1links_srt.rds')




#### Get motifs ####
library(JASPAR2020)

pfm = getMatrixSet(JASPAR2020, opts = list(tax_group = 'vertebrates', all_versions = F))
pfm %>% write_rds('~/resources/DB/JASPAR/JASPAR2020_vertebrates_motifs.rds')

pfm_list <- as.list(pfm)
motif_df <- map_dfr(pfm_list, function(x){
    tibble(
        motif = x@ID,
        tf_name = x@tags$remap_tf_name,
        tf_split = unlist(str_split(str_replace(x@name, '(.+)\\(.+\\)', '\\1'), '::')),
        name = x@name,
        class = x@matrixClass,
        collection = x@tags$collection,
        type = x@tags$type,
        symbol = x@tags$symbol
    )
})


motif_df <- motif_df %>% 
    mutate(
        tf = case_when(
            !is.na(tf_name) ~ tf_name,
            str_to_upper(tf_split) == symbol ~ symbol,
            T ~ str_to_upper(tf_split)
        )
    ) %>% select(motif, tf, everything())
    
motif_df %>% write_tsv('data/JASPAR2020_vertebrates_motif2tf.tsv')


#### Select peaks ####
pfm <- read_rds('~/resources/DB/JASPAR/JASPAR2020_vertebrates_motifs.rds')
motif_df <- read_tsv('data/JASPAR2020_vertebrates_motif2tf.tsv')

ATAC_links <- StringToGRanges(mome_atac_SCT@assays$ATAC@links$peak)
MACS_links <- StringToGRanges(mome_atac_SCT@assays$MACS@links$peak)

# Initiate GRN object and select candidate regions
mome_atac_SCT_plus = initiate_grn(
    mome_atac_SCT, 
    peak_assay = 'ATAC', 
    rna_assay = 'RNA', 
    regions = ATAC_links,
    exclude_exons = F
)

# Scan candidate regions for TF binding motifs
mome_atac_SCT_plus = find_motifs(
    mome_atac_SCT_plus, 
    pfm = pfm, 
    motif_tfs = motif_df,
    genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M
)

mome_atac_SCT_plus %>% write_rds('data/muo_glut_v1.2.1pando_ATAC_links_srt.rds')
# mome_atac_SCT_plus <- read_rds('data/muo_glut_v1.2.1pando_ATAC_links_srt.rds')


























