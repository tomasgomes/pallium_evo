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

setwd('~/projects/axolotl/')


#### Read stuff ####
mome_atac <- read_rds('data/muo_glut_v1.4annot_links10m_srt.rds')

grn_net <- read_tsv('data/grn/glm_modules.tsv')

fate_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_fate_da.tsv')
region_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_region_da.tsv')
group_da <- read_tsv('data/diff_expression/ATAC_linked_peaks_group_da.tsv')


reg_peaks_ranges <- grn_net$regions %>% StringToGRanges()
glm_coefs <- coef(mome_atac, network='glm_network') %>% 
    filter(padj<0.05) %>% 
    filter(!str_detect(target, '^(AMEX|LOC)'))


#### Plot regions ####
DefaultAssay(mome_atac) <- 'ATAC'
mome_atac_cov <- DietSeurat(mome_atac, assays=c('RNA', 'ATAC'))


#### EOMES ####
region <- FindRegion(mome_atac, 'EOMES') %>% 
    resize(width = 2, fix = 'end') %>%
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region, 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/EOMES_TSS.pdf', width=6, height=4)
ggsave('plots/tracks/EOMES_TSS.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr13q-334716878-334718671') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/EOMES_chr13q-334716878-334718671.pdf', width=6, height=4)
ggsave('plots/tracks/EOMES_chr13q-334716878-334718671.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr13q-343909677-343911142') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/EOMES_chr13q-343909677-343911142.pdf', width=6, height=4)
ggsave('plots/tracks/EOMES_chr13q-343909677-343911142.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr13q-335438834-335439733') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/EOMES_chr13q-335438834-335439733.pdf', width=6, height=4)
ggsave('plots/tracks/EOMES_chr13q-335438834-335439733.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr13q-335170478-335170721') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/EOMES_chr13q-335170478-335170721.pdf', width=6, height=4)
ggsave('plots/tracks/EOMES_chr13q-335170478-335170721.png', width=6, height=4)



#### NFIA ####
region <- FindRegion(mome_atac, 'NFIA') %>% 
    resize(width = 2, fix = 'end') %>%
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region, 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIA_TSS.pdf', width=6, height=4)
ggsave('plots/tracks/NFIA_TSS.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr1qs3-190962084-190963103') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIA_chr1qs3-190962084-190963103.pdf', width=6, height=4)
ggsave('plots/tracks/NFIA_chr1qs3-190962084-190963103.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr1qs3-198492918-198500888') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIA_chr1qs3-198492918-198500888.pdf', width=6, height=4)
ggsave('plots/tracks/NFIA_chr1qs3-198492918-198500888.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr1qs3-201931676-201934403') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIA_chr1qs3-201931676-201934403.pdf', width=6, height=4)
ggsave('plots/tracks/NFIA_chr1qs3-201931676-201934403.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr1qs3-198083625-198089074') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIA_chr1qs3-198083625-198089074.pdf', width=6, height=4)
ggsave('plots/tracks/NFIA_chr1qs3-198083625-198089074.png', width=6, height=4)



#### NFIC ####
region <- FindRegion(mome_atac, 'NFIC') %>% 
    resize(width = 2, fix = 'start') %>%
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region, 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIC_TSS.pdf', width=6, height=4)
ggsave('plots/tracks/NFIC_TSS.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr1ps1-494147302-494157648') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIC_chr1ps1-494147302-494157648.pdf', width=6, height=4)
ggsave('plots/tracks/NFIC_chr1ps1-494147302-494157648.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr1ps1-497477772-497478826') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIC_chr1ps1-497477772-497478826.pdf', width=6, height=4)
ggsave('plots/tracks/NFIC_chr1ps1-497477772-497478826.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr1ps1-493846572-493850251') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIC_chr1ps1-493846572-493850251.pdf', width=6, height=4)
ggsave('plots/tracks/NFIC_chr1ps1-493846572-493850251.png', width=6, height=4)




region <- FindRegion(mome_atac, 'chr1qs3-198083625-198089074') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NFIC_chr1qs3-198083625-198089074.pdf', width=6, height=4)
ggsave('plots/tracks/NFIC_chr1qs3-198083625-198089074.png', width=6, height=4)

#### BCL11B ####
region <- FindRegion(mome_atac, 'BCL11B') %>% 
    resize(width = 2, fix = 'start') %>%
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region, 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/BCL11B_TSS.pdf', width=6, height=4)
ggsave('plots/tracks/BCL11B_TSS.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr14p-147300753-147316598') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/BCL11B_chr14p-147300753-147316598.pdf', width=6, height=4)
ggsave('plots/tracks/BCL11B_chr14p-147300753-147316598.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr14p-147321619-147325851') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/BCL11B_chr14p-147321619-147325851.pdf', width=6, height=4)
ggsave('plots/tracks/BCL11B_chr14p-147321619-147325851.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr14p-147388548-147391157') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/BCL11B_chr14p-147388548-147391157.pdf', width=6, height=4)
ggsave('plots/tracks/BCL11B_chr14p-147388548-147391157.png', width=6, height=4)


#### MEIS1 ####
region <- FindRegion(mome_atac, 'MEIS1') %>% 
    resize(width = 2, fix = 'end') %>%
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region, 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/MEIS1_TSS.pdf', width=6, height=4)
ggsave('plots/tracks/MEIS1_TSS.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr2qs2-191117466-191117950') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/MEIS1_chr2qs2-191117466-191117950.pdf', width=6, height=4)
ggsave('plots/tracks/MEIS1_chr2qs2-191117466-191117950.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr2qs2-191554729-191555073') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/MEIS1_chr2qs2-191554729-191555073.pdf', width=6, height=4)
ggsave('plots/tracks/MEIS1_chr2qs2-191554729-191555073.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr2qs2-201615884-201617168') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/MEIS1_chr2qs2-201615884-201617168.pdf', width=6, height=4)
ggsave('plots/tracks/MEIS1_chr2qs2-201615884-201617168.png', width=6, height=4)


#### NEUROD4 ####
region <- FindRegion(mome_atac, 'NEUROD4') %>% 
    resize(width = 2, fix = 'end') %>%
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region, 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NEUROD4_TSS.pdf', width=6, height=4)
ggsave('plots/tracks/NEUROD4_TSS.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr3qs1-388732470-388733502') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NEUROD4_chr3qs1-388732470-388733502.pdf', width=6, height=4)
ggsave('plots/tracks/NEUROD4_chr3qs1-388732470-388733502.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr3qs1-390085496-390086748') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/NEUROD4_chr3qs1-390085496-390086748.pdf', width=6, height=4)
ggsave('plots/tracks/NEUROD4_chr3qs1-390085496-390086748.png', width=6, height=4)



#### GLI2 ####
region <- FindRegion(mome_atac, 'GLI2') %>% 
    resize(width = 2, fix = 'end') %>%
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region, 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI2_TSS.pdf', width=6, height=4)
ggsave('plots/tracks/GLI2_TSS.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr9qs2-302041910-302042510') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI2_chr9qs2-302041910-302042510.pdf', width=6, height=4)
ggsave('plots/tracks/GLI2_chr9qs2-302041910-302042510.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr9qs2-302129609-302130336') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI2_chr9qs2-302129609-302130336.pdf', width=6, height=4)
ggsave('plots/tracks/GLI2_chr9qs2-302129609-302130336.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr9qs2-302140878-302145286') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI2_chr9qs2-302140878-302145286.pdf', width=6, height=4)
ggsave('plots/tracks/GLI2_chr9qs2-302140878-302145286.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr9qs2-294764275-294765790') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI2_chr9qs2-294764275-294765790.pdf', width=6, height=4)
ggsave('plots/tracks/GLI2_chr9qs2-294764275-294765790.png', width=6, height=4)


#### GLI3 ####
region <- FindRegion(mome_atac, 'GLI3') %>% 
    resize(width = 2, fix = 'end') %>%
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region, 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI3_TSS.pdf', width=6, height=4)
ggsave('plots/tracks/GLI3_TSS.png', width=6, height=4)


region <- FindRegion(mome_atac, 'chr5ps2-303730325-303731442') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI3_chr5ps2-303730325-303731442.pdf', width=6, height=4)
ggsave('plots/tracks/GLI3_chr5ps2-303730325-303731442.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr5ps2-306248977-306249868') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI3_chr5ps2-306248977-306249868.pdf', width=6, height=4)
ggsave('plots/tracks/GLI3_chr5ps2-306248977-306249868.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr5ps2-306323709-306324349') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI3_chr5ps2-306323709-306324349.pdf', width=6, height=4)
ggsave('plots/tracks/GLI3_chr5ps2-306323709-306324349.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr5ps2-306236045-306237200') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI3_chr5ps2-306236045-306237200.pdf', width=6, height=4)
ggsave('plots/tracks/GLI3_chr5ps2-306236045-306237200.png', width=6, height=4)



region <- FindRegion(mome_atac, 'chr5ps2-306423351-306425922') %>% 
    Extend(upstream=10000, downstream=10000)

CoveragePlot(
    mome_atac_cov, 
    region = region[1,], 
    group.by = 'region',
    ranges = reg_peaks_ranges
) 

ggsave('plots/tracks/GLI3_chr5ps2-306423351-306425922.pdf', width=6, height=4)
ggsave('plots/tracks/GLI3_chr5ps2-306423351-306425922.png', width=6, height=4)










