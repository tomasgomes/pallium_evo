# libraries
library(Signac)
library(BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M, 
        lib.loc = "/local1/USERS/tomasgomes/multiome_analysis/")
library(Pando)
library(ggplot2)
library(igraph)


# root directory
dir = "/links/groups/treutlein/USERS/tomasgomes/projects/pallium_evo/"


# read data (CellRanger peaks and MACS3 peaks, here only using the former)
mome_atac_SCT = readRDS(file = paste0(dir,
                                      "data/processed/multiome/multiome_integATAC_SCT.RDS"))
mome_macs_SCT = readRDS(file = paste0(dir,
                                      "data/processed/multiome/multiome_integMACS_SCT.RDS"))


# annotations (not directly used here)
## cell type annotations
meta = read.csv(paste0(dir,"data/annotations/axolotl_all_umeta.csv"), 
                header = T, row.names = 1)
### these cell names need reformatting
newcellnames = rownames(meta)
newcellnames = gsub("-1_1", "-a1-1", newcellnames)
newcellnames = gsub("-1_2", "-a1-2", newcellnames)
newcellnames = gsub("-1_3", "-a3-1", newcellnames)
newcellnames = gsub("-1_4", "-a3-2", newcellnames)
rownames(meta) = newcellnames
mome_atac_SCT = AddMetaData(mome_atac_SCT, metadata = meta)

## predicted regions (top is used)
meta_reg = read.csv(paste0(dir,"data/processed/multiome/WP_region_predictions.csv"), 
                    header = T, row.names = 1)
mome_atac_SCT = AddMetaData(mome_atac_SCT, metadata = meta_reg)


# plotting most important labels
DimPlot(mome_atac_SCT, reduction = "umap", group.by = "subclasses", label = T)
DimPlot(mome_atac_SCT, reduction = "umap", group.by = "cellclusters", label = T)
DimPlot(mome_atac_SCT, reduction = "umap", group.by = "pred_regions_top", label = T)


# Load and prepare marker genes for all cell types
mk_ct_l = readRDS(paste0(dir,
                         "./results/RegionAnalysis/mk_ct_allData_subsets.RDS"))
mk_ct = mk_ct_l$all[mk_ct_l$all$padj<=0.05 & mk_ct_l$all$logFC>0.15,]

## compile markers from subsets
mk_ct_comp = rbind(mk_ct_l$GABA, mk_ct_l$glut, mk_ct_l$ependymal, mk_ct_l$npc,
                   mk_ct_l$all[mk_ct_l$all$group %in% c("microglia_8", "oligodendrocyte_10","oligodendrocyte_15", 
                                                        "endothelial_11", "endothelial_12", "endothelial_14"),])
mk_ct_comp = mk_ct_comp[mk_ct_comp$padj<=0.05 & mk_ct_comp$logFC>0.15,]


# Split ATAC into peaks from 3 types of regions CellRanger peaks
DefaultAssay(mome_atac_SCT) = "ATAC"
transc = Annotation(mome_atac_SCT)[Annotation(mome_atac_SCT)$type=="transcript",] # get transcripts
transc = granges(transc[width(transc)>600,])
proms = restrict(promoters(transc, upstream=5000, downstream=500), start = 1) # restrict fixes issues with negative genome coordinates
proms$type = "promoter"
genebodies = setdiff(transc, proms)
genebodies$type = "genebody"
feat = c(proms, genebodies)

closest_f = ClosestFeature(mome_atac_SCT, regions = rownames(mome_atac_SCT), annotation = feat)
closest_f$type = ifelse(closest_f$distance!=0, "distal", closest_f$type)

for(tt in unique(closest_f$type)){
  p = closest_f$query_region[closest_f$type==tt]
  i = paste0("ATAC_", tt)
  mome_atac_SCT[[i]] = CreateChromatinAssay(mome_atac_SCT@assays$ATAC@counts[p,])
  mome_atac_SCT = RunTFIDF(mome_atac_SCT, assay = i)
  mome_atac_SCT = FindTopFeatures(mome_atac_SCT, assay=i, min.cutoff = 300)
  print(length(mome_atac_SCT[[i]]@var.features))
}


# Basic Pando run
DefaultAssay(mome_atac_SCT) = "RNA"

# Get motif data
data(motifs)

# define genes
genes = VariableFeatures(mome_atac_SCT, assay = "RNA")
genes = unique(mk_ct_comp$feature)
length(genes)

# define peaks
mome_atac_SCT = FindTopFeatures(mome_atac_SCT, assay = "ATAC_promoter", min.cutoff = 450)
mome_atac_SCT = FindTopFeatures(mome_atac_SCT, assay = "ATAC_distal", min.cutoff = 450)
mome_atac_SCT = FindTopFeatures(mome_atac_SCT, assay = "ATAC_genebody", min.cutoff = 450)
regions = unique(c(mome_atac_SCT@assays$ATAC_promoter@var.features,
                   mome_atac_SCT@assays$ATAC_distal@var.features,
                   mome_atac_SCT@assays$ATAC_genebody@var.features))

isreg = rownames(mome_atac_SCT@assays$ATAC@meta.features) %in% regions
regions = mome_atac_SCT@assays$ATAC@ranges[isreg,]
length(regions)

# Initiate GRN object and select candidate regions
mome_atac_SCT_plus = initiate_grn(mome_atac_SCT, peak_assay = "ATAC", rna_assay = "RNA",
                                  regions = regions)

# Scan candidate regions for TF binding motifs
mome_atac_SCT_plus = find_motifs(mome_atac_SCT_plus, pfm = motifs, 
                                 genome = BSgenome.Amexicanum.axolotlomics.AmexGv6cut500M)

DefaultAssay(mome_atac_SCT_plus) = "ATAC"
# Infer gene regulatory network
cl = parallel::makeCluster(32)
doParallel::registerDoParallel(cl)
mome_atac_SCT_plus = infer_grn(mome_atac_SCT_plus, downstream = 500, parallel = T)

# Print inferred coefficients
coef(mome_atac_SCT_plus)

# save
saveRDS(mome_atac_SCT_plus, 
        file = paste0(dir, "data/processed/multiome/multiome_integATAC_SCT_plus.RDS"))
write.csv(coef(mome_atac_SCT_plus), 
          file = paste0(dir, "results/multiome/mome_atac_SCT_plus_coef.csv"), 
          col.names = T, row.names = F, quote = F)


# plot global networks
pando_atac_coef = coef(mome_atac_SCT_plus)
pando_atac_coef = pando_atac_coef[pando_atac_coef$padj<=0.05,]

fff = pando_atac_coef$estimate>0.01 & pando_atac_coef$target %in% pando_atac_coef$tf
dat = as.data.frame(pando_atac_coef[fff,c(1,2,5,9)])

igdat = igraph::graph.data.frame(dat[,1:2])
xxx = igraph::betweenness(igdat, directed = T)
ppp = igraph::page_rank(igdat, directed = T)
ddd = igraph::degree(igdat, mode = "out")
eee = ego_size(igdat, mode = "out", order = 100)
iii = ego_size(igdat, mode = "in", order = 100)
plot(xxx, eee, pch = 21, cex = (ddd+1)/4)
V(igdat)$name[eee==12]

netdat = network::as.network(unique(dat[,1:2]))
set.seed(1)
plot_df = ggnetwork(netdat, layout = "fruchtermanreingold", cell.jitter = 1, component_wise = T)

mk_cluster = lapply(unique(mk_ct_comp[mk_ct_comp$logFC>0.5,"feature"]), function(x) as.character(mk_ct_comp$group[which(mk_ct_comp$logFC==max(mk_ct_comp$logFC[mk_ct_comp$feature==x]))]))
names(mk_cluster) = unique(mk_ct_comp[mk_ct_comp$logFC>0.5,"feature"])
mk_cluster = unique(reshape2::melt(mk_cluster)[,2:1])
colnames(mk_cluster) = c("gene", "cluster")
plot_df = merge(plot_df, mk_cluster, by.x = "vertex.names", by.y = "gene", all.x = T)
plot_df$cluster[is.na(plot_df$cluster)] = "other genes"
plot_df$highlevel = unlist(lapply(strsplit(plot_df$cluster, "_"), function(x) x[1]))
plot_df$highlevel[is.na(plot_df$highlevel)] = "other genes"

labdf = head(dat[order(abs(dat$estimate), decreasing = T),], 1000)
ggplot(plot_df, aes(x = x, y = y, fill = highlevel)) +
  geom_edges(arrow = arrow(length = unit(6, "pt")), colour = "grey35",
             mapping = aes(xend = xend, yend = yend)) +
  geom_nodes() +
  geom_label(data = unique(plot_df[plot_df$vertex.names %in% c(labdf$tf, labdf$target),]), 
             mapping = aes(label = vertex.names), size=2.7, label.padding = unit(0.1, "lines"))+
  theme_blank()+
  theme(legend.position = "bottom")

ggplot(plot_df, aes(x = x, y = y, fill = cluster)) +
  geom_edges(arrow = arrow(length = unit(6, "pt")), colour = "grey35",
             mapping = aes(xend = xend, yend = yend)) +
  guides(fill = guide_legend(nrow = 2))+
  geom_nodes() +
  geom_label(data = unique(plot_df[plot_df$vertex.names %in% c(labdf$tf, labdf$target),]), 
             mapping = aes(label = vertex.names), size=2.7, label.padding = unit(0.1, "lines"))+
  theme_blank()+
  theme(legend.position = "bottom")











