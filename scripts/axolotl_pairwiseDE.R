library(Seurat)

dataset = readRDS("data/expression/axolotl/palliumres07.RDS")
meta = read.csv("data/annotations/axolotl_all_umeta.csv", header = T, row.names = 1)

dataset = AddMetaData(dataset, metadata = meta)

# normalisation
dataset = suppressWarnings(SCTransform(dataset, do.correct.umi = T, verbose = F, 
                                       seed.use = 1, vars.to.regress = "nCount_RNA",
                                       variable.features.rv.th = 1, return.only.var.genes = F,
                                       variable.features.n = NULL))

# pairwise DE
avg_exp = list()
mk_list = list()
for(g in c("subclasses", "cellclusters")){
  dataset = SetIdent(dataset, value = g)
  avg_exp[[g]] = AverageExpression(dataset, assays = "SCT")$SCT
  
  freqcl = table(dataset@meta.data[,g])
  ucl = unique(dataset@meta.data[,g])
  clpairs = combn(ucl[ucl %in% names(freqcl)[freqcl>=5]], 2)
  avg_exp[[g]] = avg_exp[[g]][,colnames(avg_exp[[g]]) %in% names(freqcl)[freqcl>=5]]
  
  mk_list[[g]] = list()
  for(i in 1:ncol(clpairs)){
    nn = paste0(clpairs[1,i],"..",clpairs[2,i])
    mk_list[[g]][[nn]] = FindMarkers(dataset, ident.1 = clpairs[1,i], ident.2 = clpairs[2,i],
                                     logfc.threshold = 0.2, only.pos = F, pseudocount.use = 0.1, 
                                     assay = "SCT")
    mk_list[[g]][[nn]] = mk_list[[g]][[nn]][mk_list[[g]][[nn]]$p_val_adj<=0.05,]
    mk_list[[g]][[nn]]$cl = ifelse(mk_list[[g]][[nn]]$avg_log2FC>0, clpairs[1,i], clpairs[2,i])
  }
}

# save results
saveRDS(mk_list, file = "results/pairwiseDE/axolotl_markers.RDS")
saveRDS(avg_exp, file = "results/pairwiseDE/axolotl_avgExp.RDS")


