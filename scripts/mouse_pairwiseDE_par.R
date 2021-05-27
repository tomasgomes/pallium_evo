library(Seurat)
library(parallel)
numCores = 20

pairwiseDESeurat = function(i, mat, ds){
  nn = paste0(mat[1,i],"..",mat[2,i])
  mk= FindMarkers(ds, ident.1 = mat[1,i], ident.2 = mat[2,i],logfc.threshold = 0.2, only.pos = F, 
                  pseudocount.use = 0.1, assay = "SCT")
  mk = mk[mk$p_val_adj<=0.05,]
  mk$cl = ifelse(mk$avg_log2FC>0, mat[1,i], mat[2,i])
  
  return(mk)
}

dataset = readRDS("data/expression/l5_all_seurat.RDS")
meta = read.csv("data/annotations/mouse_all_umeta.csv", header = T, row.names = 1)

dataset = AddMetaData(dataset, metadata = meta)

# normalisation
dataset = suppressWarnings(SCTransform(dataset, do.correct.umi = T, verbose = F, 
                                       seed.use = 1, vars.to.regress = "nCount_RNA",
                                       variable.features.rv.th = 1, return.only.var.genes = F,
                                       variable.features.n = NULL))

# pairwise DE
avg_exp = list()
mk_list = list()
for(g in c("subclasses", "cellclusters")[1]){
  dataset = SetIdent(dataset, value = g)
  avg_exp[[g]] = AverageExpression(dataset, assays = "SCT")$SCT
  
  freqcl = table(dataset@meta.data[,g])
  ucl = unique(dataset@meta.data[,g])
  clpairs = combn(ucl[ucl %in% names(freqcl)[freqcl>=5]], 2)
  avg_exp[[g]] = avg_exp[[g]][,colnames(avg_exp[[g]]) %in% names(freqcl)[freqcl>=5]]
  
  mk_list[[g]] = mclapply(1:ncol(clpairs), pairwiseDESeurat, 
                          mat = clpairs, ds = dataset, mc.cores = numCores)
  names(mk_list[[g]]) = apply(clpairs, 2, function(x) paste0(x[1], "..", x[2]))
}

# save results
saveRDS(mk_list, file = "results/pairwiseDE/mouse_markers.RDS")
saveRDS(avg_exp, file = "results/pairwiseDE/mouse_avgExp.RDS")


