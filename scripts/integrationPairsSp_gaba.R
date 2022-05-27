library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(harmony)


gaba_ort_l = readRDS("data/processed/integration_evo/gaba_ort_l.RDS")

for(n in names(gaba_ort_l)){
  print(n)
  # merge data
  srat_m = Reduce(merge, gaba_ort_l[[n]])
  srat_m = AddMetaData(srat_m, metadata = unlist(lapply(names(gaba_ort_l[[n]]), 
                                                        function(x) rep(x, ncol(gaba_ort_l[[n]][[x]])))), col.name = "species")
  integr_feat = rownames(srat_m@assays$SCT@scale.data)
  
  ##### rPCA #####
  print("rPCA")
  # calculate possible missing Pearson residuals for SCTransform
  srat_ort_int = PrepSCTIntegration(gaba_ort_l[[n]], anchor.features = integr_feat, verbose = T)
  
  srat_ort_int <- lapply(X = srat_ort_int, FUN = function(x) {
    x <- ScaleData(x, features = integr_feat, verbose = FALSE)
    x <- RunPCA(x, features = integr_feat, verbose = FALSE)
  })
  
  # finding the anchors for integration
  all_cell_anchors = FindIntegrationAnchors(srat_ort_int, normalization.method = "SCT", 
                                            anchor.features = integr_feat, reduction = "rpca",
                                            assay = rep("SCT", length(srat_ort_int)), 
                                            dims = 1:50, verbose = T)
  
  # actual data integration
  allgenes = rownames(srat_ort_int[[1]]@assays$RNA@counts)
  srat_ort_rpca = IntegrateData(all_cell_anchors, normalization.method = "SCT", dims = 1:50,
                                verbose = T, features.to.integrate = allgenes)
  
  # run PCA and UMAP to see how it looks
  srat_ort_rpca = RunPCA(srat_ort_rpca, verbose = F)
  srat_ort_rpca = RunUMAP(srat_ort_rpca, dims = 1:30)
  
  srat_ort_rpca = AddMetaData(srat_ort_rpca,  metadata = unlist(lapply(names(gaba_ort_l[[n]]), 
                                                                       function(x) rep(x, ncol(gaba_ort_l[[n]][[x]])))), col.name = "species")
  
  # save
  saveRDS(srat_ort_rpca, file = paste0("data/processed/integration_evo/", n, "_gaba_srat_ort_rpca.RDS"))
  
  
  ##### Harmony #####
  print("Harmony")
  VariableFeatures(srat_m) = integr_feat
  srat_ort_harm = ScaleData(srat_m, features = integr_feat, use.umi = T, 
                            do.scale = F, verbose = FALSE)
  srat_ort_harm = RunPCA(srat_ort_harm, verbose = FALSE, assay = "SCT", npcs = 50,
                         features = integr_feat)
  
  # Run Harmony
  srat_ort_harm = RunHarmony(srat_ort_harm, "species", tau = 30, 
                             plot_convergence = F, assay.use = "SCT")
  
  # Run UMAP on Harmony dimensions
  srat_ort_harm = RunUMAP(srat_ort_harm, reduction = "harmony", dims = 1:10)
  
  # save
  saveRDS(srat_ort_harm, file = paste0("data/processed/integration_evo/", n, "_gaba_srat_ort_harm.RDS"))
  
  
  ##### Combined #####
  print("Combined")
  # Run Harmony
  srat_ort_both = RunHarmony(srat_ort_rpca, "species", tau = 30, 
                             plot_convergence = F, assay.use = "SCT")
  
  # Run UMAP on Harmony dimensions
  srat_ort_both = RunUMAP(srat_ort_both, reduction = "harmony", dims = 1:10)
  
  # save
  saveRDS(srat_ort_both, file = paste0("data/processed/integration_evo/", n, "_gaba_srat_ort_both.RDS"))
}
