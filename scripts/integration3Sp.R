library(Seurat)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(harmony)


glut_ort_l = readRDS("data/processed/integration_evo/glut_ort_3.RDS")
names(glut_ort_l) = paste0("glut_", names(glut_ort_l), "_3")
gaba_ort_l = readRDS("data/processed/integration_evo/gaba_ort_3.RDS")
names(gaba_ort_l) = paste0("gaba_", names(gaba_ort_l), "_3")

glut_gaba = c(glut_ort_l, gaba_ort_l)

for(n in names(glut_gaba)){
  print(n)
  # merge data
  srat_m = Reduce(merge, glut_gaba[[n]])
  srat_m = AddMetaData(srat_m, metadata = unlist(lapply(names(glut_gaba[[n]]), 
                                                        function(x) rep(x, ncol(glut_gaba[[n]][[x]])))), col.name = "species")
  # select integration features
  #integr_feat = rownames(srat_m@assays$SCT@scale.data)
  integr_feat = Reduce(intersect, lapply(glut_gaba[[n]], VariableFeatures))
  
  ##### rPCA #####
  print("rPCA")
  # calculate possible missing Pearson residuals for SCTransform
  srat_ort_int = PrepSCTIntegration(glut_gaba[[n]], anchor.features = integr_feat, verbose = T)
  
  srat_ort_int <- lapply(X = srat_ort_int, FUN = function(x) {
    #x <- ScaleData(x, features = integr_feat, verbose = FALSE)
    x <- RunPCA(x, features = integr_feat, verbose = FALSE)
  })
  
  # finding the anchors for integration
  all_cell_anchors = FindIntegrationAnchors(srat_ort_int, normalization.method = "SCT", 
                                            anchor.features = integr_feat, reduction = "rpca",
                                            assay = rep("SCT", length(srat_ort_int)), 
                                            k.anchor = 25, dims = 1:15, n.trees = 200,
                                            max.features = 500, verbose = T)
  
  # actual data integration
  allgenes = rownames(srat_ort_int[[1]]@assays$RNA@counts)
  srat_ort_rpca = IntegrateData(all_cell_anchors, normalization.method = "SCT", dims = 1:15,
                                verbose = T, features.to.integrate = allgenes,
                                sample.tree = matrix(c(-1,-2,1,-3), 
                                                     nrow = 2, ncol = 2, byrow = T))
  
  # run PCA and UMAP to see how it looks
  srat_ort_rpca = RunPCA(srat_ort_rpca, verbose = F)
  srat_ort_rpca = RunUMAP(srat_ort_rpca, dims = 1:15)
  
  srat_ort_rpca = AddMetaData(srat_ort_rpca,  metadata = unlist(lapply(names(glut_gaba[[n]]), 
                                                                       function(x) rep(x, ncol(glut_gaba[[n]][[x]])))), col.name = "species")
  
  # save
  saveRDS(srat_ort_rpca, file = paste0("data/processed/integration_evo/", n, "_srat_ort_rpca.RDS"))
  
  
  ##### Harmony #####
  print("Harmony")
  VariableFeatures(srat_m) = integr_feat
  #srat_ort_harm = ScaleData(srat_m, features = integr_feat, use.umi = T, 
  #                          do.scale = F, verbose = FALSE)
  srat_ort_harm = RunPCA(srat_m, verbose = FALSE, assay = "SCT", npcs = 50,
                         features = integr_feat)
  
  # Run Harmony
  srat_ort_harm = RunHarmony(srat_ort_harm, "species", tau = 30, 
                             plot_convergence = F, assay.use = "SCT")
  
  # Run UMAP on Harmony dimensions
  srat_ort_harm = RunUMAP(srat_ort_harm, reduction = "harmony", dims = 1:15)
  
  # save
  saveRDS(srat_ort_harm, file = paste0("data/processed/integration_evo/", n, "_srat_ort_harm.RDS"))
  
  
  ##### Combined #####
  print("Combined")
  # Run Harmony
  srat_ort_both = RunHarmony(srat_ort_rpca, "species", tau = 30, 
                             reduction = "pca", dims.use = 1:15,
                             plot_convergence = F, assay.use = "SCT")
  
  # Run UMAP on Harmony dimensions
  srat_ort_both = RunUMAP(srat_ort_both, reduction = "harmony", dims = 1:15)
  
  # save
  saveRDS(srat_ort_both, file = paste0("data/processed/integration_evo/", n, "_srat_ort_both.RDS"))
}
