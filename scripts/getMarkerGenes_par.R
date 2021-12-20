library(Seurat)
library(future)
library(parallel)
numCores = 32
options(future.globals.maxSize= 2111289600)

DESeurat = function(i, ds, allct){
  plan("multiprocess", workers = 30)
  markers = FindMarkers(ds, ident.1 = allct[i], verbose = FALSE, pseudocount.use = 0.1, only.pos = T)
  return(markers)
}



# files lists
expr_f = list(axolotl = "data/expression/axolotl_nuclei_filtered_annotated.RDS",
              human = "data/expression/human_10x.RDS",
              mouse = "data/expression/l5_all_seurat.RDS",
              turtle = "data/expression/turtle_all_v3.RDS",
              lizard = "data/expression/lizard_all_v3.RDS",
              zebrafinch = "data/expression/zfinch_seurat.RDS",
              bengalesefinch = "data/expression/bfinch_seurat.RDS",
              zebrafish = "data/expression/drerio_brain_v3.RDS")[2:3]
meta_f = list(axolotl = "data/annotations/axolotl_nuc_umeta.csv",
              human = "data/annotations/human10x_all_umeta.csv",
              mouse = "data/annotations/mouse_all_umeta.csv",
              turtle = "data/annotations/turtle_all_umeta.csv",
              lizard = "data/annotations/lizard_all_umeta.csv",
              zebrafinch = "data/annotations/zfinch_all_umeta.csv",
              bengalesefinch = "data/annotations/bfinch_all_umeta.csv",
              zebrafish = "data/annotations/drerio_all_umeta.csv")[2:3]

mk_list = list()
for(sp in names(expr_f)){ #for each species
  print(sp)
  # load and format dataset
  dataset = readRDS(expr_f[[sp]])
  meta = read.csv(meta_f[[sp]], header = T, row.names = 1)
  dataset = AddMetaData(dataset, metadata = meta)
  dataset = dataset[,dataset@meta.data$cellclusters!="doublets" &
                      dataset@meta.data$subclasses!="doublets"]
  
  # normalisation
  dataset = suppressWarnings(SCTransform(dataset, do.correct.umi = T, verbose = F, 
                                         seed.use = 1, vars.to.regress = "nCount_RNA",
                                         variable.features.rv.th = 1, return.only.var.genes = F,
                                         variable.features.n = NULL))
  
  mk_list[[sp]] = list()
  for(g in c("subclasses", "cellclusters")){ # for each relevant column
    print(g)
    dataset = SetIdent(dataset, value = g)
    
    # only clusters with 5 cell minimum
    freqcl = table(dataset@meta.data[,g])
    ucl = unique(dataset@meta.data[,g])
    ucl = ucl[ucl %in% names(freqcl)[freqcl>=5]]
    
    # get markers for each cluster
    #mk_list[[sp]][[g]] = mclapply(1:length(ucl), DESeurat, ds = dataset, allct = ucl, mc.cores = numCores)
    mk_list[[sp]][[g]] = lapply(1:length(ucl), DESeurat, ds = dataset, allct = ucl)
    mk_list[[sp]][[g]] = lapply(mk_list[[sp]][[g]], function(x) if(is.data.frame(x)){x[x$p_val_adj<=0.05,]})
    names(mk_list[[sp]][[g]]) = ucl
  }
  
  # save results
  saveRDS(mk_list[[sp]], file = paste0("results/pairwiseDE/markersGenes_", sp, ".RDS"))
}

# save results
#saveRDS(mk_list, file = "results/pairwiseDE/markersGenes_allspecies.RDS")


