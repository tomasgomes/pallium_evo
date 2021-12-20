#########3
# FUNCTIONS FOR CROSS SPECIES COMPARISONS
#########3



# Calculate transcriptional similarity index using the pairwise genes (not recommended)
calcPWTSI = function(pwde_list, sp_pair, ort_l, ort_scope = "ortholog_one2one", 
                     ct_sets = NULL, fcthr = NULL){
  # which cells to use
  if(is.null(ct_sets)){
    ct_sets = list()
    for(sp in sp_pair){
      ct_sets[[sp]] = unique(unlist(lapply(strsplit(names(pwde_list[[sp]]), "..", fixed = T), 
                                           function(x) c(x[1], x[2]))))
    }
  }
  
  # get DE genes for all cell types being considered
  de_genes_pw = list()
  for(sp in sp_pair){
    comp = names(pwde_list[[sp]])
    de_genes_pw[[sp]] = list()
    for(ct in ct_sets[[sp]]){
      detabl = pwde_list[[sp]]
      if(is.null(fcthr)){
        de_genes_pw[[sp]][[ct]] = unique(unlist(lapply(detabl, 
                                                       function(x) rownames(x[x$cl==ct &
                                                                                x$p_val_adj<=0.05,]))))
      } else{
        de_genes_pw[[sp]][[ct]] = unique(unlist(lapply(detabl, 
                                                       function(x) rownames(x[x$cl==ct &
                                                                                x$p_val_adj<=0.05 &
                                                                                abs(x$avg_log2FC)>fcthr,]))))
      }
    }
  }
  
  # choose species orthologs and scope
  ort_tab = ort_l[[paste0(sp_pair, collapse = ".vs.")]]
  if(!is.null(ort_scope)){
    ort_tab = ort_tab[ort_tab$homology.type==ort_scope,]
  }
  ort_tab = unique(ort_tab[ort_tab[,2] %in% unique(unlist(de_genes_pw[[sp_pair[1]]])) |
                             ort_tab[,4] %in% unique(unlist(de_genes_pw[[sp_pair[2]]])),c(2,4)])
  
  # calculate tSI
  tsi_mat = matrix(NA, nrow = length(names(de_genes_pw[[sp_pair[1]]])),
                   ncol = length(names(de_genes_pw[[sp_pair[2]]])))
  rownames(tsi_mat) = names(de_genes_pw[[sp_pair[1]]])
  colnames(tsi_mat) = names(de_genes_pw[[sp_pair[2]]])
  for(cti in names(de_genes_pw[[sp_pair[1]]])){
    for(ctj in names(de_genes_pw[[sp_pair[2]]])){
      gcti = de_genes_pw[[sp_pair[1]]][[cti]]
      gctj = de_genes_pw[[sp_pair[2]]][[ctj]]
      
      gt = length(intersect(gcti, unique(ort_tab[ort_tab[,2] %in% gctj,1])))
      tSI = 1-sqrt((1-gt/length(gctj))*(1-gt/length(gcti)))
      tsi_mat[cti,ctj] = tSI
    }
  }
  
  return(tsi_mat)
}


# Calculate transcriptional similarity index using marker genes 
## (needs a dataframe including at least a gene and a cluster column)
calcTSI = function(mk_l, sp_pair, ort_l, ort_scope = "ortholog_one2one", ct_sets = NULL){
  # which cells to use
  if(is.null(ct_sets)){
    ct_sets = list()
    for(sp in sp_pair){
      ct_sets[[sp]] = unique(mk_l[[sp]]$cluster)
    }
  }
  
  # get DE genes for all cell types being considered
  de_genes_pw = list()
  for(sp in sp_pair){
    de_genes_pw[[sp]] = list()
    for(ct in ct_sets[[sp]]){
      de_genes_pw[[sp]][[ct]] = mk_l[[sp]]$gene[mk_l[[sp]]$cluster==ct]
    }
  }
  
  # choose species orthologs and scope
  ort_tab = ort_l[[paste0(sp_pair, collapse = ".vs.")]]
  if(!is.null(ort_scope)){
    ort_tab = ort_tab[ort_tab$homology.type==ort_scope,]
  }
  ort_tab = unique(ort_tab[ort_tab[,2] %in% unique(unlist(de_genes_pw[[sp_pair[1]]])) |
                             ort_tab[,4] %in% unique(unlist(de_genes_pw[[sp_pair[2]]])),c(2,4)])
  
  # calculate tSI
  tsi_mat = matrix(NA, nrow = length(names(de_genes_pw[[sp_pair[1]]])),
                   ncol = length(names(de_genes_pw[[sp_pair[2]]])))
  rownames(tsi_mat) = names(de_genes_pw[[sp_pair[1]]])
  colnames(tsi_mat) = names(de_genes_pw[[sp_pair[2]]])
  for(cti in names(de_genes_pw[[sp_pair[1]]])){
    for(ctj in names(de_genes_pw[[sp_pair[2]]])){
      gcti = de_genes_pw[[sp_pair[1]]][[cti]]
      gctj = de_genes_pw[[sp_pair[2]]][[ctj]]
      
      gt = length(intersect(gcti, unique(ort_tab[ort_tab[,2] %in% gctj,1])))
      tSI = 1-sqrt((1-gt/length(gctj))*(1-gt/length(gcti)))
      tsi_mat[cti,ctj] = tSI
    }
  }
  
  return(tsi_mat)
}


# Function to match cells from 2 species and compare
corrCellTypesPW = function(means_list, pwde_list, ort_l, ort_scope = "ortholog_one2one", 
                           sp_pair = c("human", "mouse"), ct_sets = NULL, 
                           gene_filter = NULL, norm_all = T, filter_ct = T){
  # which cells to use
  if(is.null(ct_sets)){
    ct_sets = list(colnames(means_list[[sp_pair[1]]]), colnames(means_list[[sp_pair[2]]]))
    names(ct_sets) = sp_pair
  }
  
  # choose species orthologs and scope
  ort_tab = ort_l[[paste0(sp_pair, collapse = ".vs.")]]
  if(!is.null(ort_scope)){
    ort_tab = ort_tab[ort_tab$homology.type==ort_scope,]
  }
  ort_tab = unique(ort_tab[ort_tab[,2] %in% rownames(means_list[[sp_pair[1]]]) &
                             ort_tab[,4] %in% rownames(means_list[[sp_pair[2]]]),c(2,4)])
  
  # get DE genes for all cell types being considered
  de_genes_pw = list()
  for(spi in 1:length(sp_pair)){
    sp = sp_pair[spi]
    comp = names(pwde_list[[sp]])
    # select those comparing to each other
    compuse = if(length(ct_sets[[sp]])>1){
      rowSums(Reduce(cbind, lapply(ct_sets[[sp]], function(x) grepl(x, comp, fixed = T))))>1
    } else{
      grepl(ct_sets[[sp]], comp, fixed = T)
    }
    de_genes_pw[[spi]] = unique(unlist(lapply(pwde_list[[sp]][compuse], rownames)))
  }
  # using marker intersection
  ort_tab = ort_tab[ort_tab[,1] %in% de_genes_pw[[1]] & ort_tab[,2] %in% de_genes_pw[[2]],]
  
  # only use group of genes
  if(!is.null(gene_filter)){
    ort_tab = ort_tab[ort_tab[,1] %in% gene_filter | ort_tab[,2] %in% gene_filter,]
  }
  
  message(paste0(nrow(ort_tab), " genes used."))
  
  if(filter_ct){ # use only ct present (will print those absent)
    if(any(!(ct_sets[[1]] %in% colnames(means_list[[sp_pair[1]]])))){
      noct = ct_sets[[1]][!(ct_sets[[1]] %in% colnames(means_list[[sp_pair[1]]]))]
      message(paste0(c("Cell types absent in", sp_pair[1], ":", noct), collapse = " "))
    }
    ct_sets[[1]] = ct_sets[[1]][ct_sets[[1]] %in% colnames(means_list[[sp_pair[1]]])]
    
    if(any(!(ct_sets[[2]] %in% colnames(means_list[[sp_pair[2]]])))){
      noct = ct_sets[[2]][!(ct_sets[[2]] %in% colnames(means_list[[sp_pair[2]]]))]
      message(paste0(c("Cell types absent in", sp_pair[2], ":", noct), collapse = " "))
    }
    ct_sets[[2]] = ct_sets[[2]][ct_sets[[2]] %in% colnames(means_list[[sp_pair[2]]])]
  }
  
  # normalise means
  if(norm_all){
    m_list_sub = lapply(means_list[sp_pair], 
                        function(x) t(apply(x, 1, function(y) y/mean(y))))
  } else{
    m_list_sub = list()
    for(spi in 1:length(sp_pair)){
      sp = sp_pair[spi]
      if(length(ct_sets[[spi]])>1){
        m_list_sub[[spi]] = t(apply(means_list[[sp]][,ct_sets[[spi]]], 
                                   1, function(y) y/mean(y)))
      } else{
        m_list_sub[[spi]] = data.frame(means_list[[sp]][,ct_sets[[spi]]])
        colnames(m_list_sub[[spi]]) = ct_sets[[spi]]
      }
      m_list_sub[[spi]][is.na(m_list_sub[[spi]])] = 0
    }
  }
  
  # get correlations, determine max per row and column
  sp1dat = m_list_sub[[1]][ort_tab[,1],ct_sets[[1]],drop=FALSE]
  sp2dat = m_list_sub[[2]][ort_tab[,2],ct_sets[[2]],drop=FALSE]
  if(any(c(dim(sp1dat), dim(sp2dat))==0)) return(NULL)
  cort = psych::corr.test(sp1dat, sp2dat, method = "sp", 
                          adjust = "fdr", alpha = 0.01, ci = F)
  cort$maxrow = apply(cort$r, 1, which.max)
  cort$maxcol = apply(cort$r, 2, which.max)
  cort[[paste0("sp1_", sp_pair[1])]] = sp1dat
  cort[[paste0("sp2_", sp_pair[2])]] = sp2dat
  return(cort)
}


# Function to match cells from 2 species and compare - pass custom genes
corrCellTypes = function(means_list, g_sets = NULL, ort_l, ort_scope = "ortholog_one2one",
                         sp_pair = c("human", "mouse"), ct_sets = NULL, filter_ct = T){
  # which cells to use
  if(is.null(ct_sets)){
    ct_sets = list(colnames(means_list[[sp_pair[1]]]), colnames(means_list[[sp_pair[2]]]))
    names(ct_sets) = sp_pair
  }
  
  # choose species orthologs and scope
  ort_tab = ort_l[[paste0(sp_pair, collapse = ".vs.")]]
  if(!is.null(ort_scope)){
    ort_tab = ort_tab[ort_tab$homology.type==ort_scope,]
  }
  ort_tab = unique(ort_tab[ort_tab[,2] %in% rownames(means_list[[sp_pair[1]]]) &
                             ort_tab[,4] %in% rownames(means_list[[sp_pair[2]]]),c(2,4)])
  
  # get DE genes for all cell types being considered
  if(is.null(g_sets)){
    g_sets = list(rownames(means_list[[sp_pair[1]]]), rownames(means_list[[sp_pair[2]]]))
    names(ct_sets) = sp_pair
  }
  # using marker intersection
  ort_tab = ort_tab[ort_tab[,1] %in% g_sets[[1]] & ort_tab[,2] %in% g_sets[[2]],]
  
  # normalise means
  m_list_sub = lapply(means_list[which(names(means_list) %in% sp_pair)], 
                      function(x) t(apply(x, 1, function(y) y/mean(y))))
  
  # get correlations, determine max per row and column
  if(filter_ct){ # use only ct present (will print those absent)
    if(any(!(ct_sets[[1]] %in% colnames(m_list_sub[[sp_pair[1]]])))){
      noct = ct_sets[[1]][!(ct_sets[[1]] %in% colnames(m_list_sub[[sp_pair[1]]]))]
      message(paste0(c("Cell types absent in", sp_pair[1], ":", noct), collapse = " "))
    }
    ct_sets[[1]] = ct_sets[[1]][ct_sets[[1]] %in% colnames(m_list_sub[[sp_pair[1]]])]
    
    if(any(!(ct_sets[[2]] %in% colnames(m_list_sub[[sp_pair[2]]])))){
      noct = ct_sets[[2]][!(ct_sets[[2]] %in% colnames(m_list_sub[[sp_pair[2]]]))]
      message(paste0(c("Cell types absent in", sp_pair[2], ":", noct), collapse = " "))
    }
    ct_sets[[2]] = ct_sets[[2]][ct_sets[[2]] %in% colnames(m_list_sub[[sp_pair[2]]])]
  }
  
  cort = psych::corr.test(m_list_sub[[sp_pair[1]]][ort_tab[,1],ct_sets[[1]]], 
                          m_list_sub[[sp_pair[2]]][ort_tab[,2],ct_sets[[2]]], method = "sp", 
                          adjust = "fdr", alpha = 0.01, ci = F)
  cort$maxrow = apply(cort$r, 1, which.max)
  cort$maxcol = apply(cort$r, 2, which.max)
  cort[[sp_pair[1]]] = m_list_sub[[sp_pair[1]]][ort_tab[,1],ct_sets[[1]]]
  cort[[sp_pair[2]]] = m_list_sub[[sp_pair[2]]][ort_tab[,2],ct_sets[[2]]]
  return(cort)
}


# Plot correlations
plotCorr = function(cort, sp1 = "sp1", sp2 = "sp2"){
  # cluster and order labels
  hcr = hclust(dist(cort$r), method = "ward.D2")
  hcc = hclust(dist(t(cort$r)), method = "ward.D2")
  hcr = hcr$labels[hcr$order]
  hcc = hcc$labels[hcc$order]
  
  # reshaping the correlations
  plot_df = reshape2::melt(cort$r)
  plot_df$Var1 = factor(plot_df$Var1, levels = rev(hcr))
  plot_df$Var2 = factor(plot_df$Var2, levels = hcc)
  
  # add pvalue and max cor infor
  plot_df$padj = -log10(reshape2::melt(cort$p.adj+min(cort$p.adj[cort$p.adj>0])/10)$value)
  plot_df$rowmax = apply(Reduce(cbind, lapply(names(cort$maxrow), 
                                              function(n) plot_df$Var1==n &
                                                plot_df$Var2==colnames(cort$r)[cort$maxrow[n]])), 
                         1, any)
  plot_df$colmax = apply(Reduce(cbind, lapply(names(cort$maxcol), 
                                              function(n) plot_df$Var2==n &
                                                plot_df$Var1==rownames(cort$r)[cort$maxcol[n]])), 
                         1, any)
  plot_df$markcol = plot_df$value>quantile(plot_df$value, 0.98)
  
  # getting a colourscale where 0 is white in the middle, and intensity leveled by max(abs(value))
  cols = colorRampPalette(c(rev(RColorBrewer::brewer.pal(9, "Blues")),
                            RColorBrewer::brewer.pal(9, "Reds")))(101)
  br = seq(-max(abs(cort$r)), max(abs(cort$r)), length.out = 101)
  cols = cols[!(br>max(cort$r) | br<min(cort$r))]
  
  corplot = ggplot()+
    geom_point(data = plot_df, mapping = aes(x = Var2, y = Var1, fill = value, size = padj), 
               shape = 21)+
    geom_point(data = plot_df[plot_df$rowmax,], mapping = aes(x = Var2, y = Var1, size = padj), 
               shape = "—", show.legend = F, colour = "grey10")+
    geom_point(data = plot_df[plot_df$colmax,], mapping = aes(x = Var2, y = Var1, size = padj), 
               shape = "|", show.legend = F, colour = "grey10")+
    scale_x_discrete(expand = c(0,0.7))+
    scale_y_discrete(expand = c(0,0.7))+
    scale_fill_gradientn(breaks = signif(c(min(cort$r)+0.005, 0, max(cort$r)-0.005),2), 
                         values = scales::rescale(c(min(br), 0, max(br))),
                         colours = cols)+
    labs(x = sp2, y = sp1, fill = "Spearman's\nrho", size = "-log10\nadj. p-value")+
    theme_classic()+
    theme(axis.title = element_text(colour = "black", face = "bold"),
          axis.text = element_text(colour = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          legend.title = element_text(size = 9),
          legend.text = element_text(size = 8))
  
  return(corplot)
}


# Calculate and plot UMAP of cells
umapCellTypes = function(means_list, pwde_list, ort_l, ort_scope = "ortholog_one2one", 
                         ct_sets = list("human" = "all", "mouse" = "all"), seed = 1){ 
  # species are detected from ct_sets names
  # which cells to use
  for(sp in names(ct_sets)){
    if(ct_sets[[sp]][1]=="all"){ 
      ct_sets[[sp]] = colnames(means_list[[sp]])
    }
  }
  
  # get orthologs for all species pairs
  sp_p_mat = combn(names(ct_sets), 2)
  for(i in 1:ncol(sp_p_mat)){
    sp_pair = sp_p_mat[,i]
    ort_tab = ort_l[[paste0(sp_pair, collapse = ".vs.")]]
    if(!is.null(ort_scope)){
      ort_tab = ort_tab[ort_tab$homology.type==ort_scope,]
    }
    # orthologs in the means matrix
    ort_tab = unique(ort_tab[ort_tab[,2] %in% rownames(means_list[[sp_pair[1]]]) &
                               ort_tab[,4] %in% rownames(means_list[[sp_pair[2]]]),c(2,4)])
    colnames(ort_tab) = gsub("genename", "", colnames(ort_tab))
    if(i==1){
      ort_all = ort_tab
    } else{
      ort_all = unique(merge(ort_all, ort_tab, by = intersect(colnames(ort_all), colnames(ort_tab))))
    }
  }
  
  # get DE genes for all cell types being considered
  de_genes_pw = list()
  for(sp in names(ct_sets)){
    comp = names(pwde_list[[sp]])
    # select those comparing to each other
    compuse = rowSums(Reduce(cbind, lapply(ct_sets[[sp]], function(x) grepl(x, comp, fixed = T))))>1
    de_genes_pw[[sp]] = unique(unlist(lapply(pwde_list[[sp]][compuse], rownames)))
  }
  keep = apply(data.frame(lapply(colnames(ort_all), function(x) ort_all[,x] %in% de_genes_pw[[x]])), 
               1, all)
  ort_all = ort_all[keep,]
  
  # normalise means
  m_list_sub = lapply(means_list[which(names(means_list) %in% names(ct_sets))], 
                      function(x) t(apply(x, 1, function(y) y/mean(y))))
  
  avg_exp_all = t(Reduce(cbind, lapply(names(m_list_sub), 
                                       function(x) m_list_sub[[x]][ort_all[,x],ct_sets[[x]]])))
  
  # define metadata
  meta = data.frame("cell_types" = unlist(lapply(m_list_sub, colnames)),
                    "species" = unlist(lapply(names(m_list_sub), 
                                              function(x) rep(x, ncol(m_list_sub[[x]])))))
  rownames(meta) = paste0(meta$species, "_", meta$cell_types)
  meta = meta[match(rownames(avg_exp_all), meta$cell_types),]
  rownames(avg_exp_all) = rownames(meta)
  
  # calculate UMAP
  set.seed(seed)
  l = uwot::umap(avg_exp_all, metric = "cosine", ret_nn = T, n_epochs = 1000)
  l_sp = data.frame(l$embedding)
  rownames(l_sp) = rownames(avg_exp_all)
  
  plot_df = cbind(l_sp, meta[rownames(l_sp),])
  
  return(plot_df)
}


# Ortholog match - get orthologs for all species pairs
ortMatch = function(ort_l, sp_vec, ort_scope = "ortholog_one2one"){
  sp_p_mat = combn(sp_vec, 2)
  for(i in 1:ncol(sp_p_mat)){
    sp_pair = sp_p_mat[,i]
    ort_tab = ort_l[[paste0(sp_pair, collapse = ".vs.")]]
    if(!is.null(ort_scope)){
      ort_tab = ort_tab[ort_tab$homology.type==ort_scope,]
    }
    
    colnames(ort_tab) = gsub("genename", "", colnames(ort_tab))
    if(i==1){
      ort_all = unique(ort_tab[,c(2,4)])
    } else{
      ort_all = unique(merge(ort_all, unique(ort_tab[,c(2,4)]), 
                             by = intersect(colnames(ort_all), colnames(ort_tab))))
    }
  }
  
  return(ort_all)
}


# Subsetting Seurat objects using the same orthologs
seuratOrthologs = function(s_l, sp_samples, o_all){
  # only orthologs present in the Seurat objects
  for(sp in names(s_l)){
    o_all = o_all[o_all[,sp_samples[sp]] %in% rownames(s_l[[sp]]),]
  }
  
  # get a unique name for the genes
  allg = unlist(o_all)
  jointnames = c()
  rep_n = list()
  for(i in 1:nrow(o_all)){
    nname = unlist(lapply(o_all[i,], function(x) sum(allg==x)))
    g = unlist(o_all[i,names(nname[nname==1])])
    # if there are acceptable unique names
    if(any(nname==1) & (!is.null(g) & !all(grepl("..", g, fixed = T)))){ # no weird axololt annot
      g = g[!grepl("..", g, fixed = T)] # no weird axololt annot
      g = c(g[!grepl("LOC", g)], g[grepl("LOC", g)])[1] #prioritize non-LOC
      jointnames = c(jointnames, g)
      rep_n[[g]] = 1
    } else{ # we'll add a number to those names
      g = unlist(o_all[i,])
      g = g[!grepl("..", g, fixed = T)] # no weird axololt annot
      g = c(g[!grepl("LOC", g)], g[grepl("LOC", g)])[1] #prioritize non-LOC
      if(g %in% names(rep_n)){
        gn = rep_n[[g]]
        rep_n[[g]] = gn+1
        jointnames = c(jointnames, paste0(g, "-", gn+1))
      } else{
        rep_n[[g]] = 1
        jointnames = c(jointnames, paste0(g, "-", 1))
      }
    }
  }
  
  # filter the Seurats by the orthologs
  ## this requires redoing the Seurat with new names, since some genes will be duplicated
  for(sp in names(s_l)){
    cnts = s_l[[sp]]@assays$RNA@counts[o_all[,sp_samples[sp]],]
    rownames(cnts) = jointnames
    m = s_l[[sp]]@meta.data
    s_l[[sp]] = CreateSeuratObject(cnts, project = sp, meta.data = m)
  }
  
  return(s_l)
}