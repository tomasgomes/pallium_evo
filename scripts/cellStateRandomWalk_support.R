# libraries
library(igraph)
library(Matrix)


rowNorm <- function(mat){
  diag_mat <- Diagonal(x = 1/rowSums(mat))
  res <- diag_mat %*% mat
  dimnames(res) <- dimnames(mat)
  return(res)
}

colNorm <- function(mat){
  diag_mat <- Diagonal(x = 1/colSums(mat))
  res <- mat %*% diag_mat
  dimnames(res) <- dimnames(mat)
  return(res)
}

# util functions
summarize_data_to_groups.Matrix <- function(object, groups, bycol = T)
{
  groups <- factor(groups)
  summfunc <- Matrix::rowMeans
  if (bycol) summfunc <- Matrix::colMeans
  
  summ_mat <- sapply(levels(groups), function(x){
    if (bycol) obj <- object[which(groups==x),]
    if (!bycol) obj <- object[,which(groups==x)]
    summfunc(obj)
  })
  if (! bycol) summ_mat <- t(summ_mat)
  return(summ_mat)
}

summarize_data_to_groups.data.frame <- function(object, groups)
{
  groups <- factor(groups)
  ident_mat_groups <- sparseMatrix(i = 1:length(groups),
                                   j = as.numeric(groups),
                                   x = 1,
                                   dims = c(length(groups), length(levels(groups))), dimnames = list(names(groups), levels(groups)))
  ident_mat_groups_norm <- colNorm(ident_mat_groups)
  
  df_summ <- do.call(cbind.data.frame, lapply(1:ncol(object), function(i){
    dat <- object[,i]
    if (is.numeric(dat)){
      dat <- matrix(dat, nrow = 1, dimnames = list(colnames(object)[i], rownames(object)))
      summ_dat <- dat %*% ident_mat_groups_norm
      return(setNames(as.numeric(summ_dat), colnames(ident_mat_groups_norm)))
    } else{
      dat <- factor(dat)
      lab_mat_groups <- sparseMatrix(i = which(!is.na(dat)),
                                     j = as.numeric(dat)[which(!is.na(dat))],
                                     x = 1,
                                     dims = c(length(dat), length(levels(dat))), dimnames = list(rownames(object), levels(dat)))
      summ_mat <- t(lab_mat_groups) %*% ident_mat_groups
      summ_dat <- factor(setNames(rownames(summ_mat)[apply(summ_mat, 2, which.max)], levels(groups)), levels = levels(dat))
      return(summ_dat)
    }
  }))
  colnames(df_summ) <- colnames(object)
  
  return(df_summ)
}

summarize_data_to_groups <- function(object, ...) {
  UseMethod(generic = 'summarize_data_to_groups', object = object)
}

data_to_h5ad.Seurat <- function(object,
                                assay = DefaultAssay(object),
                                savefile = NULL,
                                verbose = F)
{
  if (verbose)
    cat(">> start to create the anndata object...\n")
  adata <- anndata::AnnData(X = t(object[[assay]]@data),
                            obs = object@meta.data,
                            var = object[[assay]]@meta.features,
                            layers = list(count = t(object[[assay]]@counts)),
                            obsm = setNames(lapply(names(object@reductions), function(x) Embeddings(object, x)), names(object@reductions))
  )
  if (verbose)
    cat(">> done.\n")
  
  if (!is.null(savefile)){
    if (verbose)
      cat(paste0(">> saving the anndata object to file: ",savefile,"\n"))
    adata$write_h5ad(savefile)
    if (verbose)
      cat(">> done.\n")
  }
  
  return(adata)
}

data_to_h5ad.default <- function(object,
                                 vars = NULL,
                                 obs = NULL,
                                 obsm = list(),
                                 layers = list(),
                                 savefile = NULL,
                                 verbose=F)
{
  if (verbose)
    cat(">> start to create the anndata object...\n")
  if (is.null(vars))
    vars <- data.frame(id = colnames(object), row.names = colnames(object))
  if (is.null(obs))
    obs <- data.frame(id = rownames(object), row.names = rownames(object))
  adata <- anndata::AnnData(X = object,
                            obs = obs,
                            var = vars,
                            obsm = obsm,
                            layers = layers
  )
  if (verbose)
    cat(">> done.\n")
  
  if (!is.null(savefile)){
    if (verbose)
      cat(paste0(">> saving the anndata object to file: ",savefile,"\n"))
    adata$write_h5ad(savefile)
    if (verbose)
      cat(">> done.\n")
  }
  
  return(adata)
}

data_to_h5ad <- function(object, ...) {
  UseMethod(generic = 'data_to_h5ad', object = object)
}