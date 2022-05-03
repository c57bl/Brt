# =============================== wrap seurat ==================================
#' brtSeuratWraper
#' @description wrap seurat object
#' @param seuratobj seurat objects
#' @param rduction_method 'umap' or 'tsne'
#' @examples
#' foo<-brtSeuratWrap(brt1_seurat,level=1)
#' @return a metadata of seurat
brtSeuratWraper <- function(seuratobj, reduction_method = 'umap') {
  if (("Seurat" %in% is(seuratobj)) == FALSE) {
    stop("input class error")
  }
  if (isFALSE(reduction_method %in% c('umap', 'tsne'))) {
    stop('reduction_method should be "umap" or "tsne"')
  }

  if (reduction_method == 'umap') {
    reduction <- seuratobj@reductions$umap@cell.embeddings
  } else {
    reduction <- seuratobj@reductions$tsne@cell.embeddings
  }
  reduction <- data.frame(reduction)
  reduction$cell <- substr(rownames(reduction), 1, 16)
  ident <- data.frame(seuratobj@active.ident)
  ident$cell <- substr(rownames(ident), 1, 16)
  ident <- rename(ident, ident = seuratobj.active.ident)
  metadata <- seuratobj@meta.data[, 1:4]
  metadata$cell <- substr(rownames(metadata), 1, 16)
  metadata <- left_join(metadata, ident, by = 'cell')
  metadata <- left_join(metadata, reduction, by = 'cell')
  metadata <- cbind(metadata[, 5], metadata[,-5])
  colnames(metadata)[1] <- 'cell'
  return(metadata)
}

