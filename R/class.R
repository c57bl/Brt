#' @import dplyr
NULL
# ================================ brtVirus ====================================

#' the brtvirus class
#' @slot bc a data.frame describe bc and counts of each bc
#' @slot batch which batch the virus is
#' @slot U1000 the unique proportion of virus, if 1000 virus particles were sampled
setClass("brtVirus",
         slots = c(bc = "data.frame",
                   batch = "character",
                   bc_type = "data.frame",
                   U1000 = "numeric",
                   tools = "list"))
#' brtVirus
#' construct brtVirus object
#' @param x a data.frame, with colnames bc describe the sequence of each barcode
#' and counts describe the counts of each bc
#' @param batch which batch the virus is
brtVirus <- function(x,batch){
  1-(purrr::map(x$counts/sum(x$counts),~.x*(1-(1-.x)^(1000-1))) %>%
       unlist %>%
       sum) -> u1000
  x %>%
    select(bc,counts) %>%
    arrange(desc(counts)) -> x
  purrr::map(seq(0,100000,100),function(n){
    x %>%
      mutate(P = counts/sum(counts)) %>%
      mutate(PN = 1-(1-P)^n) %>%
      select(PN) %>%
      sum()
  }) %>%
    unlist() -> PN
  bc_type <- data.frame(S = seq(0,100000,100), PN = PN)
  new(
    "brtVirus",
    bc = x,
    batch = batch,
    bc_type = bc_type,
    U1000 = u1000,
    tools = list()
  )
}

setMethod("show","brtVirus",function(object){
  cat("virus ",object@batch,":",'\n',"barcodes number: ",nrow(object@bc),'\n',
      "U1000: ",object@U1000)
})
#' The brtData class
#' @slot data a list of Matrix
#' @slot focus the default data
setClass("brtData",
         slots = c(data = 'list',
                   focus = 'character'))

brtData<-function(counts = NA){
  new("brtData",
      data = list(counts = Matrix::Matrix(counts)),
      focus = "counts")
}

#' The brtExperiment class
#' @slot data brtData class
#' @slot coldata data.frame describe the features of each column-element
#' @slot rowdata data.frame describe the features of each row-element
setClass("brtExperiment",
         slots = c(data="brtData",
                   coldata = 'data.frame',
                   rowdata = 'data.frame',
                   virus = "brtVirus",
                   tools = 'list'),
         prototype = list(data = brtData(),
                          coldata = data.frame(),
                          rowdata = data.frame(),
                          tools = list()))
setMethod("show","brtExperiment",function(object){
  cat(nrow(object@rowdata),'features across',nrow(object@coldata),'samples\n')
})

# ==============================brtStarter======================================

#' The brtStarter class
#' @slot rowdata a data.frame, with conserved colnames name,bc,rvcounts.raw
#' @slot coldata a data.frame, with conserved colnames name,sample,region,cell,
#' rvcounts.raw
setClass("brtStarter",
         contains = "brtExperiment")

#' brtStarter
#' @param x an tibble returned from brtMatch_10x_bc
#' @param seuratobj the seurat object of this sample
#' @param sample sample name
#' @param source from which region
#' @return an brtStarter obj
#' @export
brtStarter <- function(x,seuratobj,sample,region,verbose = TRUE) {
  # import cell name from seurat obj
  seurat.metadata <- brtSeuratWraper(seuratobj)
  x <- subset(x, cell %in% seurat.metadata$cell)
  x <- select(x,cell,bc,counts)
  # construct sparse matrix
  x.spread <- tidyr::spread(x, cell, counts)
  x.matrix <- Matrix::Matrix(as.matrix(x.spread[, -1]),
                             sparse = TRUE)
  x.matrix[is.na(x.matrix)] <- 0
  rownames(x.matrix) <-  x.spread$bc
  # sort matrix
  x.matrix <- .sortMatrix(x.matrix, "col")
  x.matrix <- .sortMatrix(x.matrix, 'row')
  # build brtdata from counts
  x.data <- brtData(counts = x.matrix)
  # initialize coldata and rowdata
  cell <- colnames(x.matrix)
  x.coldata <-
    data.frame(
      name = paste(sample,region,cell,sep = "_"),
      sample = sample,
      region = region,
      cell = cell,
      rvcounts.raw = Matrix::colSums(x.matrix)
    )
  bc <- rownames(x.matrix)
  x.rowdata <- data.frame(
    name = bc,
    bc = bc,
    rvcounts.raw = Matrix::rowSums(x.matrix)
  )
  # update seurat metadata
  x.coldata <- left_update(x.coldata,
                           seurat.metadata,
                           by = "cell",
                           treat.na = NA)
  # record conserve colname of col/row data
  rawnames <- list(rowdata = colnames(x.rowdata),
                   coldata = colnames(x.coldata))
  new(
    "brtStarter",
    data = x.data,
    rowdata = x.rowdata,
    coldata = x.coldata,
    tools = list(rawnames = rawnames)
  )
}





