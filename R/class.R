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
#' @slot rowdata a data.frame, with conserved colnames bc,rvcounts.raw
#' @slot coldata a data.frame, with conserved colnames name,sample,region,cell,
#' rvcounts.raw
setClass("brtStarter",
         contains = "brtExperiment",
         slots = c(virus = "brtVirus"))

#' brtStarter
#' @param x an tibble returned from brtMatch_10x_bc
#' @param seuratobj the seurat object of this sample
#' @param virus a brtVirus obj
#' @param sample sample name
#' @param source from which region
#' @return an brtStarter obj
#' @export
brtStarter <- function(x,seuratobj,virus,sample,region,verbose = TRUE) {
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
    bc = bc,
    rvcounts.raw = Matrix::rowSums(x.matrix)
  )
  # update seurat metadata
  x.coldata <- left_update(x.coldata,
                           seurat.metadata,
                           by = "cell",
                           treat.na = NA)
  # build brtdata from counts
  colnames(x.matrix) <- x.coldata$name
  x.data <- brtData(counts = x.matrix)

  new(
    "brtStarter",
    data = x.data,
    rowdata = x.rowdata,
    coldata = x.coldata,
    virus = virus,
    tools = list()
  )
}


# ============================ brtUnity ========================================
#' The brtUnity class
#' @slot metadata a list contains coldata of all the coldatas of input objs.
#' support one starter region in one animal only.
#' @slot sc a list of objs records single-cell connection
#' @slot bulk a brtInputome obj record bulk-regions connection
#' @slot tools a list record some useful objs
setClass("brtUnity",
         slots = c(metadata = "list",
                   sc = 'list',
                   bulk = 'brtExperiment',
                   tools = 'list'))
#' brtUnity
#' @param starter a list of starter objs
#' @param sc a list of brtStarter objs records single cell long-range inputs
#' @param bulk a list of brtInputome objs records long-range bulk inputs
#' @return  a brtUnity obj
#' @seealso [brtUnity-class]
#' @importFrom Matrix Matrix
#' @export
brtUnity <- function(starter = list, sc = NULL, bulk = NULL){
  if (length(starter) == 0) stop("starter should not be NULL")
  samples <- purrr::map(starter,~paste0(unique(.x@coldata$sample))) %>%
    unlist
  bcmaps <- purrr::map(starter,function(x){
    filter(x@rowdata,unique.bc == TRUE) %>%
      select(bc,cell) %>%
      arrange(cell)
  })
  names(bcmaps) <- samples
  names(starter) <- samples
  # merge starter
  purrr::map2(starter,bcmaps,function(x,y){
    input.idx <- which(x@coldata$state == "input")
    bc.idx <- which(x@rowdata$unique.bc == TRUE)
    data.sub <- purrr::map(x@data@data,~.x[bc.idx,input.idx])
    counts.sub <- data.sub$counts
    counts.merge <- mergeBc(counts.sub,y)
    newObjFromMerged(counts.merge,x)
  }) -> scobjs
  names(scobjs) <- samples
  # a littl complicated here
  # merge sc input, merge was performed sample by sample. the sc inputs with
  # the same sample will be merged into the starter based sc obj in column
  # direction one by one.
  if (! is.null(sc)){
    # build idx first
    purrr::map(sc,~unique(.x@coldata$sample)) %>%
      unlist -> samples.sc
    purrr::map(sc,~unique(paste(.x@coldata$sample,.x@coldata$region,sep = "_"))) %>%
      unlist -> id.sc
    names(sc) <- id.sc
    # enter sample by sample loop
    scobjs <- purrr::map2(scobjs,names(scobjs),function(x,y){
      scidx <- which(samples.sc == y)
      # check all sc inputs and search the matched sample name
      if (length(scidx) > 0){
        for (i in 1:length(scidx)){
          sc.curr <- sc[[scidx[i]]]
          sc.curr <- subsetCol(sc.curr,state == "input")
          counts.merge <- mergeBc(sc.curr@data@data$counts,bcmaps[[y]])
          newobj <- newObjFromMerged(counts.merge,sc.curr)
          # keep updating the matched scobj
          x@data@data$counts <- cbind(x@data@data$counts,
                                      newobj@data@data$counts)
          x@coldata <- bind_rows(x@coldata,newobj@coldata)
        }
      }
      x
    })
  }
  if (is.null(bulk)) bulkobjs <- new("brtExperiment")
  # record all the coldata
  metadata <- purrr::map(list(starter = starter,sc = sc,bulk = bulk),function(x){
    purrr::map(x,~.x@coldata)
    })
  new("brtUnity",metadata = metadata,sc = scobjs, bulk = bulkobjs)
}

setMethod("show","brtUnity",function(object){
  cat("a brtUnity obj with ", length(object@sc)," samples","\n\n")
  cat("single cell connections: \n")
  show(object@sc)
  cat("bulk connections: \n")
  show(object@bulk)
})

# inner function of brtUnity, merge bc in same cells
mergeBc <- function(data,bcmap) {
  cellList <- unique(bcmap$cell)
  newdata <- purrr::map(cellList, function(starterid) {
    bclist <- subset(bcmap, cell == starterid)$bc
    bcidx <-
      purrr::map(bclist,  ~ which(.x == rownames(data))) %>%
      unlist
    if (length(bcidx) == 0) {
      newdata <- Matrix(rep(0, ncol(data)), nrow = 1)
    } else if (length(bcidx) == 1) {
      newdata <- data[bcidx, ]
    } else if (length(bcidx) > 1) {
      newdata <- colSums(data[bcidx, ])
    }
    return(newdata)
  })
  newdata <- do.call(rbind, newdata)
  rownames(newdata) <- cellList
  colnames(newdata) <- colnames(data)
  return(newdata)
}

# innner function of brtunity, construct new brtExperiment using merged bc
newObjFromMerged <- function(counts,obj){
  counts <- counts[,which(colSums(counts) >0)]
  data <- brtData(counts = counts)
  rowname.idx <- match(rownames(counts), obj@coldata$name)
  rowdata <- obj@coldata[rowname.idx,]
  colname.idx <- match(colnames(counts), obj@coldata$name)
  coldata <- obj@coldata[colname.idx,]
  new("brtExperiment",
      data = data,
      rowdata = rowdata,
      coldata = coldata)
}
