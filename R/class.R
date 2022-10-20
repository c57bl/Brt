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

brtSetIdent <- function(obj,ident) {
  if (! ident %in% colnames(obj@coldata))
    stop("can not find this ident in coldata")
  obj@coldata$active.ident <- obj@coldata[[ident]]
  obj
}
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
  x.spread$bc <- paste0(sample,"_",x.spread$bc)
  x.matrix <- Matrix::Matrix(as.matrix(x.spread[, -1]),
                             sparse = TRUE)
  x.matrix[is.na(x.matrix)] <- 0
  rownames(x.matrix) <-  x.spread$bc
  # sort matrix
  x.matrix <- .sortMatrix(x.matrix, "col")
  x.matrix <- .sortMatrix(x.matrix, 'row')
  # initialize coldata and rowdata
  cell <- colnames(x.matrix)
  x.coldata <- data.frame(
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

# =============================== brtInputome ==================================
#' The inputData class
#' @description The inputData class stores the data and describe the column feature
#' of the data
#' @slot data a brtData obj, stores the raw counts,normalized ... matrix
#' @slot coldata usually describe the region, rvcounts,and spk counts of each
#' region. brtInputome obj usually be a part of unity obj, the rowdata has
#' been described in the parent object.
setClass("inputData",
         contains = "brtData",
         slots = c(coldata = "data.frame"))
inputData <- function(counts = NA,coldata = NA,focus = "counts") {
  new("inputData",
      data = list(counts = Matrix::Matrix(counts)),
      focus = focus,
      coldata = coldata)
}
#' The bcInputome class
#' @slot rowdata a data.frame, with conserved colnames bc,rvcounts.raw
#' @slot assays a list contains inputData objects
#' @slot sample a character describe the name of sample
#' @slot active_assay a character describe which assay should be used in the
#' downstream analysis
setClass("bcInputome",
         slots = c(rowdata = "data.frame",
                   assays = "list",
                   sample = "character",
                   active_assay = "character"))
bcInputome <-
  function(tbs,
           sample,
           bc,
           spk = NULL) {
    message("start construct matrix...")
    tb.bc <- tbs[[bc]]
    tb.bc %>%
      arrange(id) %>%
      tidyr::spread(id, counts) -> x
    x$bc <- paste0(sample,"_",x$bc)
    x.matrix <- Matrix::Matrix(as.matrix(x[,-1]),
                               sparse = TRUE)
    x.matrix[is.na(x.matrix)] <- 0
    rownames(x.matrix) <-  x$bc
    # sort matrix by rowsum
    x <- .sortMatrix(x.matrix, 'row')
    message("build brtInputome obj...")
    x.rowdata <- data.frame(bc = rownames(x),
                            rvcounts.raw = Matrix::rowSums(x))
    # raw assay
    x.coldata <- data.frame(id = colnames(x),
                            rvcounts.raw = Matrix::colSums(x))
    if (!is.null(spk)) {
      tb.spk <- tbs[[spk]]
      tb.spk <- rename(tb.spk, spk = counts)
      x.coldata <- left_join(x.coldata, tb.spk, by = "id")
    }
    assay.raw <- inputData(x,x.coldata)
    y <-  new(
      "bcInputome",
      rowdata = x.rowdata,
      assays = list(raw = assay.raw),
      sample = sample,
      active_assay = "raw"
    )
    message("done")
    return(y)
  }
setMethod("show","bcInputome",function(object){
  cat("an inputome obj with ", nrow(object@rowdata)," features","\n\n")
  purrr::walk2(object@assays,names(object@assays),function(x,y){
    cat("assay :",y,"\n",
        nrow(x@coldata)," samples","\n",
        "data focus: ",x@focus,"\n")
  })
  cat("active assay: ", object@active_assay)
})
setClass("brtCluster",
         slots = c(dist = 'list',
                   objs = 'list',
                   reduction = "list"))
#' The scInputome class
#' @description a part of brtunity obj, the rowdata is shared with parent
#' rowdata
setClass("scInputome",
         contains = "bcInputome",
         slots = c(clusters = "brtCluster"))
scInputome <- function(data,rowdata,coldata){
  idx.row <- match(rownames(data),rowdata$name)
  rowdata<- rowdata[idx.row,]
  new("scInputome",
      rowdata = rowdata,
      assays = list(subregion = inputData(data,coldata)))

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
                   bulk = 'scInputome',
                   tools = 'list'))
#' brtUnity
#' @param starter a list of starter objs
#' @param sc a list of brtStarter objs records single cell long-range inputs
#' @param bulk a list of scInputome objs records long-range bulk inputs
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
      arrange(cell) %>%
      filter(cell %in% x@coldata$name)
  })
  names(bcmaps) <- samples
  names(starter) <- samples
  metadata <- purrr::map(starter,~.x@coldata)
  # merge starter
  purrr::map2(starter,bcmaps,function(x,y){
    counts <- x@data@data$counts
    input.idx <- which(x@coldata$state == "input")
    input.counts <- rowSums(counts[,input.idx])
    bc.idx <- which(x@rowdata$unique.bc == TRUE &input.counts > 0)
    counts.sub <- counts[bc.idx,input.idx]
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
  if (! is.null(bulk)) {
    samples <- unlist(purrr::map(bulk,~.x@sample))
    idxs <- purrr::map(bulk,function(x){
      y.coldata <- x@assays$raw@coldata
      idx <- match(brtParams$regions$subregion,y.coldata$subregion)
    })
    data <- purrr::pmap(list(bulk, bcmaps,idxs), function(x,y,z) {
      # merge counts_adj
      if (!"data" %in% names(x@assays$raw@data))
        stop(x@sample,
             ": please estimate and remove noise before run brtUnity")
      data <- x@assays$raw@data$data
      data.merge <- mergeBc(data, y)
      data.merge <- data.merge[,z]
      idx.remove <- which(rowSums(data.merge) == 0)
      message(x@sample," remove ", length(idx.remove),
              " cells with 0 bulk connection counts")
      data.merge <- data.merge[-idx.remove, ]
      colnames(data.merge) <- x@assays$raw@coldata$subregion[z]
      data.merge
    })
    rowdatas <- purrr::map2(starter,data,function(x,y){
      idx <- match(rownames(y),x@coldata$name)
      y.rowdata <- x@coldata[idx,]
      y.rowdata <- select(y.rowdata,name)
    }) %>%
      do.call(rbind,.)
    rownames(rowdatas) <- NULL
    data <- do.call(rbind,data)
    coldata <- data.frame(subregion = colnames(data),
                          rvcounts.adj = colSums(data))
    assays.raw <- new("inputData",
                  coldata = coldata,
                  data = list(data = data),
                  focus = "data")
    bulkobjs <- new(
      "scInputome",
      rowdata = rowdatas,
      assays = list(raw = assays.raw),
      sample = samples,
      active_assay = "raw"
    )
  } else {
    bulkobjs <- new("scInputome")
  }
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
  rowdata <- obj@coldata[rowname.idx[which(!is.na(rowname.idx))],]
  colname.idx <- match(colnames(counts), obj@coldata$name)
  coldata <- obj@coldata[colname.idx[which(!is.na(colname.idx))],]
  coldata$rvcounts.unique <- colSums(counts)
  new("brtExperiment",
      data = data,
      rowdata = rowdata,
      coldata = coldata)
}

