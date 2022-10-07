# ============================= noise bc  ======================================
brtSetNoiseBc.brtStarter <- function(obj,background = NULL,threshold = 3,
                          signal = "Neurons"){
  result <- getNoiseBc(obj,background,threshold,signal)
  obj@rowdata$noise <- FALSE
  obj@rowdata$noise[obj@rowdata$bc %in% result$noise] <- TRUE
  message("set ",length(result$noise),"/",nrow(obj@rowdata), " noise")
  return(obj)
}
setMethod("brtSetNoiseBc",signature = "brtStarter",brtSetNoiseBc.brtStarter)

# a sub function of brtSetNoiseBc
# a sub function of brtPlotNoise
# return a list
getNoiseBc<- function(obj,background = NULL,threshold = 3,signal) {
  # n, neuron; nn, non-neuron
  if (is.null(background))
    NN.idx <- which(obj@coldata$active.ident != signal)
  else
    NN.idx <- which(obj@coldata$active.ident %in% background)
  NN.data <- obj@data@data$counts[, NN.idx]
  NN.data %>%
    Matrix::rowSums() -> bc.nn
  bc.nn <- bc.nn[which(bc.nn > 0)]
  N.idx <- which(obj@coldata$active.ident == signal)
  N.data <- obj@data@data$counts[, N.idx]
  N.data %>%
    Matrix::rowSums() -> bc.n
  cutoff <- median(bc.nn) + threshold * IQR(bc.nn)
  bc.n <- data.frame(bc = names(bc.n), signal = bc.n)
  bc.nn <- data.frame(bc = names(bc.nn), Background = bc.nn)
  # for brtPlotNoise
  bcs <- left_join(bc.n,bc.nn,by = "bc") %>%
    tidyr::gather("source","counts",-bc)
  bcs$source[which(bcs$source == "signal")] <- signal
  bcs$counts[is.na(bcs$counts)] <- 0
  # noise is below cutoff bcs in signal group
  noise <- unique(subset(bcs,counts <= cutoff & source == signal)$bc)
  list(data = bcs,
       cutoff = cutoff,
       noise = noise)
}


# ========================== unique bc/virus ===================================

#' brtSetUniqueVirus
#' set unique bc particles based on umi counts
#' @description Assuming that the probability of bcs being selected from virus
#' library is linearly correlated with the umi counts of bcs, the total virus
#' particles number and the unique probability of each bc could be estimated.
#' @param obj a brtstarter obj
#' @param S The total number of particles, default to be NULL and S will be estimated
#' by the observed bcs.
#' @param P The cutoff unique probability of each bc.
#' @export
brtSetUniqueVirus <- function(obj,S = NULL,P = 0.95){
  if (! "noise" %in% colnames(obj@rowdata))
    stop("can not find noise column in rowdata, please run brtSetNoiseBc")
  if (is.null(S)){
    N <- length(which(obj@rowdata$noise == FALSE))
    uv <- brtUniqueVirus(obj@virus,N = N, P.min = P)
  } else if (S > 0 & is.integer(S)) {
    uv <- brtUniqueVirus(obj@virus,S = S, P.min = P)
  } else stop("S should be a positive integer")
  sample <- obj@coldata$sample[1]
  # bc name in each sample is sample_bc
  uvlist <- paste0(sample,"_",uv$bc)
  obj@rowdata$unique.virus <- FALSE
  obj@rowdata$unique.virus[which(obj@rowdata$bc %in% uvlist &
                                   obj@rowdata$noise == FALSE)] <- TRUE
  message("set ",sum(obj@rowdata$unique.virus),"/",N," unique virus")
  obj
}

#' brtUniqueVirus
#' @description describe unique labeling rate of the brtVirus object
#' @param x an brtvirus object
#' @param S an integer describe the number of infected cells
#' @param N an integer describe the types of barcodes, ignored if S is provided
#' @export
brtUniqueVirus <- function(x,S = NULL,N = NULL,P.min = 0.95){
  if (is.null(S) & is.null(N)) stop("S or N should be provided")
  else if (is.null(S)) {
    S.idx <- min(which(x@bc_type$PN > N))
    S <- x@bc_type$S[S.idx]
    message("estimate S: ", S)
  }
  x <- x@bc
  bc <- x$bc
  counts <- x$counts
  N <- sum(counts)
  P <- purrr::map(counts,~((N-.x)/N)^S+((N-.x)/N)^(S-1)*.x*S/N) %>%
    unlist
  idx <- which(P>=P.min)
  results <- data.frame(bc=bc[idx],vcount=counts[idx],P=P[idx])
  return(results)
}

#' brtSetUniqueBc
#' @description unique bc should be non-noise bc, unique virus bc, and is distinct
#' among starters. Thus, brtSetNoiseBc, brtSetUniqueVirus and brtSetCellState should
#' have been run.
#' @param obj a brtstarter obj
#' @param N cutoff of counts
#' @param U cutoff of proportion within cells
#' @param P cutoff of proportion across cells
#' @return update rowdate and coldata with useful columns
#' cell.p: for bc i, cell.p is max(bci)/sum(bci) and
#' cell: cell with max bci
#' cell.u: for bc i and cell j, cell.u is bci/sum(cell.u)
#' unique.bc, logical, which bc is unique and useful bc.
#' rvcounts.unique: total unique counts within one cell, useful for downstream
#' analysis
#' @importFrom Matrix rowSums colSums t
brtSetUniqueBc <- function(obj,N = 3,U = 0.3,P = 0.8){
  if (P <= 0.5)
    stop("P should > 0.5")
  data <- obj@data@data$counts
  data[data < N] <- 0
  if (! "state" %in% colnames(obj@coldata))
    stop("please run brtSetCellState")
  if (! "unique.virus" %in% colnames(obj@rowdata))
    stop("please run brtSetUniqueVirus")
  idx_cell <- which(obj@coldata$state == "starter")
  idx_bc <- which(obj@rowdata$unique.virus)
  data.sub <- data[idx_bc,idx_cell]
  # +1 here to avoid /0, but doesn't affect the result
  data.p <- data.sub / (rowSums(data.sub) + 1)
  df.r <- apply(data.p,1,function(x){
    x.max <- max(x)
    max.idx <- which(x == x.max)
    # zero counts or equal counts
    if (length(max.idx) > 1) {
      cell <- NA_character_
      cell.p <- NA_real_
    } else {
      cell <- names(x)[max.idx]
      cell.p <- x[max.idx]
    }
    data.frame(cell = cell,
               cell.p= cell.p)
  }) %>%
    do.call(rbind,.)
  df.r$bc <- rownames(data.p)
  # subset data with row-wise proportion > P
  data.p[data.p <= P] <- 0
  data.p[data.p > P] <- 1
  data.reduce <- data.sub * data.p
  # calculate cell.u should use all bc counts
  data.u <- t(t(data.reduce) / (colSums(data[,idx_cell]) + 1))
  cell.u <- apply(data.u, 1, max)
  data.u[data.u < U] <- 0
  cell.u[cell.u == 0] <- NA_real_
  df.r$cell.u <- cell.u
  # update results to rowdata
  obj@rowdata <- left_update(obj@rowdata,df.r,by = "bc")
  obj@rowdata$unique.bc <- FALSE
  obj@rowdata$unique.bc[which(obj@rowdata$cell.u > U)] <- TRUE
  # calculate unique counts of each cell
  cell.uc <- apply(data.u, 2, sum) * (colSums(data[,idx_cell]) + 1)
  obj@coldata$rvcounts.unique <- NA_real_
  idx <- match(names(cell.uc), obj@coldata$name)
  obj@coldata$rvcounts.unique[idx] <- cell.uc
  message("set ", length(which(cell.u > U))," unique bcs across ",
          length(which(cell.uc > 0))," cells")
  return(obj)
}



# =========================== cell identity ====================================
#' brtGetInfection
#' @description  The negative cutoff can not be inferred from the background noise
#' distribution, instead, the 0 counts of tva or rv was treated as ground truth for
#' the negative. Under this estimate, the infection rate of tva- cells was
#' n_rv+/n_tva-, which is the spread efficiency,se. se is quite important in
#' inferring the  initial infection efficiency of TVA+ cells, and the cutoff of
#' tva- cells.
#' @param obj an brtstarter obj
#' @param tva.c cutoff of tva
#' @param rv.c cutoff of rvcounts.raw
#' @return A list contains:
#'
#' RPT: rv proportion in tva positive cells;
#'
#' SE: spread efficiency;
#'
#' IE: infection efficiency,
#'
#' SP: estimated spread proportion at tva+ cells,
#'
#' record$IP: infection proportion at different tva counts
#'
#' record$RP: rv+ proportion at different tva counts
#' @export
brtGetInfection <- function(obj,tva.c,rv.c) {
  tva.cell <- subset(obj@coldata, tva > tva.c)$name
  tva.n.cell <- subset(obj@coldata, tva == 0)$name
  rv.cell <- subset(obj@coldata, rvcounts.raw > rv.c)$name
  # ie infection efficiency, se spread efficiency
  RPT <- length(intersect(tva.cell,rv.cell))/length(tva.cell)
  se <- length(intersect(tva.n.cell,rv.cell))/length(tva.n.cell)
  ie <- (RPT - se) / (1 - se)
  # R.x RV+ proportion at x tva counts
  purrr::map(0:tva.c, function(x) {
    tva.c.c <- subset(obj@coldata, tva == x)$name
    R.x <- length(intersect(tva.c.c, rv.cell)) / length(tva.c.c)
  }) %>% unlist -> R.x
  # P.x proportion of cell number at x tva counts
  purrr::map(0:tva.c, function(x) {
    tva.c.c <- subset(obj@coldata, tva == x)$name
    p.x <-  length(tva.c.c)/nrow(obj@coldata)
  }) %>% unlist -> P.x
  # PI.x infection proportion at x tva counts
  purrr:::map(R.x, function(x) {
    ((x - se) / (1 - se))
  }) %>%
    unlist() -> PI.x
  # PS proportion of spread at tva+ cells
  SP <- (1-ie) * se
  record <- data.frame(tva.counts = 0:tva.c,IP = PI.x,RP = R.x)
  return(list(
    RPT = RPT,
    SE = se,
    IE = ie,
    SP = SP,
    record = record
  ))
}

#' brtSetCellIdentity
#' @description set state of cells. states including naive: bc negative cells,
#' starter: bc+ tva+ cells and input: bc+ tva- cells. The logical expression
#' starter and input is free to set, butshould always include tva and rvcounts.raw.
#' @param starter 	logical expression indicating the starter rows
#' @param input logical expression indicating the input rows
#' @export
brtSetCellState <- function(obj,starter,input){
  starter <- substitute(starter)
  starter.r <- eval(starter, obj@coldata, parent.frame())
  if (!is.logical(starter.r))
    stop("'starter' must be logical")
  input <- substitute(input)
  input.r <- eval(input, obj@coldata, parent.frame())
  if (!is.logical(input.r))
    stop("'input' must be logical")
  obj@coldata$state <- "naive"
  obj@coldata$state[starter.r] <- "starter"
  obj@coldata$state[input.r] <- "input"
  message("set ", sum(starter.r)," starters and ",sum(input.r)," inputs")
  obj
}

# =============================   bcinputome =================================
#' brtAdjustSpk
#' @description adjust counts using spike-in counts
#' @param obj an bcinputome object
#' @param assay an character describe which assay will be used
#' @return an bcinputome object, with counts_adj added in specific assay
brtAdjustSpk <- function(obj,assay = "raw") {
  checkBcInputParams(obj,assay = assay)
  data <- obj@assays[[assay]]@data$counts
  coldata <- obj@assays[[assay]]@coldata
  spk <- coldata$spk
  rvcounts <- coldata$rvcounts.raw
  rvcounts.adj <- round(rvcounts / spk * max(spk))
  # # signal.f: sampling frequency
  # data.f <- as.matrix(t(t(data) / rvcounts))
  # data.n <- matrix(rep(rvcounts.adj,nrow(data)),
  #                    nrow = nrow(data),
  #                    byrow = T)
  # # use 0.05 pdf
  # data.adj <- purrr::map2(data.f,data.n,function(x,y){
  #   if (x > 0) {
  #     qbinom(0.05,y,x)
  #   } else {
  #     0
  #   }
  # })
  # data.adj <- Matrix(unlist(data.adj),nrow = nrow(data))
  data.adj <- round(t(t(data) / spk * max(spk)))
  dimnames(data.adj) <- dimnames(data)
  # add obj
  coldata$rvcounts.adj <- colSums(data.adj)
  obj@assays[[assay]]@coldata <- coldata
  obj@assays$raw@data$counts_adj <- data.adj
  obj
}

#' brtEstimateBcInputNoise
#' @description estimate the noise level of inputs via the expression in the
#' negative control regions in each slice.
#' @param obj an bcInputome obj
#' @param threshold an integer control the threshold of noise, threshold is
#' estimated as median(nc) + threshold * IQR(nc), default to be 3
#' @param min.counts an integer, bcs with total counts less than min.counts
#' were excluded
#' @details #  steps 1. remove counts less than min.counts; 2. perform proportion
#' normalization; 3 calculate nc distribution; 4. determine nc threshold,
#' all bcs are equal; 5. replace threshold with exceed threshold; 6. for bcs with
#' sum counts less than min.counts, replace the threshold with the observed
#' proportion in counts; 7. adjust threshold with spk counts
#' @return an bcinputome obj, with added slot:noise in the data in selected assay
brtEstimateBcInputNoise <- function(obj,threshold = 3,min.counts = 10) {
  checkBcInputParams(obj)
  if (! "counts_adj" %in% names(obj@assays$raw@data))
    stop("please run brtAdjustSpk before set noise")

  coldata <- obj@assays$raw@coldata
  if (! "Slice" %in% colnames(coldata) | ! "NC" %in% colnames(coldata))
    stop("coldata should contain slice and NC")
  counts <- obj@assays$raw@data$counts
  # remove counts less than min.counts, and calculate proportion
  counts.rowsum <- rowSums(counts)
  mat.rowsum <- Matrix(rep(counts.rowsum,ncol(counts)),
                       nrow = nrow(counts))
  idx.low <- which(counts.rowsum < min.counts)
  proportions <- counts / counts.rowsum
  # calculate nc distribution
  nc.idx <- which(coldata$NC == TRUE)
  signal.idx <- which(coldata$NC == FALSE)
  nc.proportions <- proportions[,nc.idx]
  nc.proportions[idx.low,] <- 0
  nc.threshold <- apply(as.matrix(nc.proportions),2,function(x){
    x <- x[which(x > 0)]
    if (length(x) == 0) 0
    else
      median(x) + threshold * IQR(x)
  })
  # nc.t.m nc threshold matrix
  nc.t.m <- Matrix(rep(nc.threshold,nrow(counts)),
                 nrow = nrow(counts),
                 byrow = T)
  idx.replace <- Matrix::which(nc.t.m < nc.proportions)
  # replace the nc.threshold with exceeded proportion
  nc.t.m[idx.replace] <- nc.proportions[idx.replace]
  # generate nm: noise matrix
  coldata$idx <- 1:nrow(coldata)
  nm.idxs <- lapply(unique(coldata$Slice), function(x) {
    coldata[coldata$Slice == x, ]$idx
  })
  nm <- Matrix(0,nrow = nrow(counts), ncol = ncol(counts))
  for (i in 1:length(nm.idxs)) {
    nm[,nm.idxs[[i]]] <- nc.t.m[,i]
  }
  nm[idx.low,] <- proportions[idx.low,]
  # trans to counts
  nm <- nm * mat.rowsum
  # adjust spk
  nm.adj <- t(t(nm) * (coldata$rvcounts.adj / coldata$rvcounts.raw))
  counts.adj <- obj@assays$raw@data$counts_adj
  for (i in nc.idx) {
    idx <- which(counts.adj[,i] > nm.adj[,i])
    nm.adj[idx,i] <- counts.adj[idx,i]
  }
  dimnames(nm.adj) <- dimnames(counts)
  nm.adj <- Matrix(nm.adj)
  obj@assays$raw@data$noise <- nm.adj
  obj
}

# plotnoise
# removenoise
#' brtRemoveBcInputNoise
#' @description remove noise of inputs
#' @param obj an bcInputome obj
#' @param assay an character describe which assay will be used
#' @return an scInputome obj, with with added slot:data in the data in selected
#'  assay, the bcs with 0 total counts were removed
brtRemoveBcInputNoise <- function(obj,assay = "raw") {
  checkBcInputParams(obj,assay = assay,focus = "noise")
  noise <- obj@assays[[assay]]@data$noise
  counts.adj <- obj@assays[[assay]]@data$counts_adj
  data <- round(counts.adj - noise)
  # bc below threshold and less than 1 will be treated as 0
  data[Matrix::which(data < 0)] <- 0
  idx <- which(rowSums(data) > 0)
  # remove bc counts with 0 bc counts
  data <- data[idx,]
  obj@assays[[assay]]@data$data <- data
  obj@assays[[assay]]@focus <- "data"
  obj
}

# normalize input----
brtNormalizeData <- function(x,method){
  if (method == "proportion") {
    x <- x / rowSums(x)
  } else if (method == "scale") {
    x <- t(as.matrix(x))
    x <- t(scale(x))
  } else {
    stop("invalid method")
  }
  x
}

brtNormalizeInput.bcInputome <- function(obj,method,assay,focus) {
  checkBcInputParams(obj,assay,focus)
  data <-  obj@assays[[assay]]@data[[focus]]
  data <- brtNormalizeData(x = data,method = method)
  obj@assays[[assay]]@data[[paste0(focus,"_",method)]] <- data
  obj
}

brtNormalizeInput.brtUnity <- function(obj,method,assay,focus) {
  checkBcInputParams(obj@bulk,assay,focus)
  data <-  obj@bulk@assays[[assay]]@data[[focus]]
  data <- brtNormalizeData(x = data,method = method)
  obj@bulk@assays[[assay]]@data[[paste0(focus,"_",method)]] <- data
  obj
}
setMethod("brtNormalizeInput","bcInputome",brtNormalizeInput.bcInputome)
setMethod("brtNormalizeInput","brtUnity",brtNormalizeInput.brtUnity)
# subset

checkBcInputParams <- function(obj,assay = NULL,focus = NULL) {
  if (! class(obj) %in% c("bcInputome","scInputome"))
    stop("obj should be a Inputome object")
  else if (!is.null(assay)) {
    if (!assay %in% names(obj@assays))
      stop("can not find assay: ", assay)
  }
  else if (!is.null(focus)) {
    if (!focus %in% names(obj@assays[[assay]]@data))
      stop("can not find data: ", focus, " in assay: ", assay)
  }
}

# =========================== merge and bin regions ============================
brtMergeRegion.inputData <- function(obj,target,focus = "data") {
  data <- obj@data[[focus]]
  data.colname <- colnames(data)
  regions <- brtParams$regions
  bin_regions <- unique(brtParams$regions[[target]])
  idxs <- lapply(bin_regions,function(x){
    subregion.current <- regions[which(regions[[target]] == x),]$subregion
    match(subregion.current,data.colname)
  })
  data.merge <- lapply(idxs,function(x){
    if (length(x) > 1)
      rowSums(data[,x])
    else
      data[,x]
  })
  data.merge <- Matrix(do.call(cbind,data.merge))
  colnames(data.merge) <- bin_regions
  coldata <- data.frame(target = bin_regions,
                        rvcounts.adj = colSums(data.merge))
  colnames(coldata)[1] <- target
  new("inputData",
      coldata = coldata,
      data = list(data = data.merge),
      focus = "data")
}
brtMergeRegion.unity <- function(obj,target,focus) {
  newassay <- list(brtMergeRegion.inputData(obj@bulk@assays$raw,
                                            target = target,
                                            focus = focus))
  names(newassay) <- target
  obj@bulk@assays <- c(obj@bulk@assays,newassay)
  obj
}
setMethod("brtMergeRegion","brtUnity",brtMergeRegion.unity)
