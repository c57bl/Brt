# ============================= description ====================================
#' brtEstimateBc
#' @param obj an brtstarter object with data from neurons and nor-neurons
brtEstimateBc <- function(obj,background = NULL,threshold = 3){
  if (is.null(background)) NN.idx <- which(obj@coldata$ident != "Neurons")
  else NN.idx <- which(obj@coldata$ident %in% background)
  NN.data <- obj@data@counts[, NN.idx]
  NN.data %>% Matrix::rowSums() -> bc.nn
  N.idx <- which(obj@coldata$ident == "Neurons")
  N.data <- obj@data@counts[, N.idx]
  N.data %>% Matrix::rowSums() -> bc.n
  bc.n <- data.frame(bc = names(bc.n), counts_n = bc.n)
  threshold <- median(bc.nn) + threshold * IQR(bc.nn)
  bc.pass <- subset(bc.n, counts_n > threshold)
  return(bc.pass)
}
# ============================= get cutoff =====================================
#' brtGetCutoff
#' In fact, the positive cutoff could be inferred form background, using threshold
#' such as 3*IQR, but the negative cutoff can never be inferred from background
#' levels using datasets with large variances
#' @param obj an brtexperiment obj
#' @param feature an colname of coldata
#' @param signal an celltype from ident column
#' @param background celltype/celltypes from ident column
#' @param threshold an integer
#' @return an list contain the cutoff and the plot of data
#' @export
brtGetCutoff <- function(obj,
                         feature = NULL,
                         signal = "Neurons",
                         background = NULL,
                         threshold = 3){
  if (is.null(background)) bg.idx <- which(obj@coldata$ident != signal)
  else bg.idx <- which(obj@coldata$ident %in% background)
  s.idx <- which(obj@coldata$ident == signal)
  data <- obj@coldata[[feature]]
  data.signal  <- data[s.idx]
  data.bg <- data[bg.idx]
  data.bg <- data.bg[which(data.bg > 0)]
  cut.line <- median(data.bg) + threshold * IQR(data.bg)
  rbind(data.frame(source = "signal",counts = data.signal),
        data.frame(source = "background", counts = data.bg)) %>%
    ggplot(aes(source,counts)) +
    geom_jitter(color = "gray",size = 1) +
    geom_violin() +
    scale_y_log10() +
    geom_hline(yintercept = cut.line + 0.5,color = "red") +
    theme_classic() +
    labs(title = feature) -> p.line
  return(list(cutoff = cut.line,
              plot = p.line))
  # step based functions
  # if (class(obj)!='brtStarter') stop('invalid obj')
  # groups <- obj@coldata[[split_by]]
  # if (is.null(background)) background <- setdiff(groups,signal)
  # if (source=='coldata'){
  #   if(is.na(signal) || is.na(feature)) stop('incomplete input')
  #   else data <- select(obj@coldata,split_by,cell,feature) %>%
  #       rename(counts = feature)
  # } else {
  #     as_tibble(t(as.matrix(obj@data@counts))) %>%
  #     mutate(cell = colnames(obj@data@counts)) %>%
  #     tidyr::gather("bc","counts",-cell) %>%
  #     subset(counts >0) %>%
  #     left_join(select(obj@coldata,split_by,cell),by = "cell") %>%
  #     select(ident,cell,counts) -> data
  # }
  # data.signal <- subset(data,ident %in% signal)
  # data.back <- subset(data,ident %in% background)
  # purrr::map(steps,~.brtgc.getp(data.signal,data.back,.x)) %>%
  #   do.call(rbind,.) -> data
  # data$f <- data$f/max(data$f)
  # data %>%
  #   ggplot(aes(step,p)) +
  #   geom_line(color = 'red') +
  #   geom_line(aes(step,f),color = 'blue') +
  #   theme_light() +
  #   geom_hline(yintercept = 0.95) +
  #   geom_hline(yintercept = data$rand[1]) +
  #   ggtitle(feature) -> plot
  # return(list(data = data,plot = plot))

}

# x signal df
# y  background df
.brtgc.getp <- function(x,y,step){
  node.signal <- length(unique(x$cell))
  node.back <- length(unique(y$cell))
  x <- subset(x,counts >= step)
  y <- subset(y,counts >= step)
  counts.signal <- nrow(x)
  counts.back <- nrow(y)
  counts.all <- counts.signal + counts.back
  node.all <-  node.signal + node.back
  p <- (counts.back*node.all)/(node.back*counts.all)
  n <- counts.signal
  rand <- node.back/(node.signal+node.back)
  data.frame(step = step,p = 1-p,f = n,rand = rand)
}

#' brtGetNoise
#' @description  Estimate noise level of inputs using NC expression level
#' @param obj an brtInputome object
#' @param threshold an integer control the threshold of noise. Signal should
#' above median(noise)+threshold*mad(noise)
#' @param plot logic
#' @export
brtGetNoise<-function(obj,method="IQR",gain=3,usespk = TRUE,
                      clean_nc = TRUE){
  if (class(obj) != "brtInputome") stop("obj should be a brtInputome object")
  if (usespk) {
    curr.scale<-max(obj@coldata$spk)
    curr.data<-obj@data@counts
    curr.data<-Matrix::t(curr.data)/obj@coldata$spk*curr.scale
    curr.data<-Matrix::t(curr.data)
  } else curr.data <- obj@data@counts
  slice.list <- unique(obj@coldata$Slice)
  obj.coldata<-obj@coldata
  obj.coldata$idx<-1:nrow(obj.coldata)
  noise.dfs<-list()
  for (slice in slice.list) {
    curr.coldata<-subset(obj.coldata,Slice==slice)
    nc.idx<-subset(curr.coldata,NC==TRUE)$idx
    if (length(nc.idx)==0) stop(paste0("can not find nc in slice ",
                                       as.character(slice)))
    sample.idx<-subset(curr.coldata,NC==FALSE)$idx
    if (length(sample.idx)==0) warning(paste0("can not find sample in slice ",
                                              as.character(slice)))
    nc.data <- curr.data[,nc.idx]
    sample.data <- curr.data[,sample.idx]
    sample.data.merge <- Matrix::rowSums(sample.data)
    data.df <- data.frame(noise=nc.data,sample=sample.data.merge,
                          bc=obj@rowdata$bc,slice=slice,
                          idx=1:length(nc.data))
    if (clean_nc) {
      data.df.noise <- subset(data.df,sample==0 & noise>0)
    } else data.df.noise <- subset(data.df,noise>0)
    noise.data <- data.df.noise$noise
    if (method=="IQR") {
      noise.quan <- quantile (noise.data)
      noise.IQR <- IQR(noise.data)
      noise.threshold <- noise.quan[2]+gain*noise.IQR
    }
    # record outer and parameters
    data.df$outer <- FALSE
    data.df$outer[which(data.df$noise>noise.threshold)] <- TRUE
    data.df$threshold <- noise.threshold

    # record current df
    noise.dfs[[slice]] <- data.df
  }
  # merge
  noise.dfs <- do.call(rbind,noise.dfs)
  obj@tools$"brtGetNoise" <- list(method=method,gain=gain,results=noise.dfs)
  # plot noise and threshold
  return(obj)
}
# ============================= cut off ========================================
#' brtCutoff
#' @description brtCutoff will generate data from counts after remove the
#' background
#' @export
setGeneric("brtCutoff",
           function(obj,cutoff=0,...) standardGeneric("brtCutoff"),
           signature = "obj")
setMethod("brtCutoff","brtStarter",function(obj,cutoff){
  if (!is.numeric(cutoff)) stop('cutoff should be a number')
  currdata<-obj@data@counts
  currdata[currdata<=cutoff]<-0
  obj@data@data<-currdata
  obj@coldata$rvcounts<- Matrix::colSums(currdata)
  obj@rowdata$rvcounts<- Matrix::rowSums(currdata)
  return(obj)
})
#' brtCutoff-brtInputome
#' @export
setMethod("brtCutoff","brtInputome",function(obj,usespk = TRUE){
  if (is.null(obj@tools$brtGetNoise)) stop("please run brtGetCutoff")
  if (usespk) {
    curr.scale<-max(obj@coldata$spk)
    curr.matrix<-obj@data@counts
    curr.matrix<-Matrix::t(curr.matrix)/obj@coldata$spk*curr.scale
    curr.matrix<-Matrix::t(curr.matrix)
  } else curr.matrix<-obj@data@counts
  obj.coldata <- obj@coldata
  obj.coldata$idx <- 1:nrow(obj.coldata)
  noise.dfs<-obj@tools$brtGetNoise$results
  curr.matrix.list <- purrr::map(unique(noise.dfs$slice),function(x){
    curr.dfs <- subset(noise.dfs,slice==x)
    curr.dfs.outer <- subset(curr.dfs,outer==T)
    curr.coldata <- subset(obj.coldata,Slice==x)
    sample.data <- curr.matrix[,curr.coldata$idx]
    noise.threshold <- curr.dfs$threshold[1]
    sample.data[sample.data<noise.threshold] <- 0
    noise.outer.mask <- Matrix::Matrix(0,nrow=nrow(curr.dfs),
                                       ncol=nrow(curr.coldata))
    rownames(noise.outer.mask) <- curr.dfs$bc
    noise.outer.mask[match(curr.dfs.outer$bc,
                           rownames(noise.outer.mask)),] <- curr.dfs.outer$noise
    sample.data <- sample.data-noise.outer.mask
    sample.data[sample.data<0] <- 0
    return(sample.data)
  })
  curr.matrix <- do.call(cbind,curr.matrix.list)
  obj@data@data <- curr.matrix
  message("noise removed")
  return(obj)
})
# ============================ set ident =======================================
#' brtSetStarter
#' @param starter an expression character
#' @param input an expression character
#' @export
brtSetStarter<-function(obj,starter=NULL,input=NULL){
  coldata<-obj@coldata
  if (! is.null(starter)){
    starterlist<-coldata[eval(parse(text = starter),coldata),]
    coldata$starter<-FALSE
    coldata$starter[coldata$cell %in% starterlist$cell]<-TRUE
  } else {
    starterlist<-data.frame()
  }
  if (! is.null(input)){
    inputlist<-coldata[eval(parse(text = input),coldata),]
    coldata$inputer<-FALSE
    coldata$inputer[coldata$cell %in% inputlist$cell]<-TRUE
  } else {
    inputlist<-data.frame()
  }
  obj@coldata<-coldata
  cat(nrow(starterlist),' starters, ',nrow(inputlist),' inputs were set \n')
  return(obj)
}

# ========================== set unique ========================================
#' brtUnique
#' @param subset an expression describe the subset of coldata
#' @param rr character vector provide the unique barcode reference.
#' Unique will be False if the barcode is not one of the reference
#' @param cutoff.p inter group cutoff
#' @param cutoff.within within group cutoff
#' @param cutoff.n counts cutoff
#' @export
setGeneric("brtUnique",
           function(obj,
                    subset,
                    cutoff.p = 0.9,
                    cutoff.n = 3,
                    rr = NULL,
                    cutoff.within = 0.3,
                    ...)
             standardGeneric("brtUnique"),
           signature = "obj")

setMethod("brtUnique","brtStarter",function(obj,subset,
                                            cutoff.p,cutoff.n,rr,
                                            cutoff.within){
  message("select unique by barcode proportion")
  objsub <- brtSubset(obj, 'col', subset)
  objsub <- brtSubset(objsub, 'row', 'rvcounts>0')
  data <- objsub@data@data %>%
    as.matrix()
  data.rs <- Matrix::rowSums(data)
  data.cc <- objsub@coldata$rvcounts
  pass.1 <- data/data.rs > cutoff.p
  pass.2 <- apply(data, 1, max, na.rm = T) > cutoff.n
  pass.3 <- t(t(data)/data.cc) > cutoff.within
  maxidx <- purrr::map(1:nrow(data),
                       ~.x+ (which.max(data[.x,])-1) * nrow(data)) %>%
    unlist()
  maxCounts <- data[maxidx]
  maxCell <- colnames(data)[apply(data,1,which.max)]
  maxProportion <- (data/data.rs)[maxidx]
  Proportion_incell <- t(t(data)/data.cc)[maxidx]
  passed <- Matrix::rowSums(pass.1 & pass.2 & pass.3) %>%
    as.logical()
  isunique <- data.frame(
    bc = rownames(data),
    cell = maxCell,
    unique = passed,
    maxPorportion = maxProportion,
    maxCounts = maxCounts,
    Proportion_incell = Proportion_incell
  )
  if (! is.null(rr)) {
    isunique$unique[! isunique$bc %in% rr] <- FALSE
  }
  obj@rowdata <- left_update(obj@rowdata,isunique,by = 'bc',treat.na = FALSE)
  message("calculate unique barcodes proportion within each cell...")
  obj <- .detectUniqueProportion(obj)
  isunique<-subset(isunique,unique==TRUE)
  cat(nrow(isunique),'unique barcodes across', length(unique(isunique$cell)),'cells')
  return(obj)
})
.detectUniqueProportion <- function(obj){
  obj.tibble <- as_tibble(obj@data@data)
  obj.tibble$bc <- rownames(obj@data@data)
  calpro <- function(xx,yy){
    sum(xx[which(yy == T)],na.rm = T)/
      sum(xx,na.rm = T)
  }
  obj.tibble %>%
    tidyr::gather("cell","counts",-bc) %>%
    subset(counts >0) %>%
    left_join(select(obj@rowdata,bc,cell,unique),
              c("cell","bc")) %>%
    group_by(cell) %>%
    summarise(unique.proportion = calpro(counts,unique)) -> obj.tibble
  obj@coldata <- left_update(obj@coldata,obj.tibble,by = "cell",treat.na = 0)
  return(obj)
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

# ============================ normalize =======================================
#' brtNormalize
#' @param method can be one of proportion,zscore,logNormalize
#' @export
setGeneric("brtNormalize",
           function(obj,method,...) standardGeneric("brtNormalize"),
           signature = "obj")
setMethod("brtNormalize","brtStarter",function(obj,method){
  if (method=="proportion"){
    currdata<-obj@data@data
    currdata<-currdata/Matrix::rowSums(currdata,na.rm = T)
    obj@data@normalized<-currdata
    return(obj)
  }
})
setMethod("brtNormalize","brtUnity",function(obj,method,slot,scale.size = 1e4){
  if (slot == "bulk") currdata <- obj@bulk@data@data
  if (method=="proportion"){
    results<-currdata/Matrix::rowSums(currdata,na.rm = T)
    results[is.na(results)] <- 0
  }
  else if (method == "zscore") {
    data.mean <- Matrix::rowMeans(currdata)
    data.sd <- apply(as.matrix(currdata),1,sd)
    results <- purrr::map(1:nrow(currdata),
                         ~(currdata[.x,]-data.mean[.x])/data.sd[.x]) %>%
      do.call(rbind,.) %>%
      Matrix::Matrix()
    results[is.na(results)] <- 0
    colnames(results) <- colnames(currdata)
    rownames(results) <- rownames(currdata)
  }
  else if (method == "logNormalize"){
    data.sum <- Matrix::rowSums(currdata)
    f <- data.sum/scale.size
    results <- currdata/f
    results[results == 0] <- NaN
    results <- log2(results)
    results[is.na(results)] <- 0
  }
  else if (method == "log_min_max"){
    data.sum <- Matrix::rowSums(currdata)
    f <- data.sum/scale.size
    results <- currdata/f
    results[results == 0] <- NaN
    results <- log2(results)
    results[is.na(results)] <- 0
    data.max <- apply(as.matrix(results),1,max)
    data.min <- apply(as.matrix(results),1,min)
    results <- (results - data.min)/(data.max - data.min)
  }
  if (slot == "bulk") obj@bulk@data@normalized <- results
  return(obj)
})
#' brtBinInput
#' @description cutoff inputome matrix by given threshold and generate an bin matrix
#' the number above n and p * max of each row was treated as TRUE.
#' @param n a number describe the numeric threshold
#' @param p a number between 0-1 describe the proportion threshold
#' @export
brtBinInput <- function(obj, n = 0, p = 0) {
  data <- obj@bulk@data@data
  data <- .brtBinInput_Matrix(data, n = n, p = p)
  obj@tools$brtBinInput <- list(data = data,
                                params = list(n = n,
                                              p = p))
  return(obj)
}
brtBinInput_brtExp <- function(obj, n = 0, p = 0) {
  data <- obj@data@data
  data <- .brtBinInput_Matrix(data, n = n, p = p)
  obj@tools$brtBinInput_brtExp <- list(data = data,
                                       params = list(n = n,
                                                     p = p))
  return(obj)
}
.brtBinInput_Matrix <- function(data,n,p){
  data <- data - n
  data[data < 0] <- 0
  idx <- which(Matrix::rowSums(data) > 0)
  idx.backup <- which(Matrix::rowSums(data) == 0)
  data[idx,] <- data[idx,] / apply(as.matrix(data[idx,]), 1 ,max)
  data[data < p] <- 0
  data[data >= p] <- 1
  data[idx.backup,] <- 0
  return(data)
}

