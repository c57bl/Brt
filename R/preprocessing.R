# ============================= noise bc  ======================================

#' brtSetNoiseBc
#' Estimate the noise level of bcs
#' @description The data is noisy, it is important to identify noise. RV seldom
#' infect non-neuronal cells from neurons, thus the bcs expressing levels of
#' non-neuronal groups reflects the noise level of experiment. By measure the
#' median + threshold * IQR, the threshold of signal can be determined.
#' @param obj an brtstarter object, the idents of coldata should contain Neurons
#' and at least one non-neuronal ident.
#' @param background a character vector, describe which non-neuronal cluster used
#' to calculate the threshold. Using all non-neuronal groups if none is provided.
#' @param threshold a numeric value, default to be 3
#' @param signal a character, describe which ident is signal, default is Neurons
#' @export
brtSetNoiseBc <- function(obj,background = NULL,threshold = 3,
                          signal = "Neurons"){
  result <- getNoiseBc(obj,background,threshold,signal)
  obj@rowdata$noise <- FALSE
  obj@rowdata$noise[obj@rowdata$bc %in% result$noise] <- TRUE
  message("set ",length(result$noise),"/",nrow(obj@rowdata), " noise")
  return(obj)
}

# a sub function of brtSetNoiseBc
# a sub function of brtPlotNoise
# return a list
getNoiseBc <- function(obj,background = NULL,threshold = 3,signal) {
  # n, neuron; nn, non-neuron
  if (is.null(background))
    NN.idx <- which(obj@coldata$ident != signal)
  else
    NN.idx <- which(obj@coldata$ident %in% background)
  NN.data <- obj@data@data$counts[, NN.idx]
  NN.data %>%
    Matrix::rowSums() -> bc.nn
  N.idx <- which(obj@coldata$ident == signal)
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
  uvlist <- uv$bc
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
  obj@coldata$rvcounts.unique <- 0
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

# ========================== inputome =========================================
