#' @include generics.R
NULL
#============================= update ==========================================
#' left_update
left_update <- function(x, y, by, treat.na = NA) {
  y <- unique(y)
  col.inter <- intersect(colnames(x), colnames(y))
  col.inter <- dplyr::setdiff(col.inter, by)
  x <- x[setdiff(colnames(x), col.inter)]
  x <- left_join(x, y, by = by)
  colnames.to.test <- setdiff(colnames(y), by)
  if (!is.null(colnames.to.test)) {
    for (i in 1:length(colnames.to.test)) {
      idx <- which(is.na(x[, colnames.to.test[i]]))
      x[idx,colnames.to.test[i]] <- treat.na
    }
  }
  return(x)
}
.treatna <- function(x, y) {
  for (i in 1:ncol(x)) {
    x[, i][which(is.na(x[, i]))] <- y
  }
  return(x)
}

# ================================ subset =====================================
setMethod("subsetRow",
          signature(obj = "brtExperiment"),
function(obj, subset, ...) {
  x <- obj@rowdata
  r <- if (missing(subset))
    rep_len(TRUE, nrow(x))
  else {
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r))
      stop("'subset' must be logical")
    r & !is.na(r)
  }
  obj@data@data <- purrr::map(obj@data@data,~.x[r,])
  obj@rowdata <- obj@rowdata[r,]
  obj
})

setMethod("subsetCol",
          signature(obj = "brtExperiment"),
function(obj, subset, ...) {
  x <- obj@coldata
  r <- if (missing(subset))
    rep_len(TRUE, nrow(x))
  else {
    e <- substitute(subset)
    r <- eval(e, x, parent.frame())
    if (!is.logical(r))
      stop("'subset' must be logical")
    r & !is.na(r)
  }
  obj@data@data <- purrr::map(obj@data@data,  ~ .x[, r])
  obj@coldata <- obj@coldata[r,]
  obj
})
setMethod("subsetCol",
          signature(obj = "bcInputome"),
  function(obj, subset,...) {
    assay <- obj@assays[[obj@active_assay]]
    assay <- subsetCol.brtdata(assay,substitute(subset))
    obj@assays[[obj@active_assay]] <- assay
    obj
  })
setMethod("subsetCol",
          signature(obj = "inputData"),
  function(obj, subset,...) {
    browser()
    x <- obj@coldata
    r <- if (missing(subset))
      rep_len(TRUE, nrow(x))
    else {
      e <- substitute(subset)
      r <- eval(e, x, parent.frame())
      if (!is.logical(r))
        stop("'subset' must be logical")
      r & !is.na(r)
    }
    obj@data <- purrr::map(obj@data,  ~ .x[, r])
    obj@coldata <- obj@coldata[r,]
    obj
  })
subsetCol.brtdata <- function(obj, condition, ...) {
  x <- obj@coldata
  r <- if (missing(condition))
    rep_len(TRUE, nrow(x))
  else {
    r <- eval(condition, x, parent.frame())
    if (!is.logical(r))
      stop("'subset' must be logical")
    r & !is.na(r)
  }
  obj@data <- purrr::map(obj@data,  ~ .x[, r])
  obj@coldata <- obj@coldata[r, ]
  obj
}
# ======================== Matrix =============================================
.sortMatrix <- function(x, by) {
  if (by == "row") {
    idx <- Matrix::rowMeans(x) %>%
      sort(decreasing = T, index.return = T)
    x <- x[idx$ix,]
  } else if (by == "col") {
    idx <- Matrix::colMeans(x) %>%
      sort(decreasing = T, index.return = T)
    x <- x[, idx$ix]
  }
}

df2mat <- function(x,key,value,rowname) {
  if (! "data.frame" %in%  class(x))
    stop("input should be a data.frame")
  x %>%
    ungroup %>%
    tidyr::spread(key = key,value = value,fill = 0) %>%
    tibble::column_to_rownames(rowname) %>%
    as.matrix()
}

mat2df <- function(x,rowname.var,colname.var,value.var) {
  x.df <- x %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "rowname") %>%
    tidyr::gather(key = colname,value = value,-rowname)
  idx <- which(colnames(x.df) == "rowname")
  colnames(x.df)[idx] <- rowname.var
  idx <- which(colnames(x.df) == "colname")
  colnames(x.df)[idx] <- colname.var
  idx <- which(colnames(x.df) == "value")
  colnames(x.df)[idx] <- value.var
  x.df
}

shuffleMat <- function(data) {
  idx <- 1:length(data)
  idx.new <- sample(idx,length(idx),replace = F)
  data.new <- data[idx.new]
  data.new <- matrix(data.new,
                     nrow = nrow(data),
                     ncol = ncol(data))
  dimnames(data.new) <- dimnames(data)
  data.new
}

shuffleMat_col <- function(data) {
  data.new <- apply(data,2,function(x){
    idx <- 1:length(x)
    idx.new <- sample(idx,length(idx),replace = F)
    x[idx.new]
  })
  dimnames(data.new) <- dimnames(data)
  data.new
}
