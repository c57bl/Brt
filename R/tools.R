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
