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
