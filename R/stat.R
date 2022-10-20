# =========================== matrix binomial ====================================
binomial2d <- function(data) {
  N <- Matrix::colSums(data)
  M <- Matrix::rowSums(data)
  S <- sum(data)
  P <-
    matrix(rep(N / S, nrow(data)), nrow = nrow(data), byrow = T) *
    matrix(rep(M / S, ncol(data)), nrow = nrow(data))
  Q <- 1 - P
  E <- S * P
  V <- (P * Q * S) ^ 0.5
  MP <- S * P > 4
  MQ <- S * Q > 4
  C <- matrix(FALSE, nrow = nrow(data), ncol = ncol(data))
  C[which(MP == TRUE & MQ == TRUE)] <- TRUE
  purrr::pmap(list(data, E, V, C, P), function(a, b, c, d, e) {
    if (d) {
      if (a < b)
        p <- 1 - stats::pnorm(abs((a - 0.5 - b) / c))
      else
        p <- 1 - stats::pnorm(abs((a + 0.5 - b) / c))
    }
    else
      p <- choose(S, a) * e ^ a * (1 - e) ^ (S - a)
    return(p)
  }) %>%
    unlist %>%
    matrix(nrow = nrow(data)) -> p.raw
  dimnames(p.raw) <- dimnames(data)
  p.raw %>%
    stats::p.adjust(method = "fdr") %>%
    matrix(nrow = nrow(data),
           ncol = ncol(data)) -> p.adj
  dimnames(p.adj) <- dimnames(data)
  dimnames(E) <- dimnames(data)
  return(list(
    E = E,
    p.value = p.raw,
    p.adj = p.adj
  ))
}

# adjust row number of a matrix, take column possibility into consideration
adjustRowNumber <- function(x) {
  I <- colSums(x)
  No <- nrow(x)
  P <- I / nrow(x)
  f <- function(Ni) Ni*(1-prod(1-I/Ni))-No
  Ni <- uniroot(f, c(No,No*5),tol = 0.1)$root
  Ni <- round(Ni)
}

binomialRowCombine <- function(x,n,p.equal = FALSE,exclude = FALSE) {
  x.nrow.adj <- adjustRowNumber(x)
  x.ncol <- ncol(x)
  if (p.equal) {
    p <- 1 - (1 - nrow(x) / x.nrow.adj) ^ (1 / x.ncol)
    p <- rep(p, ncol(x))
  } else {
    p <- colSums(x) / x.nrow.adj
  }
  pn <- lapply(n, function(i) {
    group.design <- combn(x.ncol, i)
    group.p <- lapply(1:ncol(group.design), function(ii) {
      group.curr <- group.design[, ii]
      if (exclude) {
        group.p.curr <- prod(p[group.curr], (1 - p)[-group.curr])
      } else {
        group.p.curr <- prod(p[group.curr])
      }
    }) %>%
      unlist
    group.name <- lapply(1:ncol(group.design), function(ii) {
      group.curr <- group.design[, ii]
      group.name.curr <- colnames(x)[group.curr]
      paste(group.name.curr, sep = "", collapse = "+")
    }) %>%
      unlist
    data.frame(combine = group.name,
               p = group.p)
  }) %>%
    do.call(rbind, .)
  list(pn = pn,
       nrow.adj = x.nrow.adj)
}

countRowCombine <- function(x, n) {
  x.ncol <- ncol(x)
  cn <- lapply(n, function(i) {
    group.design <- combn(x.ncol, i)
    group.c <- lapply(1:ncol(group.design), function(ii) {
      group.curr <- group.design[, ii]
      if (length(group.curr) > 1) {
        group.c.curr <- length(which(rowSums(x[, group.curr]) == i))
      } else {
        group.c.curr <- length(which(x[, group.curr] == i))
      }
    }) %>%
      unlist
    group.name <- lapply(1:ncol(group.design), function(ii) {
      group.curr <- group.design[, ii]
      group.name.curr <- colnames(x)[group.curr]
      paste(group.name.curr, sep = "", collapse = "+")
    }) %>%
      unlist
    data.frame(combine = group.name,
               counts = group.c)
  }) %>%
    do.call(rbind, .)
  cn
}
# without consider column probability
binomialRowCombine2 <- function(x,n) {
  x.nrow.adj <- adjustRowNumber(x)
  x.ncol <- ncol(x)
  p <- 1 - (1 - nrow(x) / x.nrow.adj) ^ (1 / x.ncol)
  pn <- lapply(n, function(i) {
   p^i * (1-p)^(x.ncol - i)
  }) %>%
    unlist
  names(pn) <- n
  list(p = pn,
       n = x.nrow.adj)
}


# ======================= column wised correlation =============================

colWiseCor <- function(data, method) {
  if (is.null(colnames(data)))
    stop("colname of data is required")
  comp.design <- combn(ncol(data), 2)
  data.cor <-
    purrr::map2(comp.design[1,], comp.design[2,], function(x, y) {
      cor.result <- cor.test(data.scaled[, x],
                             data.scaled[, y],
                             method = method,
                             exact = FALSE)
      data.frame(
        cor = cor.result$estimate,
        p.value = cor.result$p.value,
        row.names = NULL
      )
    }) %>%
    do.call(rbind, .) %>%
    cbind(t(comp.design), .) %>%
    data.frame() %>%
    mutate(
      X1 = factor(colnames(data)[X1], levels = rev(colnames(data))),
      X2 = factor(colnames(data)[X2], levels = rev(colnames(data))),
      p.adj = p.adjust(p = p.value, method = "fdr")
    )
  data.cor$sig <- ""
  data.cor$sig[which(data.cor$p.adj < 0.05)] <- "*"
  data.cor
}

# conditional probability
calcp <- function(a,b){
  aub <- which((a+b) > 0) %>%
    length()
  ab <- sum(a*b)
  pab <- ab / aub
  pa <- sum(a) / aub
  pb_a <-  pab / pa
  return(pb_a)
}

colWiseCp <- function(data) {
  group.design <- combn(colnames(data),2) %>%
    t() %>%
    as.data.frame()
  group.design <- rbind(group.design,
                        cbind(group.design$V2,group.design$V1),
                        cbind(colnames(data),colnames(data)))
  purrr::map2(group.design$V1,group.design$V2,function(x,y){
    idx.x <- which(colnames(data) == x)
    idx.y <- which(colnames(data) == y)
    calcp(data[,idx.x],data[,idx.y])
  }) %>%
    unlist() -> cp
  cp <- data.frame(group.design,cp = cp)
  cp
}

