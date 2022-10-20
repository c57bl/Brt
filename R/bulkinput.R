

#===============================================================================
#                             clustering bulk input
# ==============================================================================

# ========================= calculate distance ================================
brtDistScInput <- function(obj,assay = "raw",focus = "data",method = "cosin"){
  data <- obj@bulk@assays[[assay]]@data[[focus]]
  data.dist <- brtDist.matrix(data,method = method)
  obj@bulk@clusters@dist[[paste(assay,focus,method,sep = ".")]] <- data.dist
  obj
}
# calculate pairwise distance of each row
brtDist.matrix <- function(x,method){
  x <- as.matrix(x)
  if (method %in% c("euclidean", "maximum", "manhattan",
                    "canberra", "binary", "minkowski")){
    mat2 <- dist(x, method) %>%
      as.matrix()
    return(mat2)
  }
  if (is.matrix(x)) {
    if (nrow(x) < 2) {
      stop("`x` should have at least two rows.")
    }
    nr = nrow(x)
    mat2 = matrix(NA, nrow = nr, ncol = nr)
    rownames(mat2) = colnames(mat2) = rownames(x)
    for (i in 2:nr) {
      for (j in 1:(nr - 1)) {
        mat2[i, j] = brtDist.pairs(x[i,], x[j,], method = method)
      }
    }
    as.matrix(as.dist(mat2))
  } else {
    stop("`x` can be a matrix.")
  }
}

brtDist.pairs <- function(x,y,method){
  if (method == "cosin"){
    1 - sum(x*y)/sqrt(sum(x^2)*sum(y^2))
  } else {
    stop("invalid method")
  }
}
# ======================== cluster  ============================================
brtClusterScInput.search <- function(obj,focus,method = "kmeans",k.max = 20){
  focus.avaiable <- names(obj@bulk@clusters@dist)
  if (! focus %in% focus.avaiable) {
    stop("invalid focus \n", "avaiable focus: ",focus.avaiable)
  }
  data <- obj@bulk@clusters@dist[[focus]]
  if (method == "kmeans") {
    fun <- stats::kmeans
  } else if (method == "dbscan") {
    fun <- dbscan::dbscan
  } else if (method == "pam") {
    fun <- cluster::pam
  } else if (method == "cutree") {
    }else {
    stop("invalid method")
  }
  stat_method <- c("wss", "gap_stat")
  results <- lapply(stat_method, function(x) {
    factoextra::fviz_nbclust(data, fun, method = x, k.max = k.max,)
  })
  names(results) <- stat_method
  return(results)
}

brtClusterScInput <- function(obj,focus,k,method = "kmeans",n = 1,slot = "dist"){
  focus.avaiable <- names(obj@bulk@clusters@dist)
  if (!focus %in% focus.avaiable) {
    stop("invalid focus \n", "avaiable focus: ", focus.avaiable)
  }
  data <- obj@bulk@clusters@dist[[focus]]
  cluster.name <- paste(focus, method, k, sep = ".")
  if (method == "kmeans") {
    cluster <- stats::kmeans(data, k)
    ident <- cluster$cluster
  } else if (method == "pam") {
    cluster <- cluster::pam(data, k)
    ident <- cluster$clustering
  } else if (method == "cutree") {
    hc <- hclust(as.dist(data))
    cluster <- stats::cutree(hc, k)
    ident <- cluster
  } else {
    stop("invalid method")
  }
  # rank cluster by cluster size
  ident.rank <- names(sort(table(ident),decreasing = T))
  ident <- lapply(1:length(ident.rank),function(x){
    idx <- which(ident == as.character(ident.rank[x]))
    new <- ident[idx]
    new[1:length(new)] <- x
    new
  }) %>% unlist()
  ident <- factor(ident,levels = as.character(1:length(ident.rank)))
  # update to rowdata
  idx.ident <- match(obj@bulk@rowdata$name,names(ident))
  obj@bulk@rowdata[[cluster.name]] <- ident[idx.ident]
  obj@bulk@clusters@objs[[cluster.name]] <- cluster
  obj
}
# ============================= visualize clusters ============================
brtClusterScInput.plot <- function(obj,focus,reduction) {
  focus.avaiable <- setdiff(colnames(obj@bulk@rowdata), "name")
  if (!focus %in% focus.avaiable) {
    stop("invalid split \n", "avaiable split: ", focus.avaiable)
  }
  data <- left_join(obj@bulk@rowdata,
                    obj@bulk@clusters@reduction[[reduction]],
                    by = "name") %>%
          rename("ident" = focus)
  brtClusterScInput.plot.dotplot(data)
}
brtClusterScInput.splitPlot <- function(obj,
                                        splitby,
                                        type = "hm",
                                        assay = "raw",
                                        focus = "data") {
  split.avaiable <- setdiff(colnames(obj@bulk@rowdata), "name")
  if (!splitby %in% split.avaiable) {
    stop("invalid split \n", "avaiable split: ", split.avaiable)
  }
  idents <- obj@bulk@rowdata[[splitby]]
  ident.unique <- levels(obj@bulk@rowdata[[splitby]])
  data <- obj@bulk@assays[[assay]]@data[[focus]]
  data.split <- lapply(ident.unique, function(x) {
    idx <- which(idents == x)
    data[idx,]
  })
  names(data.split) <- ident.unique
  hm.split <- lapply(ident.unique, function(x) {
    data <- data.split[[x]]
    data <- as.matrix(data)
    title <-
      paste("celltype", as.character(x), "with", nrow(data), "cells")
    if (type == "hm") {
      Heatmap(
        data,
        cluster_columns = F,
        show_row_names = F,
        name = focus,
        row_title = title
      )
    } else if (type == "line") {
      as_tibble(data) %>%
        mutate(name = rownames(data)) %>%
        tidyr::gather(key = "region","value",-name) %>%
        ggplot(aes(region,value,group = name)) +
        geom_line(alpha = 0.3) +
        theme_classic() +
        labs(title = title)
    }
  })
  names(hm.split) <- ident.unique
  hm.split
}

# subfunction of brtClusterScInput.plot
brtClusterScInput.plot.dotplot <- function(x){
  idx.remove <- which(x$ident == "Others")
  x <- x[-idx.remove,]
  # calculate label axis
  x %>%
    group_by(ident) %>%
    summarise(label.x = mean(dim1),label.y = mean(dim2)) -> x.label
  # draw all
  ggplot(x, aes(dim1, dim2, color = ident)) +
    geom_point() +
    theme_classic() +
    ggrepel::geom_text_repel(
      data = x.label,
      mapping = aes(label.x, label.y, label = ident),
      color = "black",
      size = 3
    ) +
    stat_ellipse(alpha = 0.5) +
    theme_classic() +
    theme(legend.position = "none")
}
# ========================== dimension reduction ===============================
brtClusterScInput.tsne <- function(obj,assay = "raw", focus = "data") {
  reduce.name <- paste(assay,focus,"tsne",sep = ".")
  data <- as.matrix(obj@bulk@assays[[assay]]@data[[focus]])
  data.reduce <- Rtsne::Rtsne(data,check_duplicates = FALSE)
  data.reduce <- as.data.frame(data.reduce$Y)
  colnames(data.reduce) <- c("dim1","dim2")
  data.reduce$name <- rownames(data)
  obj@bulk@clusters@reduction[[reduce.name]] <- data.reduce
  obj
}
brtClusterScInput.pca <- function(obj,assay = "raw", focus = "data") {
  reduce.name <- paste(assay,focus,"pca",sep = ".")
  data <- as.matrix(obj@bulk@assays[[assay]]@data[[focus]])
  data.reduce <- stats::prcomp(data)
  data.reduce <- as.data.frame(data.reduce$x[,c(1:2)])
  colnames(data.reduce) <- c("dim1","dim2")
  data.reduce$name <- rownames(data)
  obj@bulk@clusters@reduction[[reduce.name]] <- data.reduce
  obj
}
brtClusterScInput.umap <- function(obj,assay = "raw", focus = "data") {
  reduce.name <- paste(assay,focus,"umap",sep = ".")
  data <- as.matrix(obj@bulk@assays[[assay]]@data[[focus]])
  data.reduce <- umap::umap(data)
  data.reduce <- as.data.frame(data.reduce$layout)
  colnames(data.reduce) <- c("dim1","dim2")
  data.reduce$name <- rownames(data)
  obj@bulk@clusters@reduction[[reduce.name]] <- data.reduce
  obj
}

# ======================== assign cluster ======================================
brtClusterScInput.active <- function(obj,active.ident) {
  rowdata <- obj@bulk@rowdata
  if (active.ident %in% colnames(rowdata)) {
    rowdata$ctype <- rowdata[[active.ident]]
  } else
    stop("invalid active.ident")
  obj@bulk@rowdata <- rowdata
  obj
}
brtClusterScInput.pull <- function(obj, name) {
  metadatas <- obj@metadata
  metadata.select <- lapply(obj@metadata, function(x) {
    x[c("name", name)]
  }) %>%
    do.call(rbind, .)
  obj@bulk@rowdata <-
    left_update(obj@bulk@rowdata, metadata.select, by = "name")
  obj
}
# newname, value is old, names of value is newname
brtClusterScInput.rename <- function(obj, newname) {
  rowdata <- obj@bulk@rowdata
  oldname <- rowdata$ctype
  if (! length(newname) == length(levels(oldname)))
    stop("the length of names are not equal")
  newname <- newname[match(levels(oldname),newname)]
  levels(oldname) <- names(newname)
  rowdata$ctype <- oldname
  obj@bulk@rowdata <- rowdata
  obj
}

# ========================= stat ===============================================
brtClusterScInput.binomial <- function(obj,val1,val2) {
  rowdata <- obj@bulk@rowdata
  rowdata["val1"] <- rowdata[[val1]]
  rowdata["val2"] <- rowdata[[val2]]
  data.df <- rowdata %>%
    group_by(val1, val2) %>%
    summarise(number = n(), .groups = "keep") %>%
    ungroup()
  data.mat <- df2mat(data.df,"val1","number","val2")
  data.binomial <- binomial2d(data.mat)
  # control labels
  direction.mat <- data.mat - data.binomial$E
  label.mat <- matrix("",ncol = ncol(data.mat),nrow = nrow(data.mat))
  label.mat[which(direction.mat > 0 & data.binomial$p.value < 0.05)] <- ">"
  label.mat[which(direction.mat < 0 & data.binomial$p.value < 0.05)] <- "<"
  dimnames(label.mat) <- dimnames(data.mat)
  label.mat <- mat2df(label.mat,rowname.var = "val2",
                      colname.var = "val1",
                      value.var = "label")

  # regenerate data.df to fill nan
  data.df <- mat2df(data.mat,rowname.var = "val2",
                    colname.var = "val1",
                    value.var = "counts")
  # control labels
  data.df <- left_update(data.df,label.mat,by = c("val1","val2"))
  data.df$label.count <- data.df$counts
  data.df$label.count[which(data.df$counts == 0)] <- ""
  data.df$val2 <- factor(data.df$val2,levels = levels(rowdata$val2))
  # draw
  p1 <- brtClusterScInput.binomial.plot(data.df,val1,val2)
  list(data = data.mat,
       test = data.binomial,
       plot = p1)
  }

brtClusterScInput.binomial.plot <- function(x,val1,val2) {
  ggplot(x,aes(val1,val2,fill = counts,label = label.count)) +
    geom_tile() +
    geom_text(size = 3.5,vjust = 1) +
    scale_fill_viridis_c() +
    geom_text(aes(label = label),angle = 90,hjust = -0.2,color = "gray") +
    xlab(val1) +
    ylab(val2) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90))
}


# ==============================================================================
#                        bulk input intrinsic correlation
# ==============================================================================


# ============================= bin analysis ===================================

brtBinScInput.cp <- function(obj,assay = "raw",focus = "data_bin",R = 1000) {
  data <- obj@bulk@assays[[assay]]@data[[focus]]%>%
    as.matrix()
  cp <- colWiseCp(data)
  # shuffling and cp
  cp.resample <- purrr::map(1:R, function(r){
    data.new <- shuffleMat_col(data)
    colWiseCp(data.new)
  }) %>%
    do.call(rbind,.) %>%
    as_tibble()
  # calculate p value
  cp <- purrr::map2(cp$V1, cp$V2, function(x, y) {
    idx1 <- which(cp.resample$V1 == x & cp.resample$V2 == y)
    data.resample <- cp.resample[idx1, ]$cp
    data.resample <- sort(data.resample)
    data.resample.range <- quantile(data.resample, c(0.25, 0.975))
    idx2 <- which(cp$V1 == x & cp$V2 == y)
    data.observe <- cp[idx2, ]$cp
    if (x == y) {
      p.value <- 1
    } else {
      # greater than
      idx.greater <- which(data.resample > data.observe)
      p.greater <- length(idx.greater) / R
      p.less <- 1 - p.greater
      p.value <- min(p.greater, p.less)
      p.value[p.value == 0] <- 1 / R
    }
    data.frame(
      V1 = x,
      V2 = y,
      p.value = p.value,
      observe = data.observe,
      left  = data.resample.range[1],
      right = data.resample.range[2]
    )
  }) %>%
    do.call(rbind, .)
  cp$p.adj <- p.adjust(cp$p.value,method = "fdr")
  cp$sig <- ""
  cp$sig[which(cp$observe > cp$right & cp$p.adj < 0.05)] <- ">"
  cp$sig[which(cp$observe < cp$left & cp$p.adj < 0.05)] <- "<"
  cp$V1 <- factor(cp$V1,levels = colnames(data))
  cp$V2 <- factor(cp$V2,levels = colnames(data))
  # plot
  plot.cp <- brtBinScInput.cp.plot(cp)
  # return
  list(results = cp,
       plot = plot.cp)
}

brtBinScInput.cp.plot <- function(x) {
  x %>%
    ggplot(aes(V1,V2,fill = observe,label = sig)) +
    geom_tile() +
    geom_text(angle = 90) +
    scale_fill_gradient2(low = "blue",
                         mid = "white",
                         high = "red",
                         midpoint = 0.5,name = "P(A|B)") +
    xlab("B") +
    ylab("A") +
    theme_minimal()
}

# binomial null hypothesis
brtBinScInput.binor <- function(obj,assay = "raw", focus = "data_bin",n,
                                p.equal = F,exclude = F,
                                upset.threshold = 0.1) {
  data <- obj@bulk@assays[[assay]]@data[[focus]]
  data <- as.matrix(data)
  n.observe <- countRowCombine(data,n)
  n.predict <- binomialRowCombine(data,n,p.equal = p.equal,exclude = exclude)
  n.test.tb <- left_update(n.predict$pn,n.observe,by = "combine",treat.na = 0)
  counts.predict <-
    purrr::map2(n.test.tb$p, n.test.tb$counts, function(a, b) {
      result <- binom.test(b, n = n.predict$nrow.adj, p = a)
      data.frame(
        p.value = result$p.value,
        estimated = a * n.predict$nrow.adj
      )
    }) %>%
    do.call(rbind, .)
  n.test.tb <- cbind(n.test.tb,counts.predict)
  n.test.tb <- n.test.tb %>%
    mutate(p.adj = p.adjust(p.value, method = "fdr"),
           diff = counts - estimated)
  # draw
  plot.upset <- brtBinScInput.binor.upsetplot(n.test.tb,upset.threshold)
  plot.vol <- brtBinScInput.binor.volplot(n.test.tb,esti.threshold = 0.2)
  data.lineplot <- obj@bulk@assays[[assay]]@data$data_proportion
  plot.line <- brtBinscInpput.binor.lineplot(n.test.tb,
                                             data,
                                             data.lineplot)
  list(result = n.test.tb,
       plot.upset = plot.upset,
       plot.vol = plot.vol,
       plot.line = plot.line)

}
brtBinScInput.binor.upsetplot <- function(x,upset.threshold) {
  x <- x %>%
    as_tibble %>%
    filter(p.adj < upset.threshold) %>%
    tidyr::gather("group", "number", counts, estimated)
  x$label <- ""
  x$label[which(x$p.adj < 0.05 &
                  x$group == "counts")] <- "*"
  x %>%
    mutate(combine = stringr::str_split(combine, stringr::fixed("+"))) %>%
    ggplot(aes(
      combine,
      number,
      group = group,
      fill = group,
      label = label
    )) +
    geom_col(position = position_dodge2(width = 1)) +
    geom_text(position = position_dodge2(width = 1), vjust = -0.2) +
    ggupset::scale_x_mergelist(sep = "&") +
    ggupset::axis_combmatrix(sep = "&") +
    scale_fill_manual(values = c("#FC4E07", "gray")) +
    theme_classic()

}
brtBinScInput.binor.volplot <- function(x,esti.threshold = 0.2) {
  # label control
  x$label <- x$combine
  x$label[which(x$p.adj > 0.05)] <- ""
  # color control
  x$color <- "red"
  x$color[which(x$counts < x$estimated)] <- "blue"
  x$color[which(x$p.adj > 0.05)] <- "gray"
  x <- x %>%
    filter(counts > 0, estimated > esti.threshold) %>%
    mutate(effect = log2(counts / estimated),
           p.value = - log10(p.value))
  x %>%
    ggplot(aes(effect, p.value, label = label)) +
    geom_point(color = x$color) +
    theme_classic() +
    ggrepel::geom_text_repel(size = 2.5) +
    xlab("log2(observed/expected)") +
    ylab("-log10(p.value)")
}
brtBinscInpput.binor.lineplot <- function(x, y, z) {
  x <- filter(x, p.adj < 0.05)
  if (nrow(x) == 0) {
    warning("lineplot: none of the p.adj passed")
    return(NULL)
  }
  groups <- x$combine %>%
    stringr::str_split(stringr::fixed("+"))
  idxs <- lapply(groups, function(g) {
    idx.colname <- match(g, colnames(y))
    y.sub <- y[, idx.colname]
    which(rowSums(y.sub) == length(g))
  })
  lapply(1:length(idxs), function(i) {
    idx.i <- idxs[[i]]
    if (length(idx.i) > 0) {
    name.i <- paste(groups[[i]], sep = "", collapse = "+")
    value.i <- paste("Observed",x[i,]$counts,"Expected",round(x[i,]$estimated))
    data.i <- z[idx.i, ] %>%
      as.matrix() %>%
      as_tibble() %>%
      mutate(idx = as.character(1:nrow(.))) %>%
      tidyr::gather("region", "proportion", -idx) %>%
      mutate(region = factor(region, levels = colnames(y)))
    data.i$color <- "gray"
    data.i$color[which(data.i$region %in% groups[[i]])] <- "red"
    data.i %>%
      ggplot(aes(region, proportion, group = idx,color = color)) +
      geom_line(alpha = 0.6,color = "black") +
      geom_point() +
      scale_color_manual(values = c("gray","red")) +
      ggtitle(name.i) +
      labs(subtitle = value.i) +
      theme_classic() +
      theme(legend.position = "none")
    }
  }) -> plots
  names(plots) <- x$combine
  plots
}

# ==============================================================================
#                           merge cells to celltype
# ==============================================================================
brtMergeCell.scInput <- function(obj,var,assay = "raw",focus = "data"){
  data <- obj@bulk@assays[[assay]]@data[[focus]]
  merge.by <- obj@bulk@rowdata[[var]]
  groups <- sort(unique(merge.by))
  group.n <- lapply(groups,function(x) {
    length(which(merge.by == x))
  }) %>%
    do.call(rbind,.)
  merge.data <- lapply(groups,function(x){
    idx <- which(merge.by == x)
    if (length(idx) > 1) {
      colSums(data[idx,])
    } else {
      data[idx,]
    }
  }) %>%
    do.call(rbind,.)
  rownames(merge.data) <- groups
  list(group.n = group.n,
       merge.data = merge.data)
}
