#' brtPlotNoise
#' @import ggplot2
#' @param obj an brtstarter object, the idents of coldata should contain Neurons
#' and at least one non-neuronal ident.
#' @param feature bc or a name of numeric column of coldata. bc will override bc
#' column if bc column exist.
#' @param background a character vector, describe which non-neuronal cluster used
#' to calculate the threshold. Using all non-neuronal groups if none is provided.
#' @param threshold a numeric value, default to be 3
brtPlotNoise <- function(obj,feature,background = NULL, threshold = 3,
                         signal = "Neurons"){
  # prepare data.frame with colnames source and counts,source contains
  # signal and background only
  if (length(feature) > 1) stop("multiple features are not supported")
  if (feature == "bc") {
    result <- getNoiseBc(obj,background,threshold,signal)
  } else {
    result <- getNoiseFeature(obj,feature,background,threshold,signal)
  }
  # for log visualize, 0 add to 0.5
  data <- result$data
  data$counts[data$counts == 0] <- 0.5
  ggplot(data,aes(source,counts)) +
    geom_jitter(size = 1, color = "gray") +
    geom_violin() +
    geom_hline(yintercept = result$cutoff + 0.4,
               color = "red",size = 1) +
    scale_y_log10() +
    ggtitle(feature,
            subtitle = paste0("cutoff: ",result$cutoff)) +
    theme_classic()


}

# sub function of brtPlotNoise
# reutrn a list
getNoiseFeature <- function(obj,feature,background = NULL, threshold = 3,
                            signal = "Neurons"){
  if (feature %in% colnames(obj@coldata)) {
    if (is.null(background)){
      background <- setdiff(unique(obj@coldata$ident),signal)
    }
    data <- select(obj@coldata,ident,feature) %>%
      rename(source = ident,counts= feature) %>%
      filter(source %in% c(signal,background)) %>%
      mutate(source = as.character(source))
    data$splitsource <- data$source
    data$source[data$source != signal] <- "Background"
    counts.n <- subset(data, source == signal)$counts
    counts.nn <- subset(data, source == "Background")$counts
    cutoff <- median(counts.nn) + threshold * IQR (counts.nn)
    result <- list(cutoff = cutoff,
                   data = data)
  } else {
    stop("can not find feature: ", feature)
  }
  return(result)
}

#' brtPlotInfection
#' plot the results of brtGetInfection
#' @description This function probably plot the most important features of
#' the result of brtGetInfection, especially the SP value and the correlation
#' between tva counts and bc+ proportion
#' @seealso [brtGetInfection()] for details
#' @export
brtPlotInfection <- function(obj,tva.c,rv.c){
  result <- brtGetInfection(obj,tva.c,rv.c)
  # draw bc+
  p <- c(result$record$RP,result$RPT)
  counts <- c(0:tva.c,paste0(tva.c, "+"))
  counts <- factor(counts,levels = counts)
  data.frame(p = p,
             counts = counts) %>%
    ggplot(aes(counts,p)) +
    geom_col() +
    theme_light() +
    ggtitle("bc+ proportion") -> p1
  # draw estimated
  counts <- c(1:tva.c,paste0(tva.c, "+"))
  counts <- factor(counts, levels = counts)
  p <- c(result$record$IP[-1],result$IE)
  data.frame(p = p,
             counts = counts) %>%
    ggplot(aes(counts,p)) +
    geom_col() +
    theme_light() +
    ggtitle("Estimated initially infeated proportion",
            subtitle = paste0("SP: ",sprintf("%.3f",result$SP))) -> p2
  cowplot::plot_grid(p1,p2,nrow = 2)
}

#' brtPlotUniqueBc
#' plot the results of brtSetUniqueBc
#' @param obj a brtstaretr obj
#' @param nrow control the row number of plots
#' @seealso [brtSetUniqueBc()] for details
#' @export
brtPlotUniqueBc <- function(obj,nrow = 2){
  if (! "unique.bc" %in% colnames(obj@rowdata))
    stop("please ruun brtSetUniqueBc")
  obj@rowdata %>%
    filter(!is.na(cell.p)) %>%
    select(cell.p) %>%
    arrange(cell.p) %>%
    mutate(id = 1:length(cell.p)) %>%
    ggplot(aes(id,cell.p)) +
    geom_point(size = 1) +
    ylab("max(BC) / sum(BC)") +
    xlab("Unique virus in staretrs") +
    ggtitle("Normalize across cells") +
    theme_classic() -> p1
  obj@rowdata %>%
    filter(!is.na(cell.u)) %>%
    select(cell.u,cell) %>%
    group_by(cell) %>%
    summarise(cell.u = sum(cell.u)) %>%
    arrange(cell.u) %>%
    mutate(id = 1:length(cell)) %>%
    ggplot(aes(id,cell.u)) +
    geom_point(size = 1) +
    ylab("sum(unique BC) / sum(BC)") +
    xlab("Starters with unique bc") +
    ggtitle("Normalize within each cell") +
    theme_classic() -> p2
  cowplot::plot_grid(p1,p2, nrow = nrow)
}

# vlnplot
setMethod("brtVlnPlot", "brtStarter", function(obj,
                                               feature,
                                               logscale,
                                               title,
                                               size = 0.1 ,
                                               alpha = 0.5,
                                               color = FALSE,
                                               ...) {
  coldata <- obj@coldata
  # select feature
  if (class(feature) == 'character') {
    colordata = coldata[[feature]]
  } else
    colordata <- feature
  coldata %>%
    select(ident, feature) %>%
    rename(data = feature) -> coldata
  # title
  if (is.na(title) & class(feature) == 'character')
    mytitle = feature
  else if (is.na(title) & !class(feature) == 'character')
    mytitle = NULL
  else
    mytitle = title
  if (logscale) coldata$data[coldata$data == 0] <- 0.5
  brtVlnPlot(
    obj = coldata,
    logscale = logscale,
    title = mytitle,
    size = size,
    alpha = alpha,
    color = color
  )
})

setMethod("brtVlnPlot", "data.frame", function(obj,
                                               feature,
                                               logscale,
                                               title = NULL,
                                               size = 0.1 ,
                                               alpha = 0.5,
                                               color = color,
                                               ...) {
  coldata <- obj
  coldata$ax1 <- coldata$ident
  coldata$ax2 <- coldata$data
  if (color) plot_map <- ggplot(data = coldata, aes(ax1, ax2,fill = ident))
  else plot_map <- ggplot(data = coldata, aes(ax1, ax2))
  plot_vln <- geom_violin()
  plot_jitter <- geom_jitter(size = size, alpha = alpha,color = "black")
  plot_log <- scale_y_log10()
  plot_title <- labs(title = title)
  plot_theme <- theme_classic() +
    theme(axis.text.x = element_text(angle = 45),legend.position = "none")
  plotfinal <- plot_map +
    plot_vln +
    plot_title +
    plot_jitter +
    xlab(NULL) +
    ylab("counts") +
    plot_theme

  if (logscale)
    plotfinal <- plotfinal + plot_log
  return(plotfinal)
})
