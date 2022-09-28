#' subsetRow
#' @export
setGeneric("subsetRow",function(obj,subset,...){
  standardGeneric("subsetRow")
})

#' subsetCol
#' @export
setGeneric("subsetCol",function(obj,subset,...){
  standardGeneric("subsetCol")
})

#' brtVlnPlot
#' @export
setGeneric('brtVlnPlot',
           function(obj,
                    feature,
                    logscale = TRUE,
                    title = NA,
                    ...)
             standardGeneric("brtVlnPlot"))
#' brtMergeRegion
#' @param reference a data.frame describe the previous region name
#' (colname old_name), and the corresponding new region name(colname new_name).
#' @export
setGeneric("brtMergeRegion",
           function(obj,
                    reference)
             standardGeneric("brtMergeRegion"))
