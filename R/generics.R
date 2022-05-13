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
