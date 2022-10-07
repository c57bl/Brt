
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
                    target,
                    focus = "data")
             standardGeneric("brtMergeRegion"))
#' brtSetNoiseBc
#' Estimate the noise level of bcs
#' @description The data is noisy, it is important to identify noise. RV seldom
#' infect non-neuronal cells from neurons, thus the bcs expressing levels of
#' non-neuronal groups reflects the noise level of experiment. By measure the
#' median + threshold * IQR, the threshold of signal can be determined.
#' @param obj an brtstarter/bcInputome object
#' @param background a character vector, describe which non-neuronal cluster used
#' to calculate the threshold. Using all non-neuronal groups if none is provided.
#' @param threshold a numeric value, default to be 3
#' @param signal a character, describe which ident is signal, default is Neurons
#' @export
setGeneric("brtSetNoiseBc",
           function(obj,background = NULL,threshold = 3,
                    signal = "Neurons")
             standardGeneric("brtSetNoiseBc"))
#' brtNormalizeInput
#' @description normalize each bc by selected method
#' @param obj an bcinputome/brtunity obj
#' @param method an character describe the normalize method
#' @param assay an character describe which assay to use
#' @param focus an character describe which data to use
#' @return an obj with added slot:focus_method  in the data in
#' selected assay
#' @details 1. proportion , normalize each bc to 0-1 proportion
#' 2. scale, center and scale data by calling base::scale
setGeneric("brtNormalizeInput",
           function(obj,method,assay = "raw",
                    focus = "data")
             standardGeneric("brtNormalizeInput"))

