# ============================ match 10x bc =============================
#' brtMatch_10x_bc
#' @description brtMatch_10x_bc will match and the counts the umi grouped by
#' 10x barcode and barcodes provided from library made from 10x cdna. The reads2
#' contains the cell barcode and the reads1 contains the target gene information
#' @param rvBarcode_pattern A character or character vector
#' @param quality_cutoff a integer specify low quality string
#' @param cellBelow a integer, the read will be drooped if the number of low quality
#' string in cell barcode region exceeds it.
#' @param rvBelow a integer, the read will be drooped if the number of low quality
#' string in rvbarcode region(the region could be future modified by rvBarcode_roi)
#' exceeds it.
#' @param rvBarcode_max_mismatch a integer, the pattern match were performed by
#' biostring::vmatchpattern and the rvBarcode_max_mismatch control the max_mismatch
#' parameter
#' @param blockSize a numeric value specify the number of reads to be load in
#' recursive loading
#' @param rvBarcode_row a numeric vector to extract pure barcode from barcode region
#' @param seq_library a character,'cdna'/'10x'
#' @param ncore a numeric value specify the cores used in time consuming steps
#' @param cell_reference a character vector provide the white list of cell barcode,
#' the unmatched cell will be droped
#' @param rvBarcode_reference a character vector provide the white list of rvBarcode,
#' the unmatched barcode will be droped
#' @return a tibble with colnames cell,bc counts
#' @export
brtMatch_10x_bc <- function(path,rvBarcode_pattern,quality_cutoff=25,
                            cellBelow=2,rvBelow=2,blockSize=2e6,
                            rvBarcode_max_mismatch=0,rvBarcode_roi=NULL,
                            rvBarcode_cutoff=3,cell_cutoff=3,ncore=1,
                            rvBarcode_reference=NULL,cell_reference=NULL,
                            total = NULL,return_raw = FALSE){
  paths <- .brtmatch.parse_path(path)
  loop <- 0
  results <- list()
  message('start match')
  for (i in 1:length(paths$r1)){
    message("file: ",i)
    con1<-file(paths$r1[i])
    con2<-file(paths$r2[i])
    f1<- ShortRead::FastqStreamer(con1,blockSize)
    f2<- ShortRead::FastqStreamer(con2,blockSize)
    while(TRUE){
      reads1 <- ShortRead::yield(f1)
      reads2 <- ShortRead::yield(f2)
      if (length(reads1) == 0) break
      loop <- loop+1
      if (! is.null(total)) if(loop > total) break
      main <- function(){
        curr_cellbarcode <- Biostrings::subseq(reads2@sread,start=1,end=16)
        curr_umi <- Biostrings::subseq(reads2@sread,start=17,end=28)
        curr_cellqu <- Biostrings::subseq(reads2@quality@quality,start=1,end=16)
        curr_rv <- .brtmatch.matched(reads1,
                                     p = rvBarcode_pattern,
                                     m = rvBarcode_max_mismatch)
        quality_pass.cell <- .brtmatch.pass_quality_idx(x = curr_cellqu,
                                                        q = quality_cutoff,
                                                        n = cellBelow)
        quality_pass.rv <- .brtmatch.pass_quality_idx(x = curr_rv$quality,
                                                      q = quality_cutoff,
                                                      n = rvBelow,
                                                      r = rvBarcode_roi)
        idx.cell <- intersect(curr_rv$idx[quality_pass.rv],
                              quality_pass.cell)
        idx.bc <- which(curr_rv$idx %in% idx.cell)
        curr_cellbarcode[idx.cell] %>%
          as.character() -> curr_cellbarcode
        curr_umi[idx.cell] %>%
          as.character() -> curr_umi
        curr_rv$reads[idx.bc] %>%
          as.matrix() -> curr_bc
        curr_bc[,rvBarcode_roi] %>%
          apply(1,paste0,collapse='') -> curr_bc
        if (! return_raw) {
          results.sub <- .brtmatch_10x_bc.count_character(cbind(curr_cellbarcode,
                                                                curr_bc,
                                                                curr_umi),
                                                          cr = cell_reference,
                                                          rr = rvBarcode_reference)
        } else {
          idx <- which(curr_cellbarcode %in% cell_reference)
          results.sub <- cbind(curr_cellbarcode[idx],
                               curr_bc[idx],
                               curr_umi[idx])
        }
        results[[loop]] <- results.sub
        message("loops: ",loop," ",length(curr_cellbarcode),' matched')
        return(results)
      }
      tryCatch(results <- main(),
               error = function(e){
                 message("loops: ",loop," error...")
                 return(NULL)
               })
    }
    close(con1)
    close(con2)
  }
  results <- do.call(rbind,results)
  if (! return_raw) {
    results <- .brtmatch_10x_bc.count_n(results)
  }
  return(results)
}

# x, 3 column character matrix,with columnnames cellbarcode,rvbarcode,umi
# cr cellbarcode reference
# rr rvbarcode reference
.brtmatch_10x_bc.count_character <- function(x,cr = NULL,rr = NULL){
  unique(x) %>%
    as_tibble() -> x
  colnames(x) <- c("cell","bc","umi")
  if (! is.null(cr)) x <- subset(x,cell %in% cr)
  if (! is.null(rr)) x <- subset(x, bc %in% rr)
  x %>%
    group_by(cell,bc) %>%
    summarise(counts = n()) %>%
    ungroup() -> x
  return(x)
}
# x 3 column tibble with cell,bc,counts
.brtmatch_10x_bc.count_n <- function(x){
  x %>%
    group_by(cell,bc) %>%
    summarise(counts = sum(counts)) %>%
    ungroup() -> x
  return(x)
}

# ========================== match inputs ======================================
#' brtMatch_input
#' @description match a list of patterns, typically different verisons of barcodes,and
#' return a brtInput class
#' @param sampels a character or character vector describe the path of each sample
#' @param ncore a integer describe how many cores to be used
#' @param Rvpatterns a n*3 data.frame. $pattern describe the patterns of RvBarcodes,
#' $name describe the name ('rvbVn' or 'spkVn') of each pattern.
#' @param umiPattern a character vector, describe the R2 pattern of UMI
#' @return a list contain barcode  and spike-in umi counts tibble, with colnames
#' id,bc,counts,name
#' @export
brtMatch_input <- function(path,indexs,RvPatterns,umiPattern,quality_cutoff=25,
                           rvBelow=2,umiBelow=2,blockSize=2e6,rvBarcode_roi=NA,
                           umi_roi=NA,index_roi=NA,Batch=NA_character_,
                           ID=NA_character_,ncore=1,rvBarcode_max_mismatch = 0,
                           umi_max_mismatch = 0,rvBarcode_reference = NULL,
                           total = NULL,return_raw = FALSE){
  paths <- .brtmatch.parse_path(path)
  loop <- 0
  results <- list()
  message(Sys.time(),'\n',
          ' start match...')
  for (i in 1:length(paths$r1)){
    message("file: ",i)
    con1<-file(paths$r1[i])
    con2<-file(paths$r2[i])
    f1<- ShortRead::FastqStreamer(con1,blockSize)
    f2<- ShortRead::FastqStreamer(con2,blockSize)
    while(TRUE){
      loop <- loop +1
      if (! is.null(total)) if(loop > total) break
      reads1 <- ShortRead::yield(f1)
      reads2 <- ShortRead::yield(f2)
      if (length(reads1) == 0) break
      runsub <- function(r1,r2,rvp,name){
        curr_rv <- .brtmatch.matched(r1,
                                     p = rvp,
                                     m = rvBarcode_max_mismatch)
        curr_umiplus <- .brtmatch.matched(r2,
                                          p = umiPattern,
                                          m = umi_max_mismatch)
        quality_pass.rv <- .brtmatch.pass_quality_idx(x = curr_rv$quality,
                                                      q = quality_cutoff,
                                                      n = rvBelow,
                                                      r = rvBarcode_roi)
        quality_pass.umi <- .brtmatch.pass_quality_idx(x = curr_umiplus$quality,
                                                       q = quality_cutoff,
                                                       n = umiBelow,
                                                       r = umi_roi)
        idx_pass <- intersect(curr_rv$idx[quality_pass.rv],
                              curr_umiplus$idx[quality_pass.umi])
        idx.rv <- which(curr_rv$idx %in% idx_pass)
        idx.umiplus <- which(curr_umiplus$idx %in% idx_pass)
        curr_umiplus <- as.matrix(curr_umiplus$reads[idx.umiplus])
        curr_bc <- as.matrix(curr_rv$reads[idx.rv])
        curr_bc <- apply(curr_bc[,rvBarcode_roi],1,paste0,collapse='')
        curr_idx <- apply(curr_umiplus[,index_roi],1,paste0, collapse='')
        curr_umi <- apply(curr_umiplus[,umi_roi],1,paste0, collapse='')
        if (!return_raw) {
          if (!stringr::str_detect(name, "spk")) {
            results.sub <- .brtmatch_10x_bc.count_character(cbind(curr_idx,
                                                                  curr_bc,
                                                                  curr_umi),
                                                            cr = index$index,
                                                            rr = rvBarcode_reference)
            results.sub$name <- name
          } else {
            results.sub <- .brtmatch_input.count_character_umi(cbind(curr_idx,
                                                                     curr_bc,
                                                                     curr_umi),
                                                               cr = index$index)
            results.sub$name <- name
          }
          results.sub %>%
            rename(index = cell) %>%
            left_join(indexs, by = "index") %>%
            select(id, bc, counts, name) -> results.sub
        } else{
          idx <- which(curr_idx %in% indexs$index)
          results.sub <- data.frame(id = curr_idx[idx],
                                    bc = curr_bc[idx],
                                    counts = curr_umi[idx],
                                    name = name)
        }
        message(Sys.time(),'\n',
                name," loops: ",loop," ",length(idx_pass),' matched')
        return(results.sub)
      }
      tryCatch(results[[loop]] <- do.call(rbind,
                                          purrr::map2(RvPatterns$pattern,RvPatterns$name,
                                                      ~runsub(reads1,reads2,.x,.y))),
               error = function(e){
                 message(Sys.time(),"\nloops: ",loop," error, ",
                         "probably due to limited reads")
                 return(NULL)
               })
    }
    close(con1)
    close(con2)
  }
  message(Sys.time(),'\n',
          "match complete, start merge results...")
  results <- do.call(rbind,results)
  if (!return_raw) {
    results <- rename(results, cell = id)
    results <- purrr::map(RvPatterns$name, function(x) {
      y <- subset(results, name == x) %>%
        .brtmatch_10x_bc.count_n() %>%
        rename(id = cell)
    })
  } else {
    results <- purrr::map(RvPatterns$name,~ subset(results, name == .x))
  }
  names(results) <- RvPatterns$name
  return(results)
}
.brtmatch_input.count_character_umi <- function(x,cr = NULL,rr = NULL){
  unique(x) %>%
    as_tibble() -> x
  colnames(x) <- c("cell","bc","umi")
  if (! is.null(cr)) x <- subset(x,cell %in% cr)
  if (! is.null(rr)) x <- subset(x, bc %in% rr)
  x %>%
    group_by(cell) %>%
    summarise(counts = n()) %>%
    ungroup() -> x
  x$bc <- NA
  x <- select(x,cell,bc,counts)
  return(x)
}
# =========================== match 10x features ===============================
#' brtMatch_10x_target
#' @description return a dataframe record the umi counts of a given pattern
#' @export
brtMatch_10X_target <- function(path,pattern,quality_cutoff=25,
                                cellBelow=2,patternBelow=2,blockSize=2e6,
                                pattern_max_mismatch=0,cell_reference = NULL,
                                total = NULL,return_raw = FALSE){
  paths <- .brtmatch.parse_path(path)
  loop <- 0
  results <- list()
  message(Sys.time(),'\n',
          ' start match...')
  for (i in 1:length(paths$r1)){
    message("file: ",i)
    con1<-file(paths$r1[i])
    con2<-file(paths$r2[i])
    f1<- ShortRead::FastqStreamer(con1,blockSize)
    f2<- ShortRead::FastqStreamer(con2,blockSize)
    while(TRUE){
      loop <- loop +1
      if (!is.null(total)) {
        if (loop > total) break
      }
      reads1 <- ShortRead::yield(f1)
      reads2 <- ShortRead::yield(f2)
      if (length(reads1) == 0) break
      main <- function(){
        curr_cellbarcode <- Biostrings::subseq(reads2@sread,start=1,end=16)
        curr_umi <- Biostrings::subseq(reads2@sread,start=17,end=28)
        curr_cellqu <- Biostrings::subseq(reads2@quality@quality,start=1,end=16)
        curr_target <- .brtmatch.matched(reads1,
                                         p = pattern,
                                         m = pattern_max_mismatch)
        quality_pass.cell <- .brtmatch.pass_quality_idx(x = curr_cellqu,
                                                        q = quality_cutoff,
                                                        n = cellBelow)
        quality_pass.target <- .brtmatch.pass_quality_idx(x = curr_target$quality,
                                                          q = quality_cutoff,
                                                          n = patternBelow)
        idx.cell <- intersect(curr_target$idx[quality_pass.target],
                              quality_pass.cell)
        idx.target <- which(curr_target$idx %in% idx.cell)
        curr_cellbarcode[idx.cell] %>%
          as.character() -> curr_cellbarcode
        curr_umi[idx.cell] %>%
          as.character() -> curr_umi
        if (! return_raw) {
          results.sub <- .brtmatch_10x_target.count_character(cbind(curr_cellbarcode,
                                                                    curr_umi),
                                                              cr = cell_reference)
        } else {
          idx <- which(curr_cellbarcode %in% cell_reference)
          results.sub <- cbind(curr_cellbarcode[idx],
                               curr_umi[idx])
        }
        results[[loop]] <- results.sub
        message(Sys.time(),"\nloops: ",
                loop," ",length(curr_cellbarcode),' matched')
        return(results)
      }
      tryCatch(
        results <- main(),
        error = function(e) {
          message(Sys.time(),
                  "\nloops: "
                  ,
                  loop,
                  " error, ",
                  "probably due to limited reads")
          return(NULL)
        }
      )
    }
    close(con1)
    close(con2)
  }
  results <- do.call(rbind,results)
  if (! return_raw) {
    results <- results %>%
      group_by(cell) %>%
      summarise(counts = sum(counts))
  } else {
    colnames(results) <- c("cell","umi")
  }
  return(results)

}
.brtmatch_10x_target.count_character <- function(x,cr = NULL){
  unique(x) %>%
    as_tibble() -> x
  colnames(x) <- c("cell","umi")
  if (! is.null(cr)) x <- subset(x,cell %in% cr)
  x %>%
    group_by(cell) %>%
    summarise(counts = n()) %>%
    ungroup() -> x
  return(x)
}

# ========================== inner functions ===================================
.brtmatch.parse_path <- function(path){
  files_r1 <- list.files(path=path,pattern ='*[Rr]1.fastq',
                         full.names = T,recursive = T)
  files_r2 <- list.files(path=path,pattern ='*[Rr]2.fastq',
                         full.names = T,recursive = T)
  files_md5 <- list.files(path=path,pattern ='*.md5',
                          full.names = T,recursive = T)
  files_r1 <- setdiff(files_r1,files_md5)
  files_r2 <- setdiff(files_r2,files_md5)
  if (length(files_r1) == 0 | length(files_r2) == 0) stop('can not find r1/r2 files')
  if (length(files_r1) != length(files_r2)) stop('the number of r1/r2 files unequal')
  return(list(r1 = files_r1,r2 = files_r2))
}
# x shortread vector
# q quality threshold
# n number threshold
# r roi
.brtmatch.pass_quality_idx <- function(x,q,n,r = NULL){
  x.width <- Biostrings::width(x)[1]
  x %>%
    Biostrings::PhredQuality() %>%
    as('IntegerList') -> x.quality
  x.quality@unlistData %>%
    matrix(ncol = x.width,byrow = T) -> x.quality
  if (!is.null(r)) x.quality <- x.quality[,r]
  x.quality %>%
    apply(1,function(i) length(which(i<q)) <= n ) %>%
    which() -> idx
  return(idx)
}
# x shortread vector
# p pattern character
# m max mismatch
.brtmatch.matched <- function(x,p,m){
  p <- Biostrings::DNAString(p)
  matched <- Biostrings::vmatchPattern(pattern= p,
                                       subject = x@sread,
                                       fixed = "subject",
                                       max.mismatch = m)
  read.width <- max(x@sread@ranges@width)
  idx <- matched@ends %>%
    purrr::map(function(x) {
      if (is.null(x))
        FALSE
      else x[1] > 0 & x[1] < read.width & length(x) == 1
    }) %>%
    unlist()%>%
    which()
  x.read <- x@sread[idx][matched[idx]]
  x.qu <- x@quality[idx][matched[idx]]
  return(list(reads = x.read,quality = x.qu,idx = idx))
}
