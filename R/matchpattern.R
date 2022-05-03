# ============================ match 10x bc ====================================
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
                            total = NULL){
  paths <- .brtmatch.parse_path(path)
  loop <- 0
  results <- list()
  message('start match')
  for (i in 1:length(paths$r1)){
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
        results.sub <- .brtmatch_10x_bc.count_character(cbind(curr_cellbarcode,
                                                              curr_bc,
                                                              curr_umi),
                                                        cr = cell_reference,
                                                        rr = rvBarcode_reference)
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
  }
  results <- do.call(rbind,results)
  results <- .brtmatch_10x_bc.count_n(results)
  close(con1)
  close(con2)
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

# =========================== match 10x features ===============================
#' brtMatch_10x_target
#' @description return a dataframe record the umi counts of a given pattern
#' @export
brtMatch_10X_target <- function(path,pattern,quality_cutoff=25,
                                cellBelow=2,patternBelow=2,blockSize=2e6,
                                pattern_max_mismatch=0,cell_reference = NULL,
                                total = NULL){
  paths <- .brtmatch.parse_path(path)
  loop <- 0
  results <- list()
  message(Sys.time(),'\n',
          ' start match...')
  for (i in 1:length(paths$r1)){
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

        results.sub <- .brtmatch_10x_target.count_character(cbind(curr_cellbarcode,
                                                            curr_umi),
                                                        cr = cell_reference)
        results[[loop]] <- results.sub
        message(Sys.time(),"\nloops: ",
                loop," ",length(curr_cellbarcode),' matched')
        return(results)
      }
      tryCatch(results <- main(),
               error = function(e){
                 message(Sys.time(),"\nloops: "
                         ,loop," error, ","probably due to limited reads")
                 return(NULL)
               })
    }
    results <- do.call(rbind,results)
    results <- results %>%
      group_by(cell) %>%
      summarise(counts = sum(counts))
    close(con1)
    close(con2)
    return(results)
  }
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
  idx <- matched@ends %>%
         purrr::map(~length(.x)==1) %>%
         unlist()%>%
         which()
  x.read <- x@sread[idx][matched[idx]]
  x.qu <- x@quality[idx][matched[idx]]
  return(list(reads = x.read,quality = x.qu,idx = idx))
}
