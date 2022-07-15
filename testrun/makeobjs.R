# dpca1
rvBarcode_pattern='NNNACNNNGTNNNCGNNNTANNNCANNNTGNNN'
rvBarcode_roi=c(1:3,6:8,11:13,16:18,21:23,26:28,31:33)
path='D:/BRT/BRT/Data/xc1/DPCA1_BC'
seurat <- readRDS('D:/BRT/BRT/Data/xc1/DPCA1L1.rds')
cell_reference<-substr(colnames(seurat),1,16)
rr <- read.csv('D:/BRT/BRT/Data/brv/b5virus.csv')$rvBarcode
brtMatch_10x_bc(path=path,
                rvBarcode_pattern =rvBarcode_pattern,
                rvBarcode_roi = rvBarcode_roi,rvBarcode_cutoff = 3,
                cell_cutoff = 3,ncore=6,
                cell_reference = cell_reference,
                blockSize = 1e6,
                rvBarcode_reference = rr) -> dpca1.bcs

g.path<-'D:/BRT/BRT/Data/xc1/DPCA1_G'
g.pattern<-'catcagcagttgggaaagtcataaatctgggggcgagacacggctctaa'
g.patterns<-substr(g.pattern,1,25)
g <- brtMatch_10X_target(
  path = g.path,
  pattern = g.patterns,
  pattern_max_mismatch = 3,
  cell_reference = cell_reference,
  blockSize = 1e6
)
g <- rename(g,g = counts)
write.csv(dpca1.bcs,"testdata/dpca1_bcs.csv")
write.csv(g,"testdata/dpca1_g.csv")
# aca3

path='D:/BRT/BRT/Data/xc1/ACA3_BC'
seurat <- readRDS('D:/BRT/BRT/Data/xc1/ACA3L1.rds')
cell_reference<-substr(colnames(seurat),1,16)
brtMatch_10x_bc(path=path,
                rvBarcode_pattern =rvBarcode_pattern,
                rvBarcode_roi = rvBarcode_roi,rvBarcode_cutoff = 3,
                cell_cutoff = 3,ncore=1,
                cell_reference = cell_reference,
                blockSize = 1e6,
                rvBarcode_reference = rr) -> aca3.bcs
write.csv(aca3.bcs,"testdata/aca3_bcs.csv")
