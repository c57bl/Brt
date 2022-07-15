# ============================  starter  ========================================================
rvBarcode_pattern='NNNACNNNGTNNNCGNNNTANNNCANNNTGNNN'
rvBarcode_roi=c(1:3,6:8,11:13,16:18,21:23,26:28,31:33)
path='F:/BRT/brt8/bc'
brt8seurat <- readRDS('D:/BRT/BRT/Data/brt8/objs/brt8seurat.rds')
brt8seurat_neurons <- readRDS('D:/BRT/BRT/Data/brt8/objs/brt8seurat_neurons.rds')
cell_reference<-substr(colnames(brt8seurat),1,16)
rr <- read.csv('D:/BRT/BRT/Data/brv/b5virus.csv')$rvBarcode
brtMatch_10x_bc(path=path,
                rvBarcode_pattern =rvBarcode_pattern,
                rvBarcode_roi = rvBarcode_roi,rvBarcode_cutoff = 3,
                cell_cutoff = 3,ncore=6,
                cell_reference = cell_reference,
                blockSize = 1e6,
                rvBarcode_reference = rr) -> a
brt8starter <- brtStarter_from_tibble(a,seuratobj = brt8seurat,
                               sample = "brt8",source = "mpfc")
# ======================== re-cluster ===================================================================
(brt8starter@coldata %>%
    subset(ident != "Neurons") %>%
    select(rvcounts.raw))[[1]] -> rv.background
rv.threshold <- median(rv.background) + 3 * IQR(rv.background)
brtRecluster_knn(seurat_obj = brt8seurat_neurons,
                 drop = subset(brt8starterL2@coldata,rvcounts.raw > rv.threshold)$cell,
                 p_train = 0.7,
                 resolution = 0.7,
                 findmarker = T) -> brt8seurat_neurons_2
brt8seurat_neurons@active.ident <- brt8seurat_neurons_2$clusters
saveRDS(list(cluster = brt8seurat_neurons_2$clusters,
             markers = brt8seurat_neurons_2$seurat.sub.markers),
        "../data/brt8/cacher/brt8_reclustering.rds")
# =========================== re-cluster remove degs ====================================
# brt8seurat_neurons@meta.data <- brt8seurat_neurons@meta.data %>%
#   mutate(cell = rownames(brt8seurat_neurons@meta.data))
# cell.p <- subset(brt8starterL2@coldata, rvcounts.raw > (5*rv.threshold))$cell
# cell.n <- subset(brt8starterL2@coldata, rvcounts.raw < 3)$cell
# rbind(data.frame(cell = paste0(cell.p, "-1"), rv = "high"),
#       data.frame(cell = paste0(cell.n, "-1"), rv = "low")) -> cell.rv
# brt8seurat_neurons@meta.data %>%
#   left_join(cell.rv,
#             by = "cell") -> brt8seurat_neurons@meta.data
# FindMarkers(brt8seurat_neurons,
#             ident.1 = "high",
#             ident.2 = "low",
#             group.by = "rv") -> rvdegs
# brt8seurat_neurons@meta.data <- select(brt8seurat_neurons@meta.data,-rv)
# rownames(brt8seurat_neurons@meta.data) <- brt8seurat_neurons@meta.data$cell
# runSeurat_dropgene(brt8seurat_neurons,
#                    0.6,
#                    F,
#                    rownames(subset(rvdegs,p_val_adj < 0.01))) -> brt8seurat_neurons_sub
# ========================== helpers ===================================================================
tva.path<-'F:/BRT/brt8/tva'
g.path<-'F:/BRT/brt8/g'
tva.pattern<-'agctgctggtgcctgacaagagccaggcagacttgttctccggtagtgga'
tva.patterns<-substr(tva.pattern,1,25)
g.pattern<-'tcatcagcagctgggagagccacaagagcggcggcgagaccagactgtga'
g.patterns<-substr(g.pattern,1,25)
brt8.tva<-brtMatch_10X_target(path=tva.path,
                              pattern =tva.patterns,
                              pattern_max_mismatch = 3,
                              cell_reference = cell_reference,
                              blockSize = 1e6)
brt8.g<-brtMatch_10X_target(path=g.path,
                            pattern =g.patterns,
                            pattern_max_mismatch = 3,
                            cell_reference = cell_reference,
                            blockSize = 1e6,
                            total = 8)
brt8.tva <- rename(brt8.tva,tva = counts)
brt8.g <- rename(brt8.g, g = counts)
write.csv(brt8.tva,"brt8_tva.csv")
write.csv(brt8.g,"brt8_g.csv")
brt8starter@coldata <- left_update(brt8starter@coldata,brt8.tva,by = "cell",treat.na = 0)
brt8starter@coldata <- left_update(brt8starter@coldata,brt8.g,by = "cell",treat.na = 0)
saveRDS(brt8starter,"../Data/brt8/objs/brt8starter.rds")
brt8starterL2 <- brtStarter_from_tibble(a,seuratobj = brt8seurat_neurons,
                                        sample = "brt8",source = "mpfc")
brt8starterL2@coldata <- left_update(brt8starterL2@coldata,brt8.tva,by = "cell",treat.na = 0)
brt8starterL2@coldata <- left_update(brt8starterL2@coldata,brt8.g,by = "cell",treat.na = 0)
saveRDS(brt8starterL2,"../Data/brt8/objs/brt8starterL2.rds")
# ========================= inputome ==============================================================
rr <- read.csv('../Data/brv/b5virus.csv')$rvBarcode
index <- read.csv('../Data/indexs.csv')
RvPatterns<-read.csv('../Data/brv/RvBarcodes.csv')
RvPatterns<-RvPatterns[c(2,4),]
rvBarcode_roi=c(1:3,6:8,11:13,16:18,21:23,26:28,31:33)
umi_pattern='NNNNNNNNNNNNNNNNNNNNNNNCTCGACTGAAAAGCT'
inputs_table<-brtMatch_input(path ='F:/BRT/brt8/inputome',
                             RvPatterns = RvPatterns,
                             umiPattern=umi_pattern,
                             indexs=index,
                             rvBarcode_roi = rvBarcode_roi,
                             umi_roi = 1:15,
                             index_roi = 16:23,
                             Batch="brt8",ncore = 1,
                             blockSize = 1e6,
                             rvBarcode_reference = rr)
saveRDS(inputs_table,"../Data/brt8/cacher/inputs_table.rds")
brt8inputome <- brtInputome_from_tibble(inputs_table,sample = "brt8",
                                        source = "mpfc",bc = "bcV1",
                                        spk = "spkV2")
brt8coldata <-read.csv("../Data/brt8/cacher/brtcoldata.csv")
brt8inputome@coldata <- left_join(brt8inputome@coldata,brt8coldata,by = "id")
saveRDS(brt8inputome,"../Data/brt8/objs/brt8inputome.rds")

# brtconfig <- list(starter = list(rvBarcode_pattern='NNNACNNNGTNNNCGNNNTANNNCANNNTGNNN',
#                                  rvBarcode_roi=c(1:3,6:8,11:13,16:18,21:23,26:28,31:33),
#                                  tva.pattern = "agctgctggtgcctgacaagagcca",
#                                  g.pattern = "tcatcagcagctgggagagccacaa",
#                                  ),
#                   inputs = list(rv_in_use = c(2,4),
#                                 rvBarcode_roi=c(1:3,6:8,11:13,16:18,21:23,26:28,31:33),
#                                 index_in_use = 1:91))


