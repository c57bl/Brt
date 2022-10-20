# rvbc
rr <- read.csv('data/b5virus.csv')$rvBarcode
RvPatterns <- read.csv('data/RvBarcodes.csv')
rvBarcode_roi=c(1:3,6:8,11:13,16:18,21:23,26:28,31:33)
# starter tva
tva.pattern <- 'agctgctggtgcctgacaagagccaggcagacttgttctccggtagtgga'
tva.patterns <- substr(tva.pattern, 1, 25)
# starter g
g.pattern <- 'tcatcagcagctgggagagccacaagagcggcggcgagaccagactgtga'
g.patterns <- substr(g.pattern, 1, 25)
# inputome
index <- read.csv('data/indexs.csv')
index.roi = 16:23
umi_pattern = 'NNNNNNNNNNNNNNNNNNNNNNNCTCGACTGAAAAGCT'
umi.roi = 1:15
regions <- read.csv("data/regionsV2.csv")

neuron.subclass <- c("Car3","L2/3 IT ","L4/5 IT ","L5 IT ","L5 PT ",
                        "L5/6 NP ","L6 CT ","L6 IT ","L6b ","Lamp5",
                        "Meis2","Pvalb","Sncg","Sst","Sst Chodl","Vip")
# build parmas
brtParams <- list(rr = rr,
               rvpatterns = RvPatterns,
               rvroi = rvBarcode_roi,
               tva.pattern = tva.patterns,
               g.pattern = g.patterns,
               index = index,
               indexroi = index.roi,
               umipattern = umi_pattern,
               umiroi = umi.roi,
               neuron.subclass = neuron.subclass,
               regions = regions
               )
save(brtParams,file = "data/brtParams.rdata")
