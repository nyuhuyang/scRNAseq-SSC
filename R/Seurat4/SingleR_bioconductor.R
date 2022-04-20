#====== 3.1 Create Singler Object  ==========================================
# conda activate r4.0.3 linux
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(magrittr)
library(TCGAbiolinks)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/SingleR_functions.R")

# ====== load single cell =============

object = readRDS("data/SSC_SCT_20210826.rds")
sce <- SingleCellExperiment(list(logcounts=object[["SCT"]]@data),
                                colData=DataFrame(object@meta.data))
rm(object);GC()

# ====== load reference =============
(load("data/SSCs_20180926.Rda"))
SSCs = UpdateSeuratObject(SSCs)
SSCs$label.fine = SSCs$Cell.Types

SSC_sce <- SingleCellExperiment(list(logcounts=SSCs[["RNA"]]@data),
                            colData=DataFrame(SSCs@meta.data))
rm(SSCs);GC()

common <- Reduce(intersect, list(rownames(sce),
                                 rownames(SSC_sce)
))
length(common)

table(SSC_sce$label.fine)
system.time(trained <- trainSingleR(ref = SSC_sce[common,],
                                    labels=SSC_sce$label.fine))
system.time(pred <- classifySingleR(sce[common,], trained))
# elapsed 4872.846 sec
saveRDS(object = pred, file = "output/SSC_SCT_20210826_singleR_pred.rds")