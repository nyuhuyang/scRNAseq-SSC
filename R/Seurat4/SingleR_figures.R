# conda activate r4.0.3
library(Seurat)
library(magrittr)
library(kableExtra)
library(dplyr)
library(tidyr)
library(ggpubr)
library(S4Vectors)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")

path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)
#====== 3.2 SingleR specifications ==========================================

##############################
# create singleR data frame
###############################
pred = readRDS("output/SSC_SCT_20210826_singleR_pred.rds")
object = readRDS("data/SSC_SCT_20210826.rds")

singlerDF = data.frame("cell.types" = pred$pruned.labels,
                       row.names = rownames(pred))
table(is.na(singlerDF$cell.types))
singlerDF$cell.types[is.na(singlerDF$cell.types)]= "unknown"


# combine cell types


object@meta.data %<>% cbind(singlerDF)


lapply(c(TSNEPlot.1,UMAPPlot.1), function(fun)
    fun(object = object, label = T, label.repel = T,group.by = "cell.types",
        no.legend = T,
        pt.size = 0.1,label.size = 3,
        do.print = T,do.return = F,
        title ="labeling by blue_encode and TCGA-BLCA RNA-seq"))

saveRDS(object, file = "data/SSC_SCT_20210826.rds")