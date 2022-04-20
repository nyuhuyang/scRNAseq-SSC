########################################################################
#
#  0 setup environment, install libraries if necessary, load libraries
#
# ######################################################################
# conda activate r4.0.3
#devtools::install_github("immunogenomics/harmony", ref= "ee0877a",force = T)
invisible(lapply(c("Seurat","dplyr","ggplot2","cowplot","pbapply","harmony","sctransform"), function(x) {
    suppressPackageStartupMessages(library(x,character.only = T))
}))
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

########################################################################
#
#  1 Seurat
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20200730_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
#======1.2 load  Seurat =========================
object = readRDS(file = "data/SSC_SCT_20210826.rds")
object %<>% subset(subset = cell.types %in% "Spermatogonia")
object@meta.data %<>% cbind(object[["umap"]]@cell.embeddings)
keep <- object$UMAP_1 > -5.5 & object$UMAP_1 < 10 & object$UMAP_2 < 6.4
object$discard = !keep
object %<>% subset(subset = discard == FALSE)
UMAPPlot(object, group.by = "SCT_snn_res.0.08",label = T)
object$Cell_subtype = plyr::mapvalues(object$SCT_snn_res.0.08, 
                                      from = 0:4,
                                      to = c("Quiescent SSCs","Active SSCs",
                                             "Early progenitors","Late progenitors",
                                             "Transitional"))
object$Cell_subtype %<>% factor(levels = c("Quiescent SSCs","Active SSCs",
                                           "Transitional",
                                           "Early progenitors","Late progenitors"
))
s.genes <- cc.genes$s.genes %>% gsub("MLF1IP","CENPU",.)
g2m.genes <- cc.genes$g2m.genes %>% plyr::mapvalues(from = c("FAM64A", "HN1"),
                                                    to = c("PIMREG","JPT1"))

s.genes <- Hmisc::capitalize(tolower(s.genes))
g2m.genes <- Hmisc::capitalize(tolower(g2m.genes))
object %<>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

object %<>% ScaleData(vars.to.regress = c("S.Score", "G2M.Score"))
object %<>% RunPCA(features = VariableFeatures(object), nfeatures.print = 10)

#======1.8 UMAP from harmony =========================
object = readRDS("data/SSC_SCT_20211129.rds")


npcs = 25
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
UMAPPlot.1(object, do.print = T,group.by = "Phase")

resolutions = c(seq(0.01,0.09, by = 0.01))
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}
object[["SCT"]]@scale.data = matrix(0,0,0)
saveRDS(object, file = "data/SSC_SCT_20211129.rds")