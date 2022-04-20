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
#  1 Seurat Alignment
#
# ######################################################################
#======1.1 Setup the Seurat objects =========================
# read sample summary list
df_samples <- readxl::read_excel("doc/20200730_scRNAseq_info.xlsx")
df_samples = as.data.frame(df_samples)
colnames(df_samples) %<>% tolower()
#======1.2 load  Seurat =========================
object1 = readRDS(file = "data/SSC_20210826.rds")
meta.data = object@meta.data
for(i in 1:length(df_samples$sample.id)){
    cells <- meta.data$orig.ident %in% df_samples$sample[i]
    print(df_samples$sample.id[i])
    print(table(cells))
    meta.data[cells,"patient"] = df_samples$patient[i]
    meta.data[cells,"Mean.Reads.per.Cell"] = df_samples$`mean reads per cell`[i]
    meta.data[cells,"Number.of.Reads"] = df_samples$`number of cells`[i]
    meta.data[cells,"Sequencing.Saturation"] = df_samples$`median genes per cell`[i]
}
meta.data$orig.ident %<>% factor(levels = df_samples$sample)
table(rownames(object@meta.data) == rownames(meta.data))
table(colnames(object) == rownames(meta.data))
object@meta.data = meta.data

#======1.7 UMAP from raw pca =========================
object %<>% SCTransform(method = "glmGamPoi", vars.to.regress = "percent.mt", verbose = TRUE)

object <- FindVariableFeatures(object = object, selection.method = "vst",
                               num.bin = 20, nfeatures = 2000,
                               mean.cutoff = c(0.1, 8), dispersion.cutoff = c(1, Inf))
object %<>% ScaleData(verbose = FALSE)
object %<>% RunPCA(verbose = T,npcs = 100)

jpeg(paste0(path,"S1_ElbowPlot_SCT.jpeg"), units="in", width=10, height=7,res=600)
ElbowPlot(object, ndims = 100)
dev.off()

saveRDS(object, file = "data/SSC_20210826.rds")

#======1.8 UMAP from harmony =========================
DefaultAssay(object) = "SCT"


npcs = 25
object %<>% RunUMAP(reduction = "pca", dims = 1:npcs)
system.time(object %<>% RunTSNE(reduction = "pca", dims = 1:npcs))
object %<>% FindNeighbors(reduction = "umap",dims = 1:2)
object %<>% FindClusters(resolution = 0.8)

resolutions = c(seq(0.01,0.09, by = 0.01),seq(0.1,0.9, by = 0.1))
for(i in 1:length(resolutions)){
    object %<>% FindClusters(resolution = resolutions[i], algorithm = 1)
    Progress(i,length(resolutions))
}

resolutions = seq(0.1,0.9, by = 0.1)
for(i in 1:length(resolutions)){
    object[[paste0("SCT_snn_res.",resolutions[i])]] = NULL
}
saveRDS(object, file = "data/SSC_20210826.rds")


#=======1.9 save SCT only =======================================
format(object.size(object),unit = "GB")

format(object.size(object@assays$RNA),unit = "GB")
object[['RNA']] <- NULL
object[["SCT"]]@scale.data = matrix(0,0,0)
format(object.size(object),unit = "GB")
saveRDS(object, file = "data/SSC_SCT_20210826.rds")


