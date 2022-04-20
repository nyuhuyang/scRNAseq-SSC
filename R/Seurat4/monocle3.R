# conda activate r4.1.1
library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(Matrix)
library(ggplot2)
library(patchwork)
library(data.table)
source("https://raw.githubusercontent.com/nyuhuyang/SeuratExtra/master/R/Seurat4_functions.R")
path <- paste0("output/",gsub("-","",Sys.Date()),"/")
if(!dir.exists(path)) dir.create(path, recursive = T)

object = readRDS("data/SSC_SCT_20211129.rds")

UMAPPlot(object,group.by = "SCT_snn_res.0.07")
object$cell.labels = plyr::mapvalues(object$SCT_snn_res.0.07,
                                     from = 0:4,
                                     to =c("late SSCs",
                                           "Spg3",
                                           "Spg1",
                                           "early SSCs",
                                           "Spg2"))
UMAPPlot(object,group.by = "cell.labels",label = T)

#============== re-run UMAP ==============================
# Building trajectories with Monocle 3

cds <- as.cell_data_set(object)
cds <- cluster_cells(cds = cds, reduction_method = "UMAP",k = 13,verbose = TRUE)
cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)


get_earliest_principal_node <- function(cds, cell.type="early SSCs"){
    cell_ids <- which(colData(cds)[, "cell.labels"] == cell.type)
    
    closest_vertex <-
        cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <-
        igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                                  (which.max(table(closest_vertex[cell_ids,]))))]
    
    root_pr_nodes
}
#cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
cds <- order_cells(cds)
#=============== prepare figures ===================

for(label in c("subtype_raw","subtype","subtype_label")){
    g <- plot_cells(
        cds = cds,
        label_cell_groups = switch(label,
                                   "subtype_raw" = TRUE,
                                   "subtype" = TRUE,
                                   "subtype_label" = TRUE),
        color_cells_by = "cell.labels",
        show_trajectory_graph = switch(label,
                                       "subtype_raw" = FALSE,
                                       TRUE),
        group_label_size = switch(label,
                                  "subtype_raw" = 5,
                                  "subtype_label" = 5,
                                  0))
    jpeg(paste0(path,"SSC_",label,".jpeg"), units="in", width=7, height=7,res=600)
    print(g)
    dev.off()
}

for(label in c("pseudotime","pseudotime_label")){
    g1 <-     plot_cells(
                cds = cds,
                color_cells_by = "pseudotime",
                show_trajectory_graph = switch(label,
                                               "pseudotime" = FALSE,
                                               "pseudotime_label" = TRUE
                                               )
                )
    jpeg(paste0(path,"SSC_",label,".jpeg"), units="in", width=7, height=7,res=600)
    print(g1)
    dev.off()
}
    


# save coordinates
object$pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
object$monocle_clusters <- cds@clusters[["UMAP"]]$clusters
saveRDS(object, file = "data/SSC_SCT_20211129.rds")


for(k in c(7,8,9,10,11,12,13)){
    print(k)
    cds <- as.cell_data_set(object)
    cds <- cluster_cells(cds = cds, k=k,reduction_method = "UMAP")
    cds <- learn_graph(cds, use_partition = TRUE,close_loop = TRUE)
    g <- plot_cells(
        cds = cds,
        label_cell_groups = TRUE,
        color_cells_by = "cluster",
        group_cells_by = c("cluster", "partition"),
        show_trajectory_graph = TRUE,
        group_label_size = 5)
    jpeg(paste0(path,"SSC_",k,".jpeg"), units="in", width=7, height=7,res=600)
    print(g)
    dev.off()
}
