library(Seurat)
library(dplyr)
library(cowplot)
library(scCustomize)
library(ggplot2)
library(RPresto)
library(celldex)
library(SingleCellExperiment)
library(SingleR)




# Reading th config file

config <- read.table("config.tsv",
                     sep = "\t",
                     header = T)

keys <- config$key
values <- config$value
names(values) <- keys
###############
rsdir <- values["rsdir"]
datadir <- values["datadir"]
outdir <- values["outdir"]


source("R/Functions.R")



###


data <- readRDS(paste0(rsdir,"objects/data.Clusterized1.rds"))


dataM <- subset(data, subset = Cluster == "Mki67 Mac"
)


png(paste0(outdir,"/Subclustering.Mi67Mac/umap1.png"), width=1800, height=800)
DimPlot(dataM, reduction = "umap", split.by = "group") 
dev.off()


PCA.dataM <- PCA(dataM)
PCA.dataM[1]
PCA.dataM[2]
#18

dataM <- FindNeighbors(dataM, dims = 1:25)
dataM <- FindClusters(dataM, resolution = 0.2)

dataM <- SetIdent(dataM, value = "seurat_clusters")
png(paste0(outdir,"/Subclustering.Mi67Mac/umap.subC1.png"), width=1800, height=800)
DimPlot(dataM, reduction = "umap") 
dev.off()




markers <- FindAllMarkers(
  object = dataM,
  only.pos = TRUE,      # Solo genes sobreexpresados positivos
  min.pct = 0.25,       # Mínimo porcentaje de células que expresan el gen
  logfc.threshold = 0.25 # Umbral de log2 fold change
)


top30_markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50) %>%   # Selecciona los top 30 por cluster
  ungroup()


as.data.frame(top30_markers)

dataM@meta.data$Mki67C <- NA
dataM@meta.data$Mki67C[dataM$seurat_clusters == 0]<- "Mki67 Mac 1"
dataM@meta.data$Mki67C[dataM$seurat_clusters == 1] <- "Mki67 Mac 2"
dataM@meta.data$Mki67C[dataM$seurat_clusters == 2]<- "Mki67 Mac 2"




dataM <- SetIdent(dataM, value = "Mki67C")
png(paste0(outdir,"/Subclustering.Mi67Mac/umap.subC2.png"), width=1800, height=800)
DimPlot(dataM, reduction = "umap") 
dev.off()







dataM <- SetIdent(dataM, value = "Mki67C")
markers <- FindAllMarkers(
  object = dataM,
  only.pos = TRUE,      # Solo genes sobreexpresados positivos
  min.pct = 0.25,       # Mínimo porcentaje de células que expresan el gen
  logfc.threshold = 0.25 # Umbral de log2 fold change
)


top30_markers <- markers %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 50) %>%   # Selecciona los top 30 por cluster
  ungroup()

  as.data.frame(top30_markers)



genes_to_plot <- c("Mki67", "Cdkn1a")
png(paste0(outdir, "/Subclustering.Mi67Mac/Violin_Markers.png"), width = 1400, height = 1000)
VlnPlot(
  object = dataM,
  features = genes_to_plot,
  group.by = "Mki67C"
)
dev.off()


genes_to_plot <- c("Mki67", "Cdkn1a")
png(paste0(outdir, "/Subclustering.Mi67Mac/Feature_Markers.png"), width = 1400, height = 1000)
FeaturePlot_scCustom(dataM, features= c("Ctsk","Mmp9", "S100a4", "S100a10"), reduction = "umap",
split.by = "group", ncol=4)
dev.off()




dataM@meta.data$Mki67.Clusters <- NA
dataM@meta.data$Mki67.Clusters[dataM@meta.data$Mki67C == "Mki67 Mac 1"]<- "Mki67 IFN Mac"
dataM@meta.data$Mki67.Clusters[dataM@meta.data$Mki67C == "Mki67 Mac 2"] <- "Mki67|Cstk|Mmp9|S100a4 Mac"


umap <- Embeddings(dataM, "umap")

set.seed(123)

km <- kmeans(
  umap[,1:2],
  centers = 2
)


data$UMAP_kmeans2 <- as.factor(km$cluster)


png(paste0(outdir, "/Subclustering.Mi67Mac/kmeans.png"), width = 1400, height = 1000)
DimPlot(data, group.by = "UMAP_kmeans2")
dev.off()

## Saving annotation

mon.annotation <- data.frame(
  barcode = colnames(dataM),
  Subclustering.Mki67 = dataM@meta.data$Mki67.Clusters
)

write.csv(
  mon.annotation,
  file = paste0(outdir,"/Subclustering.Mi67Mac/subclustering_Mki67.csv"),
  row.names = FALSE
)

