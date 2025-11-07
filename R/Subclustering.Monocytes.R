library(Seurat)
library(dplyr)
library(cowplot)
library(scCustomize)
library(ggplot2)
library(RPresto)





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


data <- readRDS(paste0(rsdir,"objects/data.Clusterized.rds"))


dataM <- subset(data, subset = c(Cluster == "Ly6c|Ms4ac Monocytes"|
Cluster == "Slamf Monocytes"))


png(paste0(outdir,"/Subclustering.Monocytes/umap1.png"), width=1800, height=800)
DimPlot(dataM, reduction = "umap", split.by = "group") 
dev.off()



PCA.dataM <- PCA(dataM)
PCA.dataM[1]
PCA.dataM[2]

dataM <- FindNeighbors(dataM, dims = 1:14)
dataM <- FindClusters(dataM, resolution = 0.2)

dataM <- SetIdent(dataM, value = "seurat_clusters")
png(paste0(outdir,"/Subclustering.Monocytes/umap.subC1.png"), width=1800, height=800)
DimPlot(dataM, reduction = "umap") 
dev.off()


png(paste0(outdir,"/Subclustering.Monocytes/Umap_Markers.IFN.png"), width=1400, height=3000)
FeaturePlot_scCustom(dataM, features= c("Ifit3","Ifit1", 
"Ifit2", "Ly6c2", "C1qa","Adgre4", "Ace", "Pglyrp1"), reduction = "umap",
split.by = "group", ncol=4)
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




pdf(paste0(outdir,"/Subclustering.Monocytes/Umap_Markers.pdf"), width=14, height=30)
FeaturePlot_scCustom(data, features= c("Ly6c2","C1qa", "C1qc", "Adgre4"), reduction = "umap",
split.by = "group", ncol=4)
dev.off()




## SingleR Cell identification


sce <- as.SingleCellExperiment(DietSeurat(dataM))

ref <- celldex::ImmGenData()

ref.main <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.main)

ref.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.fine)


dataM@meta.data$Cell_type.Image <- ref.main$pruned.labels
dataM@meta.data$Cell_type.Image.fine <- ref.fine$pruned.labels




dataM <- SetIdent(dataM, value = "Cell_type.Image.fine")
png(paste0(outdir,"/Subclustering.Monocytes/Umap_SingleR.Image.Fine.png"), width=2200, height=2400)
p0 <- DimPlot_scCustom(dataM, reduction = "umap", split.by="group", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
p1 <- DimPlot_scCustom(dataM, reduction = "umap", split.by="orig.ident", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
p2 <- DimPlot(dataM, reduction = "umap", split.by="orig.ident",group.by = "group", raster = FALSE)
p3 <- DimPlot(dataM, reduction = "umap", split.by="orig.ident",group.by = "tag", raster = FALSE)
plot_grid(p0,p1,p2,p3, ncol=1)
dev.off() 


pdf(paste0(outdir,"/Subclustering.Monocytes/dotplot.pdf"), width=14, height=14)
DotPlot(dataM, features = c("Ly6c2","S100a8","S100a9","Adgre1","Mrc1","C1qa","C1qc","Ifit1","Csfr1","Cxcl9","Ccl5")) + RotatedAxis()
dev.off()

pdf(paste0(outdir,"/Subclustering.Monocytes/violin.pdf"), width=14, height=14)
VlnPlot(dataM, features = c("Ly6c2","Ifit2","Ifit3","Adgre4","Mrc1","C1qa","C1qc","Ifit1","Csf1r","Cxcl9","Ccl5"), group.by="seurat_clusters")
dev.off()


dataM$Subclustering.Mon <- NA
dataM$Subclustering.Mon [dataM$seurat_clusters %in% c("0")] <- "Early IFN-TAMs"
dataM$Subclustering.Mon [dataM$seurat_clusters %in% c("1")] <- "Ly6cHi Monocytes"
dataM$Subclustering.Mon [dataM$seurat_clusters %in% c("2")] <- "Ly6cLo Monocytes"
dataM$Subclustering.Mon [dataM$seurat_clusters %in% c("3")] <- "Ly6cHi Monocytes"
dataM$Subclustering.Mon [dataM$seurat_clusters %in% c("4")] <- "Early IFN-TAMs"


dataM <- SetIdent(dataM, value = "Subclustering.Mon")
png(paste0(outdir,"/Subclustering.Monocytes/umap.clusterized.png"), width=1800, height=800)
DimPlot(dataM, reduction = "umap") 
dev.off()



## Saving annotation

mon.annotation <- data.frame(
  barcode = colnames(dataM),
  Subclustering.Mon = dataM$Subclustering.Mon
)

write.csv(
  mon.annotation,
  file = paste0(outdir,"/Subclustering.Monocytes/subclustering_monocytes_annotations.csv"),
  row.names = FALSE
)
