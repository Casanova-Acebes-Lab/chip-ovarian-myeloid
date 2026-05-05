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


data <- readRDS(paste0(rsdir,"/objects/data.Clusterized1.rds"))


dataM <- subset(data, subset = c(Cluster == "Monocytes Ly6c2Hi" |
  Cluster == "Monocytes Ly6c2Lo" |
  Cluster == "Early TAMs" |
  Cluster == "Nlrp3|Vegfa Mac"
))

png(paste0(outdir,"/Subclustering.Monocytes/umap1.png"), width=1800, height=800)
DimPlot(dataM, reduction = "umap", split.by = "group") 
dev.off()


PCA.dataM <- PCA(dataM)
PCA.dataM[1]
PCA.dataM[2]
#18

dataM <- FindNeighbors(dataM, dims = 1:25)
dataM <- FindClusters(dataM, resolution = 0.3)

dataM <- SetIdent(dataM, value = "seurat_clusters")
png(paste0(outdir,"/Subclustering.Monocytes/umap.subC1.png"), width=1800, height=800)
DimPlot(dataM, reduction = "umap") 
dev.off()


png(paste0(outdir,"/Subclustering.Monocytes/Umap_Markers.IFN.png"), width=1400, height=3000)
FeaturePlot_scCustom(dataM, features= c("Ifit3","Ifit1", 
"Ifit2", "Ly6c2", "C1qa", "Ly6c2", "Ccr2", "Nlrp3", "percent.mt"), reduction = "umap",
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




png(paste0(outdir,"/Subclustering.Monocytes/Umap_Markers.png"), width=1400, height=3000)
FeaturePlot_scCustom(dataM, features= c("Ly6c2","C1qa", "C1qc", "Adgre4"), reduction = "umap",
split.by = "group", ncol=4)
dev.off()



# Número de células por cluster
table(Idents(dataM))

# Número total de counts por cluster
library(dplyr)
counts_per_cluster <- dataM@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(total_counts = sum(nCount_RNA), n_cells = n())
counts_per_cluster



genes_per_cluster <- dataM@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(mean_genes = mean(nFeature_RNA), median_genes = median(nFeature_RNA))
genes_per_cluster


library(dplyr)

cluster_summary <- dataM@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_cells = n(),
    total_counts = sum(nCount_RNA),
    mean_genes = mean(nFeature_RNA),
    median_genes = median(nFeature_RNA),
    mean_percent_mito = mean(percent.mt)
  )

cluster_summary


## Removing cluster 4 because it has low number of cells and low number of genes


# Filtrar las células que NO pertenecen al cluster 4
dataM_filtered <- subset(dataM, idents = c("0", "1", "2", "3", "5"))

# Normalización y análisis (opcional según necesites)
dataM_filtered <- NormalizeData(dataM_filtered)
dataM_filtered <- FindVariableFeatures(dataM_filtered)
dataM_filtered <- ScaleData(dataM_filtered,vars.to.regress = "percent.mt")
dataM_filtered <- RunPCA(dataM_filtered)
dataM_filtered <- FindNeighbors(dataM_filtered, dims = 1:20)
dataM_filtered <- FindClusters(dataM_filtered, resolution = 0.3)


png(paste0(outdir,"/Subclustering.Monocytes/umap2.png"), width=1800, height=800)
DimPlot(dataM_filtered, reduction = "umap", split.by = "group", pt.size = 0.7) 
dev.off()



png(paste0(outdir,"/Subclustering.Monocytes/Umap_Markers.filtered.png"), width=1400, height=3000)
FeaturePlot_scCustom(dataM_filtered, features= c("Ly6c2","C1qa", "C1qc", "percent.mt"), reduction = "umap",
split.by = "group", ncol=4)
dev.off()




markers <- FindAllMarkers(
  object = dataM_filtered,
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




## SingleR Cell identification


sce <- as.SingleCellExperiment(DietSeurat(dataM_filtered))

ref <- celldex::ImmGenData()

ref.main <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.main)

ref.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.fine)


dataM_filtered@meta.data$Cell_type.Image <- ref.main$pruned.labels
dataM_filtered@meta.data$Cell_type.Image.fine <- ref.fine$pruned.labels




dataM <- SetIdent(dataM_filtered, value = "Cell_type.Image.fine")
png(paste0(outdir,"/Subclustering.Monocytes/Umap_SingleR.Image.Fine.png"), width=2200, height=2400)
p0 <- DimPlot_scCustom(dataM_filtered, reduction = "umap", split.by="group", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
p1 <- DimPlot_scCustom(dataM_filtered, reduction = "umap", split.by="orig.ident", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
p2 <- DimPlot(dataM_filtered, reduction = "umap", split.by="orig.ident",group.by = "group", raster = FALSE)
p3 <- DimPlot(dataM_filtered, reduction = "umap", split.by="orig.ident",group.by = "tag", raster = FALSE)
plot_grid(p0,p1,p2,p3, ncol=1)
dev.off() 




dataM_filtered <- SetIdent(dataM_filtered, value = "Cell_type.Image.fine")
png(paste0(outdir,"/Subclustering.Monocytes/umap3.png"), width=1800, height=800)
DimPlot(dataM_filtered, reduction = "umap", split.by = "group", label=T, repel=T, label.size=4) 
dev.off()

dataM_filtered <- SetIdent(dataM_filtered, value = "seurat_clusters")
pdf(paste0(outdir,"/Subclustering.Monocytes/dotplot.pdf"), width=14, height=14)
DotPlot(dataM_filtered, features = c("Ly6c2","S100a8","S100a9","Adgre1","Mrc1","C1qa","C1qc","Ifit1","Csf1r","Cxcl9","Ccl5")) + RotatedAxis()
dev.off()

pdf(paste0(outdir,"/Subclustering.Monocytes/violin.pdf"), width=14, height=14)
VlnPlot(dataM_filtered, features = c("Ly6c2","Ifit2","Ifit3","Adgre4","Mrc1","C1qa","C1qc","Ifit1","Csf1r","Cxcl9","Ccl5", "Eno3", "Cd300e", "Cd9", "Ace",
"Fabp1", "Krt19", "Ccr2", "S100a8"), group.by="seurat_clusters")
dev.off()





dataM_filtered$Subclustering.Mon <- NA
dataM_filtered$Subclustering.Mon [dataM_filtered$seurat_clusters %in% c("0")] <- "Nlrp3|Vegfa Mac"
dataM_filtered$Subclustering.Mon [dataM_filtered$seurat_clusters %in% c("1")] <- "Ly6cHi Monocytes"
dataM_filtered$Subclustering.Mon [dataM_filtered$seurat_clusters %in% c("2")] <- "Ly6cHi Monocytes"
dataM_filtered$Subclustering.Mon [dataM_filtered$seurat_clusters %in% c("3")] <- "Early IFN TAMs"
dataM_filtered$Subclustering.Mon [dataM_filtered$seurat_clusters %in% c("4")] <- "Ly6cHi Monocytes"
dataM_filtered$Subclustering.Mon [dataM_filtered$seurat_clusters %in% c("5")] <- "Ly6cLo Monocytes"
dataM_filtered$Subclustering.Mon [dataM_filtered$seurat_clusters %in% c("6")] <- "Remove"
dataM_filtered$Subclustering.Mon [dataM_filtered$seurat_clusters %in% c("7")] <- "Remove"


# Cluster 7 looks a T cell contamination

dataM_filtered <- SetIdent(dataM_filtered, value = "Subclustering.Mon")
dataM_filtered <- subset(dataM_filtered, idents = "Remove", invert=TRUE)


dataM_filtered <- SetIdent(dataM_filtered, value = "Subclustering.Mon")
png(paste0(outdir,"/Subclustering.Monocytes/umap.clusterized.png"), width=1800, height=800)
DimPlot(dataM_filtered, reduction = "umap", pt.size = 0.) 
dev.off()



## Saving annotation

mon.annotation <- data.frame(
  barcode = colnames(dataM_filtered),
  Subclustering.Mon = dataM_filtered$Subclustering.Mon
)

write.csv(
  mon.annotation,
  file = paste0(outdir,"/Subclustering.Monocytes/subclustering_monocytes_annotations.csv"),
  row.names = FALSE
)
