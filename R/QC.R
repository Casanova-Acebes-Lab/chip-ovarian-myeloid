
library(Seurat)
library(dplyr)
library(cowplot)
library(scCustomize)
library(harmony)
library(ggplot2)
library(RPresto)
library(celldex)
library(SingleR)
library(SingleCellExperiment)




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



# Reading data 

smpR065_DsRedN_KO2 <- read.data("SoupX/smR065/Seurat_smpR065_DsRedN-KO2.2.rds",
"DsRedN-KO2",
"KO",
"DsRedN", 
"smpR065",
"C"
)

smpR065_DsRedN_KO3 <- read.data("SoupX/smR065/Seurat_smpR065_DsRedN-KO3.2.rds",
"DsRedN-KO3",
"KO",
"DsRedN", 
"smpR065",
"C"
)


smpR065_DsRedP_KO2 <- read.data("SoupX/smR065/Seurat_smpR065_DsRedP-KO2.2.rds",
"DsRedP-KO2",
"KO",
"DsRedP", 
"smpR065",
"C"
)

smpR065_DsRedP_KO3 <- read.data("SoupX/smR065/Seurat_smpR065_DsRedP-KO3.2.rds",
"DsRedP-KO3",
"KO",
"DsRedP", 
"smpR065",
"C"
)


smR039a_WT1_DsRedn <- read.data("SoupX/SmR039a.9/Seurat_WT1-DsRedn.2.rds",
"WT1-DsRedN",
"WT",
"DsRedN", 
"smR039a",
"A"
)

smR039a_WT1_DsRedp <- read.data("SoupX/SmR039a.9/Seurat_WT1-DsRedp.2.rds",
"WT1-DsRedP",
"WT",
"DsRedP", 
"smR039a",
"A"
)


smR039b_WT2_DsRedn <- read.data("SoupX/SmR039b.9/Seurat_WT2-DsRedn.2.rds",
"WT2_DsRedN",
"WT",
"DsRedN", 
"smR039b",
"B"
)

smR039b_WT2_DsRedp <- read.data("SoupX/SmR039b.9/Seurat_WT2-DsRedp.2.rds",
"WT2_DsRedP",
"WT",
"DsRedP", 
"smR039b",
"B"
)


# Merge

object_list <- list(smpR065_DsRedN_KO2, smpR065_DsRedN_KO3, smpR065_DsRedP_KO2, 
smpR065_DsRedP_KO3, smR039a_WT1_DsRedn, smR039a_WT1_DsRedp,
smR039b_WT2_DsRedn, smR039b_WT2_DsRedp)
data <- Merge_Seurat_List(list_seurat = object_list)



data <- JoinLayers(data)


# Removing doublets

data <- subset(data, subset = DF.class == "Doublet", invert=TRUE)

data$batch <- factor(data$batch, levels = c("A", "B", "C"))

# QC 

data$percent.mt <- PercentageFeatureSet(data, pattern = "^mt-")

data <- SetIdent(data, value = "batch")

png(paste0(outdir,"/QC/QC.Before_filtering.Complete.png"), width=1200, height=600)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0,raster=FALSE)
dev.off()


plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)


png(paste0(outdir,"/QC/ScatterQC.Before_filtering.Complete.png"), width=1200, height=800)
plot_grid(plot1,plot2, ncol=2)
dev.off()


## By sample 

data <- SetIdent(data, value = "tag")

png(paste0(outdir,"/QC/QC.Before_filtering.Complete.tag.png"), width=1200, height=600)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, raster = FALSE)
dev.off()


plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)


png(paste0(outdir,"/QC/ScatterQC.Before_filtering.Complete.tag.png"), width=1200, height=800)
plot_grid(plot1,plot2, ncol=2)
dev.off()


data <- subset(data, subset = nFeature_RNA > 300 & nFeature_RNA < 8500 & percent.mt < 10)


data <- SetIdent(data, value = "batch")

png(paste0(outdir,"/QC/QC.After_filtering.Complete.png"), width=1200, height=800)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0 , raster = FALSE)
dev.off()


plot3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)
plot4 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)


png(paste0(outdir,"/QC/ScatterQC.After_filtering.Complete.png"), width=1200, height=800)
plot_grid(plot1,plot2,plot3,plot4, ncol=2)
dev.off()


# Analysis

##### Normalize

data<- NormalizeData(data)
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)



### Scale

data <- ScaleData(data, features = VariableFeatures(data),vars.to.regress = c("nCount_RNA"))


# Linear Reduction

data <- RunPCA(data, features = VariableFeatures(object = data))


png(paste0(outdir,"/QC/ElbowPlot.Complete.png"), width=800, height=600)
ElbowPlot(data)
dev.off()


##### PCA

PCA.data <- PCA(data)
PCA.data[1]
PCA.data[2]
#17


# Haz un plot con los PCA por muestras y por tipos



png(paste0(outdir,"/QC/PCA.groups.png"), width=1800, height=800)
DimPlot(data, reduction = "pca", group.by=c("orig.ident", "tag", "group", "batch"), raster = FALSE) 
dev.off()

#UMAP

data <- FindNeighbors(data, dims = 1:17)
data <- FindClusters(data, resolution = 0.6)

data <- RunUMAP(data, dims = 1:17)

png(paste0(outdir,"/QC/Umap.QC.png"), width=2400, height=800)
p1 <- FeaturePlot(data, features = "nCount_RNA", label=TRUE, raster = FALSE) & theme(plot.title = element_text(size=10))
p2 <- FeaturePlot(data, features = "nFeature_RNA", label=TRUE, raster = FALSE) & theme(plot.title = element_text(size=10))
p3 <- FeaturePlot(data, features = "percent.mt", label=TRUE, raster = FALSE) & theme(plot.title = element_text(size=10))
plot_grid(p1,p2,p3, ncol=3)
dev.off()



png(paste0(outdir,"/QC/Violin.QC.png"), width=1800, height=800)
plot1 <- VlnPlot(data, features = "nCount_RNA", pt.size = 0, group.by = "seurat_clusters", raster = FALSE) +
  xlab("cluster_id") +
  NoLegend()
plot2 <- VlnPlot(data, features = "nFeature_RNA", pt.size = 0, group.by = "seurat_clusters", raster = FALSE) +
  xlab("cluster_id") +
  NoLegend()
plot3 <- VlnPlot(data, features = "percent.mt", pt.size = 0, group.by = "seurat_clusters", raster = FALSE) +
  xlab("cluster_id") +
  NoLegend()
plot_grid(plot1, plot2, plot3, ncol = 1)
dev.off()

png(paste0(outdir,"/QC/Umap.data.png"), width=1800, height=2400)
p1 <- DimPlot(data, reduction = "umap", label = T, raster = FALSE)
p2 <- DimPlot(data, reduction = "umap", group.by = "group", raster = FALSE)
p3 <- DimPlot(data, reduction = "umap", group.by = "tag", raster = FALSE)
p4 <-DimPlot(data, reduction = "umap", group.by = "batch", raster = FALSE)
p5 <-DimPlot(data, reduction = "umap", group.by = "DsRed_status", raster = FALSE)
plot_grid(p1,p2,p3,p4,p5, ncol=2)
dev.off()


png(paste0(outdir,"/QC/PCA.clusters.png"), width=1200, height=600)
DimPlot(data, reduction = "pca", group.by = "seurat_clusters", label.size = 4, repel = TRUE, label = TRUE, raster = FALSE)
dev.off()


png(paste0(outdir,"/QC/Umap.tag.png"), width=2400, height=1200)
DimPlot(data, reduction = "umap", split.by = "tag", ncol=4, raster = FALSE)
dev.off()



png(paste0(outdir,"/QC/Umap.png"), width=1800, height=1200)
DimPlot(data, reduction = "umap", ncol=4, raster = FALSE)
dev.off()


# SingleR



## SingleR Cell identification


sce <- as.SingleCellExperiment(DietSeurat(data))

ref <- celldex::ImmGenData()

ref.main <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.main)

ref.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.fine)


data@meta.data$Cell_type.Image <- ref.main$pruned.labels
data@meta.data$Cell_type.Image.fine <- ref.fine$pruned.labels




data <- SetIdent(data, value = "Cell_type.Image")
png(paste0(outdir,"/QC/Umap_SingleR.Image.Fine.png"), width=2200, height=2400)
data <- SetIdent(data, value = "seurat_clusters")
p0 <- DimPlot_scCustom(data, reduction = "umap", split.by="orig.ident", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
data <- SetIdent(data, value = "Cell_type.Image")
p1 <- DimPlot_scCustom(data, reduction = "umap", split.by="orig.ident", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
p2 <- DimPlot(data, reduction = "umap", split.by="orig.ident",group.by = "group", raster = FALSE)
p3 <- DimPlot(data, reduction = "umap", split.by="orig.ident",group.by = "tag", raster = FALSE)
plot_grid(p0,p1,p2,p3, ncol=1)
dev.off() 




png(paste0(outdir,"/QC/Umap_SingleR.png"), width=1800, height=800)
DimPlot_scCustom(data, reduction = "umap", split.by="orig.ident", label=T, repel=T,
label.size=4, ggplot_default_colors=T) + NoLegend()
dev.off()

png(paste0(outdir,"/QC/Umap_Markers.png"), width=1800, height=1400)
FeaturePlot_scCustom(data, features= c("Ptprc","Csf3r", "S100a8", "C1qa", "Ace","Csf1r", "Cd3e", "Cd4",
"Cd8a", "Il1b", "Pi16","Pdgfra", "Pecam1", "Col1a1", "Col3a1", "Trem2", "Arg1", "Mrc1"), reduction = "umap")
dev.off()



saveRDS(data, paste0(rsdir,"objects/data.rds"))
data <- readRDS(paste0(rsdir,"objects/data.rds"))


pdf(paste0(outdir,"/QC/Umap_total.pdf"), width=12, height=8)
DimPlot_scCustom(data, reduction = "umap", label=T, repel=T,
label.size=4, ggplot_default_colors=T) + NoLegend()
dev.off()

# Highlight KO3-DsRedp cells in the UMAP plot
# Asegúrate de que la columna 'tag' esté presente en los metadatos

# Añade una nueva columna lógica para marcar las células KO3-DsRedp
data$highlight_KO3 <- ifelse(data$tag == "KO3-DsRedp", "KO3-DsRedp", "Other")

# Especificar los niveles del factor para controlar el orden y el color
data$highlight_KO3 <- factor(data$highlight_KO3, levels = c("Other", "KO3-DsRedp"))


# UMAP coloreado por esa nueva columna
png(paste0(outdir,"/QC/Umap_Markers2.png"), width=1800, height=1400)
DimPlot(data, group.by = "highlight_KO3", cols = c("gray", "red")) +
  ggtitle("KO3-DsRedp Cells")
dev.off()