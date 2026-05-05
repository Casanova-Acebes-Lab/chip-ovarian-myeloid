library(Seurat)
library(dplyr)
library(cowplot)
library(scCustomize)
library(harmony)
library(ggplot2)
library(RPresto)
library(celldex)
library(SingleCellExperiment)
library(ggrepel)
library(dplyr)
library(tibble)
library(clustree)
library(presto)
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
soupx_dir <- "/storage/scratch01/users/dcaceres/Sarai_Soupx"


source("R/Functions.R")



# Read the data  


data<- readRDS(paste0(rsdir,"objects/data.Clusterized1.rds"))



tam <- subset(data, Cluster %in% c(  "Monocytes",
  "IFN Mac",
  "Trem1|Ptgs2|Plaur|Clec4e Mac",
  "Mrc1|Sparc Mac",
  "Arg1|Spp1|Mmp12|Il1a Mac",
  "Nrp2|Actn1 Mac",
  "Mmp9|Ctsk Mac",
  "Fn1|Vegfa Mac",
  "MHCII|Siglec Mac",
  "MHCII|Mgl2 Mac",
  "Gas6|Folr2 Mac",
  "MHCII|Ccl8 Mac"))



options(future.globals.maxSize = 24 * 1024^3)
tam <- SCTransform(
  tam,
  assay = "RNA",
  vars.to.regress = "percent.mt",
  verbose = TRUE
)

tam <- RunPCA(tam)

PCA.data <- PCA(tam)
PCA.data[1]
PCA.data[2]
#17

tam <- FindNeighbors(tam, dims = 1:17)
tam <- FindClusters(tam, resolution = 0.5)
tam <- RunUMAP(tam, dims = 1:17)

tam <- SetIdent(tam, value = "seurat_clusters")
png(paste0(outdir,"/Macros/Umap1.png"), width=1800, height=1800)
DimPlot(tam, reduction = "umap", group.by=c("orig.ident", "DsRed", "group", "batch"), raster = FALSE) 
dev.off()



png(paste0(outdir,"/Macros/Umap2.png"), width=1200, height=1200)
DimPlot(tam, reduction = "umap", raster = FALSE, label=T, repel=T, label.size=4) + NoLegend()
dev.off()

png(paste0(outdir,"/Macros/Umap3.png"), width=2800, height=1800)
DimPlot(tam, reduction = "umap", raster = FALSE, label=T, repel=T, label.size=4, 
split.by = "tag", ncol = 4) + NoLegend()
dev.off()


png(paste0(outdir,"/Macros/Umap_Markers.png"), width=1800, height=1800)
FeaturePlot_scCustom(tam, features= c("Trem1","Arg1", "S100a8", "Emp1", "Ly6c2","Fn1", "Vegfa", "Cd4",
"Cd8a", "Il1b", "Nrp2","Vegfa", "Mmp9", "Ifit3", "Mrc1", "Gas6", "Ciita", "Siglece", "Ccl12", "H2-D1", "Col1a1", "Ccl8",
"Hsph1", "Sparc", "Folr2"), reduction = "umap")
dev.off()






# Encontrar marcadores
markers <- FindAllMarkers(
  object =tam,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)

head(markers)




top30_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 30, with_ties = FALSE) %>%
  ungroup()


as.data.frame(top30_per_cluster)

tam$Macros <- NA
tam$Macros[tam$seurat_clusters == 0] <- "Trem1|Ptgs2|Plaur|Clec4e Mac"
tam$Macros[tam$seurat_clusters == 1] <- "MHCII|Hsp Mac"
tam$Macros[tam$seurat_clusters == 2] <- "IFN Mac"
tam$Macros[tam$seurat_clusters == 3] <- "MHCII|Siglec Mac"
tam$Macros[tam$seurat_clusters == 4] <- "Emp1|Vegfa Mac"
tam$Macros[tam$seurat_clusters == 5] <- "Gas6|Folr2 Mac"
tam$Macros[tam$seurat_clusters == 6] <- "Arg1|Spp1|Mmp12|Il1a Mac"
tam$Macros[tam$seurat_clusters == 7] <- "Emp1|Vegfa Mac"
tam$Macros[tam$seurat_clusters == 8] <- "Monocytes"
tam$Macros[tam$seurat_clusters == 9] <- "Mmp9|Ctsk Mac"
tam$Macros[tam$seurat_clusters == 10] <- "Monocytes"
tam$Macros[tam$seurat_clusters == 11] <- "Gas6|Folr2 Mac"
tam$Macros[tam$seurat_clusters == 12] <- "Arg1|Spp1|Mmp12|Il1a Mac"





# Crear un data.frame usando los barcodes correctos
barcodes_clusters <- data.frame(
  Barcode = colnames(tam),  # las celdas
  Macros = tam$Macros,      # el cluster asignado a cada celda
  stringsAsFactors = FALSE
)

# Ver los primeros registros
head(barcodes_clusters)

# Guardar a archivo
write.table(barcodes_clusters, 
            file = paste0(outdir, "/Macros/tam_barcodes_macros.tsv"), 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)