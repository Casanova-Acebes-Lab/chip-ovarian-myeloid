
library(Seurat)
library(dplyr)
library(cowplot)
library(scCustomize)
library(harmony)
library(ggplot2)
library(RPresto)
library("scProportionTest")
library(ggrepel)
library(stringr)




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


# Read data

data <- readRDS(paste0(rsdir,"objects/data.Clusterized.Round2.rds"))


# Subsetting T cells

tcell <- subset(data, Clustering.Round2 == "Cd4 T cells" |
                        Clustering.Round2 == "Cd8 T cells" |
                        Clustering.Round2 == "NK")





png(paste0(outdir,"/Tcells/Umap.png"), width=1200, height=800)
 DimPlot(tcell, reduction = "umap", split.by = "group", label = TRUE)

dev.off()



# Subslustering



# Identificación de genes variables ya está hecha por SCTransform
# PCA
tcell <- RunPCA(tcell, features = VariableFeatures(tcell))

# Vecindad y clustering
tcell <- FindNeighbors(tcell, dims = 1:30) %>%
         FindClusters(resolution = 0.2)

# UMAP
tcell <- RunUMAP(tcell, dims = 1:30)

png(paste0(outdir,"/Tcells/Umap.subclustering.png"), width=1200, height=800)
DimPlot(tcell, reduction = "umap", label = TRUE)
dev.off()



# Establecer identidades con Seurat >=4
Idents(tcell) <- "seurat_clusters"

# Encontrar marcadores
markers <- FindAllMarkers(
  object = tcell,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

head(markers)




top30_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 20, with_ties = FALSE) %>%
  ungroup()


as.data.frame(top30_per_cluster)


png(paste0(outdir,"/Tcells/Umap.cd.png"), width=1600, height=3200)
FeaturePlot_scCustom(tcell, features = c("Cd8a", "Cd4" , "Fcer1g", "Pdcd1", "Tox", "Gzmb"), split.by = "DsRed")
dev.off()



## SingleR Cell identification


sce <- as.SingleCellExperiment(DietSeurat(tcell))

ref <- celldex::ImmGenData()

ref.main <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.main)

ref.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.fine)


tcell@meta.data$Cell_type.Image <- ref.main$pruned.labels
tcell@meta.data$Cell_type.Image.fine <- ref.fine$pruned.labels




tcell <- SetIdent(tcell, value = "Cell_type.Image")
png(paste0(outdir,"/Tcells/Umap_SingleR.Image.Fine.png"), width=2200, height=2400)
tcell <- SetIdent(tcell, value = "Clustering.Round2")
p0 <- DimPlot_scCustom(tcell, reduction = "umap", split.by="orig.ident", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
tcell <- SetIdent(tcell, value = "Cell_type.Image")
p1 <- DimPlot_scCustom(tcell, reduction = "umap", split.by="orig.ident", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
p2 <- DimPlot(tcell, reduction = "umap", split.by="orig.ident",group.by = "group", raster = FALSE)
p3 <- DimPlot(tcell, reduction = "umap", split.by="orig.ident",group.by = "tag", raster = FALSE)
plot_grid(p0,p1,p2,p3, ncol=1)
dev.off() 



tcell <- SetIdent(tcell, value = "Cell_type.Image.fine")
png(paste0(outdir,"/Tcells/Umap_SingleR.Image.Fine2.png"), width=2200, height=1400)
DimPlot_scCustom(tcell, reduction = "umap", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
dev.off()



iltck_genes <- c("Gzmb","Gzma","Prf1","Nkg7",
                 "Fcer1g","Tyrobp","Klrk1","Klrb1c",
                 "Il2rb","Il15ra",
                 "Xcl1","Ccl5")





tcell <- AddModuleScore(tcell, features = list(iltck_genes), name = "ILTCK_score")

png(paste0(outdir,"/Tcells/Umap_Signature.Iltck.png"), width=1200, height=800)
FeaturePlot_scCustom(tcell, features = "ILTCK_score1")

dev.off()



all_markers <- c(
  # CD8 exhausted / effector
  "Cd8a", "Cd8b1", "Pdcd1", "Lag3", "Havcr2", "Tox", "Gzmk", "Nkg7", "Prf1", "Ifng", "Ccl5",
  
  # Treg
  "Foxp3", "Il2ra", "Ikzf2", "Ccr8", "Lrrc32", "Itgae",
  
  # Naive / TCM
  "Tcf7", "Il7r", "Klf2", "S1pr1",
  
  # T activadas / TNF
  "Ccl1", "Cd70", "Tnfsf4", "Tnfsf8", "Tnfsf11",
  
  # NK-like / NKT
  "Klrb1c", "Klrk1", "Klre1", "Xcl1",
  
  # Proliferación
  "Mki67", "Cdk1", "Ccna2", "Birc5", "Stmn1",
  
  # NK puras
  "Ncr1", "Klrb1b", "Gzmc", "Styk1",
  
  # Th17 / γδ17
  "Il17a", "Il17f", "Rorc", "Il23r", "Ccr6"

)



png(paste0(outdir,"/Tcells/Umap.markers.png"), width=1800, height=3000)
FeaturePlot_scCustom(tcell, features = all_markers)
dev.off()



tcell@meta.data$Subcluster <- NA
tcell@meta.data$Subcluster[tcell$seurat_clusters == 0] <- "Cd8 Effector"
tcell@meta.data$Subcluster[tcell$seurat_clusters == 1] <- "Cd4 Treg"
tcell@meta.data$Subcluster[tcell$seurat_clusters == 2] <- "Cd4 Naive"
tcell@meta.data$Subcluster[tcell$seurat_clusters == 3] <- "Cd4 Activated"
tcell@meta.data$Subcluster[tcell$seurat_clusters == 4] <- "NK"
tcell@meta.data$Subcluster[tcell$seurat_clusters == 5] <- "Mki67+ Tcell"
tcell@meta.data$Subcluster[tcell$seurat_clusters == 6] <- "Cd8 Effector"
tcell@meta.data$Subcluster[tcell$seurat_clusters == 7] <- "NK"
tcell@meta.data$Subcluster[tcell$seurat_clusters == 8] <- "Cd4 Th17"


tcell <- SetIdent(tcell, value = "Subcluster")
pdf(paste0(outdir,"/Tcells/Umap.Clusterized.pdf"), width=12, height=8)
 DimPlot(tcell, reduction = "umap", label = TRUE)

dev.off()


tcell <- SetIdent(tcell, value = "Subcluster")
pdf(paste0(outdir,"/Tcells/Umap.Clusterized.group.pdf"), width=14, height=8)
 DimPlot(tcell, reduction = "umap", split.by = "group", label = TRUE)

dev.off()


tcell <- SetIdent(tcell, value = "Subcluster")
pdf(paste0(outdir,"/Tcells/Umap.Clusterized.DsRed.pdf"), width=14, height=8)
 DimPlot(tcell, reduction = "umap", split.by = "DsRed", label = TRUE)

dev.off()


## Saving annotation

t.annotation <- data.frame(
  barcode = colnames(tcell),
  Subclustering.T = tcell$Subcluster
)

write.csv(
  t.annotation,
  file = paste0(outdir,"/Tcells/subclustering_tcells_annotations.csv"),
  row.names = FALSE
)