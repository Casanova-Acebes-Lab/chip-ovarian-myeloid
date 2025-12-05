
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
library(clustree)
library("scProportionTest")
library(ggrepel)
library(dplyr)
library(tibble)






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


smpR065_DsRedP_KO3 <- read.data("SoupX/2.0/smR065/Seurat_smpR065_DsRedP-KO3.2.rds",
"DsRedP-KO3",
"KO",
"DsRedP", 
"smpR065",
"C",
"KO3"
)


smR039a_WT1_DsRedn <- read.data("SoupX/2.0/SmR039a.9/Seurat_WT1-DsRedn.2.rds",
"WT1-DsRedN",
"WT",
"DsRedN", 
"smR039a",
"A",
"WT1"
)

smR039a_WT1_DsRedp <- read.data("SoupX/2.0/SmR039a.9/Seurat_WT1-DsRedp.2.rds",
"WT1-DsRedP",
"WT",
"DsRedP", 
"smR039a",
"A",
"WT1"
)

smR039a_DsRedN_KO2 <- read.data("SoupX/2.0/SmR039a.9/Seurat_KO2-DsRedn.2.rds",
"DsRedN_KO2",
"KO",
"DsRedP", 
"smR039a",
"A",
"KO2"
)

smR039a_DsRedP_KO2 <- read.data("SoupX/2.0/SmR039a.9/Seurat_KO2-DsRedp.2.rds",
"DsRedP_KO2",
"KO",
"DsRedP", 
"smR039a",
"A",
"KO2"
)


smR039b_WT2_DsRedn <- read.data("SoupX/SmR039b.9/Seurat_WT2-DsRedn.2.rds",
"WT2_DsRedN",
"WT",
"DsRedN", 
"smR039b",
"B",
"WT2"
)

smR039b_WT2_DsRedp <- read.data("SoupX/SmR039b.9/Seurat_WT2-DsRedp.2.rds",
"WT2_DsRedP",
"WT",
"DsRedP", 
"smR039b",
"B",
"WT2"
)

smR039b_DsRedN_KO3 <- read.data("SoupX/2.0/SmR039b.9/Seurat_KO3-DsRedn.2.rds",
"DsRedN_KO3",
"KO",
"DsRedN", 
"smR039b",
"B",
"KO3"
)



# Merge

object_list <- list(smpR065_DsRedP_KO3, smR039a_WT1_DsRedn, smR039a_WT1_DsRedp,
smR039a_DsRedN_KO2, smR039a_DsRedP_KO2,
smR039b_WT2_DsRedn, smR039b_WT2_DsRedp, smR039b_DsRedN_KO3)
data <- Merge_Seurat_List(list_seurat = object_list)



data <- JoinLayers(data)


# Removing lowly expressed genes

gene_counts <- rowSums(GetAssayData(data, layer = "counts") > 0)
genes_keep <- names(gene_counts[gene_counts >= 20])
data <- subset(data, features = genes_keep)

# Removing doublets

data <- subset(data, subset = DF.class == "Doublet", invert=TRUE)

data$batch <- factor(data$batch, levels = c("A", "B", "C"))

# QC2 

data$percent.mt <- PercentageFeatureSet(data, pattern = "^mt-")

data <- SetIdent(data, value = "batch")

png(paste0(outdir,"/QC2/QC.Before_filtering.Complete.png"), width=1200, height=600)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0,raster=FALSE)
dev.off()


plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)


png(paste0(outdir,"/QC2/ScatterQC.Before_filtering.Complete.png"), width=1200, height=800)
plot_grid(plot1,plot2, ncol=2)
dev.off()


## By sample 

data <- SetIdent(data, value = "tag")

png(paste0(outdir,"/QC2/QC.Before_filtering.Complete.tag2.png"), width=1200, height=600)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0, raster = FALSE)
dev.off()


plot1 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)
plot2 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)


png(paste0(outdir,"/QC2/ScatterQC.Before_filtering.Complete.tag.png"), width=1200, height=800)
plot_grid(plot1,plot2, ncol=2)
dev.off()


data <- subset(data, subset = nFeature_RNA > 300 & nFeature_RNA < 7500 & percent.mt < 10)


data <- SetIdent(data, value = "batch")

png(paste0(outdir,"/QC2/QC2.After_filtering.Complete.png"), width=1200, height=800)
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size=0 , raster = FALSE)
dev.off()


plot3 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "percent.mt", raster = FALSE)
plot4 <- FeatureScatter(data, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", raster = FALSE)


png(paste0(outdir,"/QC2/ScatterQC2.After_filtering.Complete.png"), width=1200, height=800)
plot_grid(plot1,plot2,plot3,plot4, ncol=2)
dev.off()


# Analysis

##### Normalize

#data<- NormalizeData(data)
#data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 2000)

#data <- ScaleData(data, vars.to.regress = c("nCount_RNA"))



### Scale

options(future.globals.maxSize = 16 * 1024^3)  # 8 GiB

# Ahora correr SCTransform
data <- SCTransform(data, method = "glmGamPoi",
 verbose = TRUE, 
 vars.to.regress = c("percent.mt"),
 variable.features.n = 3000
 )



# Linear Reduction

data <- RunPCA(data, features = VariableFeatures(object = data))


png(paste0(outdir,"/QC2/ElbowPlot.Complete.png"), width=800, height=600)
ElbowPlot(data)
dev.off()


##### PCA

PCA.data <- PCA(data)
PCA.data[1]
PCA.data[2]
#14


# Haz un plot con los PCA por muestras y por tipos



png(paste0(outdir,"/QC2/PCA.groups.png"), width=1800, height=800)
DimPlot(data, reduction = "pca", group.by=c("orig.ident", "tag", "group", "batch"), raster = FALSE) 
dev.off()

#UMAP

data <- FindNeighbors(data, dims = 1:30)

# Eliminar clusters viejos de ese prefijo (si existen)
cols_old <- grep("^SCT_snn_res", colnames(data@meta.data), value = TRUE)
if(length(cols_old) > 0){
  data@meta.data[, cols_old] <- NULL
}

resolutions <- c(0.2, 0.4, 0.6, 0.8, 1.0)

for (res in resolutions) {
  data <- FindClusters(
    data,
    resolution = res,
    verbose = FALSE,
    graph.name = "SCT_snn"   # grafo correcto
  )
}

# Ahora sí puedes listar las columnas generadas
grep("^SCT_snn_res", colnames(data@meta.data), value = TRUE)


png(paste0(outdir,"/QC2/clustree.png"), width=2200, height=800)
clustree(data, prefix = "SCT_snn_res.")
dev.off()


data <- FindClusters(data, resolution = 0.8)

data <- RunUMAP(data, dims = 1:30)

png(paste0(outdir,"/QC2/Umap.QC2.png"), width=2400, height=800)
p1 <- FeaturePlot(data.harmony, features = "nCount_RNA", label=TRUE, raster = FALSE) & theme(plot.title = element_text(size=10))
p2 <- FeaturePlot(data.harmony, features = "nFeature_RNA", label=TRUE, raster = FALSE) & theme(plot.title = element_text(size=10))
p3 <- FeaturePlot(data.harmony, features = "percent.mt", label=TRUE, raster = FALSE) & theme(plot.title = element_text(size=10))
plot_grid(p1,p2,p3, ncol=3)
dev.off()



png(paste0(outdir,"/QC2/Violin.QC2.png"), width=1800, height=800)
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

png(paste0(outdir,"/QC2/Umap.data.png"), width=1800, height=2400)
p1 <- DimPlot(data, reduction = "umap", label = T, raster = FALSE)
p2 <- DimPlot(data, reduction = "umap", group.by = "group", raster = FALSE)
p3 <- DimPlot(data, reduction = "umap", group.by = "tag", raster = FALSE)
p4 <-DimPlot(data, reduction = "umap", group.by = "batch", raster = FALSE)
p5 <-DimPlot(data, reduction = "umap", group.by = "DsRed", raster = FALSE)
plot_grid(p1,p2,p3,p4,p5, ncol=2)
dev.off()


png(paste0(outdir,"/QC2/PCA.clusters.png"), width=1200, height=600)
DimPlot(data, reduction = "pca", group.by = "seurat_clusters", label.size = 4, repel = TRUE, label = TRUE, raster = FALSE)
dev.off()



data <- SetIdent(data, value = "seurat_clusters")
png(paste0(outdir,"/QC2/Umap.tag.png"), width=2400, height=1200)
DimPlot(data, reduction = "umap", split.by = "tag", ncol=4, raster = FALSE,
pt.size = 0.5)       
 
dev.off()



png(paste0(outdir,"/QC2/Umap.png"), width=1200, height=1200)
DimPlot(data, reduction = "umap", ncol=4, raster = FALSE)
dev.off()



#### Checking batch effect


 
    # Aseguramos que las columnas existen en el metadata
    meta <- data@meta.data %>%
      dplyr::select(
        nCount_RNA,
        nFeature_RNA,
        percent.mt,
        batch,
        tag,
        group,
        chimera
      )
 
    # Crear función para resumir
    summ_fun <- function(df, group_vars) {
      df %>%
        group_by(across(all_of(group_vars))) %>%
        summarise(
          nCells = n(),
          mean_nCount = mean(nCount_RNA),
          median_nCount = median(nCount_RNA),
          sd_nCount = sd(nCount_RNA),
 
          mean_nFeature = mean(nFeature_RNA),
          median_nFeature = median(nFeature_RNA),
          sd_nFeature = sd(nFeature_RNA),
 
          mean_percentMT = mean(percent.mt),
          median_percentMT = median(percent.mt),
          sd_percentMT = sd(percent.mt),
 
          .groups = "drop"
        )
    }
 
    # Tabla agrupada por batch, tag, group y chimera
    summary_all <- summ_fun(meta, c("batch", "tag", "group", "chimera"))
    as.data.frame(summary_all)



saveRDS(data, paste0(rsdir,"objects/data.PreSingleR.rds"))

data <-readRDS(paste0(rsdir,"objects/data.PreSingleR.rds"))


data.harmony <- RunHarmony(
data, 
group.by.vars = "batch", 
dims = 1:20,
max.iter.harmony = 50)  



data.harmony <- RunUMAP(data.harmony, reduction = "harmony", dims = 1:20)
data.harmony <- FindNeighbors(data.harmony, reduction = "harmony", dims = 1:20)
data.harmony <- FindClusters(data.harmony, resolution = 0.8)



png(paste0(outdir,"/QC2/Umap.harmony.png"), width=1200, height=1200)
DimPlot(data.harmony, raster = FALSE)
dev.off()


png(paste0(outdir,"/QC2/Umap.harmony.tag.png"), width=2800, height=1200)
DimPlot(data.harmony, raster = FALSE,split.by = "tag", ncol=4)
dev.off()




# Downsampling

set.seed(123)

n_cells <- 8000
cells_to_keep <- c()

for (group in unique(data.harmony$tag)) {
  # Subset de células por metadata
  group_cells <- rownames(data.harmony@meta.data)[data.harmony@meta.data$tag == group]
  
  if (length(group_cells) <= n_cells) {
    sampled_cells <- group_cells
  } else {
    sampled_cells <- sample(group_cells, n_cells)
  }
  
  cells_to_keep <- c(cells_to_keep, sampled_cells)
}

# Ahora sí crear el subset
data_downsampled <- subset(data.harmony, cells = cells_to_keep)

# Revisar número de células por grupo
table(data_downsampled$tag)



png(paste0(outdir,"/QC2/Umap.harmony.tag.downsampled.8k.png"), width=2800, height=1200)
DimPlot(data_downsampled, raster = FALSE,split.by = "tag", ncol=4)
dev.off()


data_downsampled <- SetIdent(data_downsampled, value = "Clustering.Round2")
pdf(paste0(outdir,"/QC2/Umap.harmony.chimera.downsampled.8k.pdf"), width=18, height=12)
DimPlot(data_downsampled, raster = FALSE,split.by = "chimera", ncol=2, cols= nora.colors)
dev.off()

png(paste0(outdir,"/QC2/Umap.harmony.group.downsampled.8k.png"), width=1000, height=800)
DimPlot(data_downsampled, raster = FALSE,group.by = "group")
dev.off()



png(paste0(outdir,"/QC2/violin.Il1b.down.png"), width=3800, height=1000)
VlnPlot(
  data_downsampled,
  features = "Il1b",
  group.by = "tag",
  split.by = "seurat_clusters",
  pt.size = 0
)
dev.off()


 pdf(paste0(outdir,"/QC2/Umap_Markers.Il.pdf"), width=16, height=16)
FeaturePlot_scCustom(data_downsampled, 
features = c("Il1b", "Il1a" , "Il16"), split.by="group")
dev.off()



 pdf(paste0(outdir,"/QC2/Umap_Markers.Il.ch.pdf"), width=32, height=16)
FeaturePlot_scCustom(data_downsampled, 
features = c("Il1b", "Il1a" , "Il16"), split.by="chimera")
dev.off()



# SingleR



## SingleR Cell identification


sce <- as.SingleCellExperiment(DietSeurat(data.harmony, assays = c("SCT")))
ref <- celldex::ImmGenData()

ref.main <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.main)

ref.fine <- SingleR(test = sce,assay.type.test = 1,ref = ref,labels = ref$label.fine)


data.harmony@meta.data$Cell_type.Image <- ref.main$pruned.labels
data.harmony@meta.data$Cell_type.Image.fine <- ref.fine$pruned.labels




data.harmony <- SetIdent(data.harmony, value = "Cell_type.Image")
png(paste0(outdir,"/QC2/Umap_SingleR.Image.Fine.png"), width=2200, height=2400)
data.harmony <- SetIdent(data.harmony, value = "seurat_clusters")
p0 <- DimPlot_scCustom(data.harmony, split.by="orig.ident", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
data <- SetIdent(data.harmony, value = "Cell_type.Image")
p1 <- DimPlot_scCustom(data.harmony, split.by="orig.ident", label=T, repel=T,
label.size=4, ggplot_default_colors=T, raster = FALSE) + NoLegend()
p2 <- DimPlot(data.harmony, split.by="orig.ident",group.by = "group", raster = FALSE)
p3 <- DimPlot(data.harmony, split.by="orig.ident",group.by = "tag", raster = FALSE)
plot_grid(p0,p1,p2,p3, ncol=1)
dev.off() 



saveRDS(data.harmony, paste0(rsdir,"objects/data.QC.harmony.SingleR.rds"))
data.harmony <- readRDS(paste0(rsdir,"objects/data.QC.harmony.SingleR.rds"))



data.harmony <- SetIdent(data.harmony, value = "Cell_type.Image.fine")
png(paste0(outdir,"/QC2/Umap_SingleR.png"), width=1800, height=800)
DimPlot_scCustom(data.harmony, reduction = "umap", label=T, repel=T,
label.size=4, ggplot_default_colors=T) + NoLegend()
dev.off()

png(paste0(outdir,"/QC2/Umap_Markers.png"), width=1800, height=1400)
FeaturePlot_scCustom(data.harmony, features= c("Ptprc","Csf3r", "S100a8", "C1qa", "Ace","Csf1r", "Cd3e", "Cd4",
"Cd8a", "Il1b", "Pi16","Pdgfra", "Pecam1", "Col1a1", "Ly6c2", "Trem2", "Arg1", "Mrc1"), reduction = "umap")
dev.off()




png(paste0(outdir,"/QC2/Umap_Markers.split.png"), width=3200, height=900)
FeaturePlot_scCustom(data.harmony, features= "Il1b", reduction = "umap", split.by="chimera")
dev.off()



pdf(paste0(outdir,"/QC2/Umap_Markers.pdf"), width=18, height=14)
FeaturePlot_scCustom(data.harmony, features= c("Ptprc","Csf3r", "Cd9", "Spp1", "Cd163","Csf1r", "H2-Aa", "Adgre1",
"H2-Eb1", "Il1b", "Mki67","Mmp9", "Ccr2", "Cd3e", "Vcan", "Ly6c2", "Arg1", "Mrc1", "Adgre4", "Gpnmb"), reduction = "umap")
dev.off()



pdf(paste0(outdir,"/QC2/Umap_Markers.IFN.pdf"), width=14, height=24)
FeaturePlot_scCustom(data.harmony, features= c("Ifit3","Ifit1", "Ifit2", "Il1b", "Oas2", "Irf7"), reduction = "umap",
split.by="group")
dev.off()


 pdf(paste0(outdir,"/QC2/Umap_Markers.tcells.pdf"), width=10, height=10)
FeaturePlot_scCustom(data.harmony, 
features = c("Cd8a", "Gzmb" , "Cd8b1", "Prf1", "Cd3e", "Cd4")) & theme(plot.title = element_text(size=10))
 dev.off()





 pdf(paste0(outdir,"/QC2/Umap_Markers.tcells2.pdf"), width=12, height=22)
FeaturePlot_scCustom(data.harmony, 
features = c("Cd8a", "Gzmb" , "Cd8b1", "Prf1", "Cd3e", "Cd4"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()





 png(paste0(outdir,"/QC2/Umap_Markers.scratch.png"), width=1200, height=2600)
FeaturePlot_scCustom(data.harmony, 
features = c("Anxa1", "Actn1", "Npr2", "Arg1"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()


 png(paste0(outdir,"/QC2/Umap_Markers.3.png"), width=1200, height=2800)
FeaturePlot_scCustom(data.harmony, 
features = c("Arg1", "Mmp12" , "Inhba", "Slc7a2"
, "Spp1", "Sdc1", "Trem1", "Nrp1", "Il1a"
), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()


## DEG


# Establecer identidades con Seurat >=4
Idents(data.harmony) <- "seurat_clusters"

# Encontrar marcadores
markers <- FindAllMarkers(
  object =data.harmony,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

head(markers)




top30_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 50, with_ties = FALSE) %>%
  ungroup()


as.data.frame(top30_per_cluster)


write.table(top30_per_cluster, paste0(rsdir,"table.clusters.markers.top30.seurat_clusters.harmony.tsv"), sep='\t')






data.harmony$Cluster <- NA
data.harmony$Cluster[data.harmony$seurat_clusters == 0] <- "MHCII|H2-D Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 1] <- "MHCII|Siglec Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 2] <- "Mrc1|C1qc|Cbr2|Gas6 Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 3] <- "Mrc1|C1qc|Cbr2|Gas6 Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 4] <- "Cd8 T cells"
data.harmony$Cluster[data.harmony$seurat_clusters == 5] <- "Trem1|Ptgs2|Plaur|Celc4e Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 6] <- "Fn1|Vegfa Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 7] <- "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 8] <- "IFN Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 9] <- "Mrc1|C1qc|Cbr2|Gas6 Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 10] <- "Ly6c|Ms4ac Monocytes"
data.harmony$Cluster[data.harmony$seurat_clusters == 11] <- "Cd4"
data.harmony$Cluster[data.harmony$seurat_clusters == 12] <- "Ly6c|Ms4ac Monocytes"
data.harmony$Cluster[data.harmony$seurat_clusters == 13] <- "DCs" 
data.harmony$Cluster[data.harmony$seurat_clusters == 14] <- "Neutrophils"
data.harmony$Cluster[data.harmony$seurat_clusters == 15] <- "Npr2|Actn1 Mac" 
data.harmony$Cluster[data.harmony$seurat_clusters == 16] <- "Cd4"
data.harmony$Cluster[data.harmony$seurat_clusters == 17] <- "NK" 
data.harmony$Cluster[data.harmony$seurat_clusters == 18] <- "Remove"
data.harmony$Cluster[data.harmony$seurat_clusters == 19] <- "B cells" 
data.harmony$Cluster[data.harmony$seurat_clusters == 20] <- "Cd8" 
data.harmony$Cluster[data.harmony$seurat_clusters == 21] <- "Remove"
data.harmony$Cluster[data.harmony$seurat_clusters == 22] <- "Mmp9|Ctsk Mac"  
data.harmony$Cluster[data.harmony$seurat_clusters == 23] <- "DCs"
data.harmony$Cluster[data.harmony$seurat_clusters == 24] <- "Mmp9|Ctsk Mac"
data.harmony$Cluster[data.harmony$seurat_clusters == 25] <- "Activated B cells"
data.harmony$Cluster[data.harmony$seurat_clusters == 26] <- "Mastocytes"




data.harmony <- subset(data.harmony, subset = (Cluster == "Remove"), invert=TRUE)

cluster_order <- c(
  # --- MONOCITOS & TAMs (ROJOS) ---
  "Ly6cHi Monocytes",
  "Ly6cLo Monocytes",
  "Early IFN|MHCII-TAMs",
  "Nrg1|Cdh1 Mac",
  "Trem1|Ptgs2|Plaur|Celc4e Mac",
  "Mrc1|C1qc|Cbr2|Gas6 Mac",   # amarillo brillante (TAM residentes)
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac", # verde TAM reparadores
  "Npr2|Actn1 Mac",   # verde lima (TAM estructurales)
  "Mmp9|Ctsk Mac",   # verde botella (remodelado matriz)
  "IFN Mac",   # morado claro (TAM interferón maduros)
  "Fn1|Vegfa Mac",   # naranja angiogénico
  "MHCII|Siglec Mac",   # azul cobalto (antigen presenting)
  "MHCII|H2-D Mac",   # azul acero

  
  "Cd8",   # turquesa
  "Cd4",   # rosa fuerte
  "NK",   # violeta oscuro
  "Activated B cells",   # fucsia intenso
  "B cells",   # rojo intenso

  # --- INNATAS ---
  "Neutrophils",   # azul fuerte
  "DCs",   # azul claro
  "Mastocytes"    # verde saturado
)




nora.colors <- c(
  # --- MONOCITOS & TAMs (ROJOS) ---
  "Ly6cHi Monocytes"     = "#FF0000",  # Rojo inflamatorio puro
  "Ly6cLo Monocytes"     = "#FF6A6A",  # Rojo salmón transición
  "Early IFN|MHCII-TAMs"       = "#B22222",  # Rojo vino (pre-TAM IFN)
  "Nrg1|Cdh1 Mac"              =  "#5F9EA0",
  "Trem1|Ptgs2|Plaur|Celc4e Mac" = "#EE7942",   # naranja vivo (inflam. activados)
  "Mrc1|C1qc|Cbr2|Gas6 Mac"      = "#FFD92F",   # amarillo brillante (TAM residentes)
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac" = "#4DAF4A", # verde TAM reparadores
  "Npr2|Actn1 Mac"       = "#A6D854",   # verde lima (TAM estructurales)
  "Mmp9|Ctsk Mac"        = "#00723F",   # verde botella (remodelado matriz)
  "IFN Mac"              = "#C080FF",   # morado claro (TAM interferón maduros)
  "Fn1|Vegfa Mac"        = "#FFA500",   # naranja angiogénico
  "MHCII|Siglec Mac"     = "#1E90FF",   # azul cobalto (antigen presenting)
  "MHCII|H2-D Mac"      = "#4682B4",   # azul acero

  # --- CÉLULAS LINFOIDES ---
  "Cd8"          = "#00BFC4",   # turquesa
  "Cd4"          = "#FF69B4",   # rosa fuerte
  "NK"                   = "#AB82FF",   # violeta oscuro
  "Activated B cells"    = "#FF1493",   # fucsia intenso
  "B cells"              = "#DC143C",   # rojo intenso

  # --- INNATAS ---
  "Neutrophils"          = "#4876FF",   # azul fuerte
  "DCs"                  = "#87CEEB",   # azul claro
  "Mastocytes"           = "#3CB371"    # verde saturado
)




data.harmony$Cluster <- factor(data.harmony$Cluster, levels = cluster_order)




saveRDS(data.harmony, paste0(rsdir,"objects/data.harmony.Clusterized.rds"))
data.harmony <- readRDS(paste0(rsdir,"objects/data.harmony.Clusterized.rds"))



### Adding monocytes subclustering

mon.annotation <- read.csv(paste0(outdir,"/Subclustering.Monocytes2/subclustering_monocytes_annotations.csv"))

# Asegura que el índice son los barcodes
rownames(mon.annotation) <- mon.annotation$barcode
mon.annotation$barcode <- NULL

# Agregar metadata al objeto global
data.harmony<- AddMetaData(data.harmony, mon.annotation)



# Crear la nueva columna como character para que permita nuevos valores
data.harmony$Clustering.Round2 <- as.character(data.harmony$Cluster)

# Reasignar según tus subclusters refinados
data.harmony$Clustering.Round2[data.harmony$Subclustering.Mon == "Ly6cHi Monocytes"] <- "Ly6cHi Monocytes"
data.harmony$Clustering.Round2[data.harmony$Subclustering.Mon == "Ly6cLo Monocytes"] <- "Ly6cLo Monocytes"
data.harmony$Clustering.Round2[data.harmony$Subclustering.Mon == "Early IFN-TAMs"] <- "Early IFN|MHCII-TAMs"
data.harmony$Clustering.Round2[data.harmony$Subclustering.Mon == "Nrg1|Cdh1 Macs"] <- "Trem1|Ptgs2|Plaur|Celc4e Mac"


# Convertimos de vuelta a factor, ya con todos los niveles presentes
data.harmony$Clustering.Round2 <- factor(data.harmony$Clustering.Round2)


#####


### Subseting macros for Trayectory analysis



macros_cells <- WhichCells(
 data.harmony,
  expression = Clustering.Round2 %in% c(
    "Ly6cHi Monocytes", "Ly6cLo Monocytes", "Early IFN|MHCII-TAMs","Nrg1|Cdh1 Mac",
    "Trem1|Ptgs2|Plaur|Celc4e Mac",
    "Mrc1|C1qc|Cbr2|Gas6 Mac",
    "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac",
    "Npr2|Actn1 Mac",
    "Mmp9|Ctsk Mac",
    "IFN Mac",
    "Fn1|Vegfa Mac",
    "MHCII|Siglec Mac",
    "MHCII|H2-D Mac"
  )
)

macros <- subset(data.harmony, cells = macros_cells)


macros <- SetIdent(macros, value = "Clustering.Round2")

pdf(paste0(outdir,"/QC2/Clustering.Macros.pdf"), width=14, height=8)
DimPlot_scCustom(macros, reduction = "umap", group.by = "Clustering.Round2", pt.size = 0.5, label=FALSE) +
  scale_color_manual(values = nora.colors)
dev.off()

# Extraer las coordenadas de UMAP
umap_coords <- Embeddings(macros, "umap")

# Seleccionar células a conservar: solo aquellas con UMAP1 <= 5
celdas_a_conservar <- rownames(umap_coords)[umap_coords[,1] <= 5]

# Crear objeto filtrado
macros_filtrado <- subset(macros, cells = celdas_a_conservar)

# Resetear identidades a la columna de clusters
macros_filtrado <- SetIdent(macros_filtrado, value = "Clustering.Round2")

# Plot UMAP
pdf(paste0(outdir,"/QC2/Clustering.Macros2.pdf"), width=14, height=8)
DimPlot_scCustom(macros_filtrado, group.by = "Clustering.Round2",
                 reduction = "umap", pt.size = 0.5, label=FALSE) +
  scale_color_manual(values = nora.colors)
dev.off()



# Extraer coordenadas de UMAP
umap_coords <- Embeddings(macros_filtrado, "umap")

# Obtener las células que cumplen la condición de estar a la derecha de 3 y por encima de 7
celdas_a_eliminar <- rownames(umap_coords)[umap_coords[,1] > 3 & umap_coords[,2] > 7]

# Seleccionar las células que NO están en esa región
celdas_a_conservar <- setdiff(Cells(macros_filtrado), celdas_a_eliminar)

# Crear objeto filtrado
macros_filtrado <- subset(macros_filtrado, cells = celdas_a_conservar)

# Resetear identidades
macros_filtrado <- SetIdent(macros_filtrado, value = "Clustering.Round2")

# Graficar UMAP
pdf(paste0(outdir,"/QC2/Clustering.Macros3.pdf"), width=14, height=8)
DimPlot_scCustom(macros_filtrado, group.by = "Clustering.Round2",
                 reduction = "umap", pt.size = 0.5, label=FALSE) +
  scale_color_manual(values = nora.colors)
dev.off()





# Extraer coordenadas UMAP
umap_coords <- Embeddings(macros_filtrado, "umap")  # usa el objeto filtrado actual

# Seleccionar células a conservar: fuera de la nueva zona
celdas_a_conservar <- rownames(umap_coords)[
  !(umap_coords[,1] > 2.5 & umap_coords[,2] > 2)
]

# Crear un nuevo objeto Seurat filtrado
macros_filtrado <- subset(macros_filtrado, cells = celdas_a_conservar)

# Ajustar identidades si es necesario
macros_filtrado <- SetIdent(macros_filtrado, value = "Clustering.Round2")

# Graficar y guardar en PDF
pdf(paste0(outdir, "/QC2/Clustering.Macros4.pdf"), width=14, height=8)
DimPlot_scCustom(macros_filtrado, reduction = "umap", pt.size = 0.5, label = FALSE) +
  scale_color_manual(values = nora.colors)
dev.off()


# Extraer coordenadas UMAP
umap_coords <- Embeddings(macros_filtrado, "umap")  # usa el objeto filtrado actual

# Seleccionar células a conservar: fuera de la nueva zona
celdas_a_conservar <- rownames(umap_coords)[
  !(umap_coords[,1] > 2.5 & umap_coords[,2] < -6)
]

# Crear un nuevo objeto Seurat filtrado
macros_filtrado <- subset(macros_filtrado, cells = celdas_a_conservar)

# Ajustar identidades si es necesario
macros_filtrado <- SetIdent(macros_filtrado, value = "Clustering.Round2")

# Graficar y guardar en PDF
pdf(paste0(outdir, "/QC2/Clustering.Macros5.pdf"), width=14, height=8)
DimPlot_scCustom(macros_filtrado, reduction = "umap", pt.size = 0.5, label = FALSE) +
  scale_color_manual(values = nora.colors)
dev.off()




library(ggplot2)
library(Seurat)

# Extraer coordenadas
umap_coords <- Embeddings(macros_filtrado, "umap")
umap_df <- as.data.frame(umap_coords)
umap_df$Cluster <- Idents(macros_filtrado)

# DimPlot estilo ggplot con cuadrícula
pdf(paste0(outdir, "/QC2/Clustering.Macros5.pdf"), width=14, height=8)
ggplot(umap_df, aes(x = umap_1, y = umap_2, color = Cluster)) +
  geom_point(size = 0.5) +
  theme_bw() +
  geom_hline(yintercept = seq(floor(min(umap_df$umap_2, na.rm = TRUE)), 
                              ceiling(max(umap_df$umap_2, na.rm = TRUE)), by = 1), 
             color = "gray80", linetype = "dashed") +
  geom_vline(xintercept = seq(floor(min(umap_df$umap_1, na.rm = TRUE)), 
                              ceiling(max(umap_df$umap_1, na.rm = TRUE)), by = 1), 
             color = "gray80", linetype = "dashed")


dev.off()



# Extraer coordenadas UMAP
umap_coords <- Embeddings(macros_filtrado, "umap")

# Extraer identidades actuales
clusters <- Idents(macros_filtrado)

# Identificar las células en la zona definida
zona_a_eliminar <- rownames(umap_coords)[
  umap_coords[,1] >= 1 & umap_coords[,1] <= 3 &
  umap_coords[,2] >= 2 & umap_coords[,2] <= 5 &
  !(clusters %in% c("Ly6cHi Monocytes", "Ly6cLo Monocytes", "Early IFN|MHCII-TAMs", "IFN Mac"))
]


# Seleccionar las células a conservar (todas las demás)
celdas_a_conservar <- setdiff(rownames(umap_coords), zona_a_eliminar)

# Crear un nuevo objeto Seurat filtrado
macros_filtrado <- subset(macros_filtrado, cells = celdas_a_conservar)

# Ajustar identidades si es necesario
macros_filtrado <- SetIdent(macros_filtrado, value = "Clustering.Round2")

# Graficar y guardar en PDF
pdf(paste0(outdir, "/QC2/Clustering.Macros7.pdf"), width = 14, height = 8)
DimPlot_scCustom(macros_filtrado, reduction = "umap", pt.size = 0.5, label = FALSE) +
  scale_color_manual(values = nora.colors)
dev.off()




saveRDS(macros_filtrado, paste0(rsdir,"objects/data.macrophages.Clusterized2.rds"))
macros <- readRDS(paste0(rsdir,"objects/data.macrophages.Clusterized.rds"))


saveRDS(data, paste0(rsdir,"objects/data.Clusterized.Round2.rds"))
data <- readRDS(paste0(rsdir,"objects/data.Clusterized.Round2.rds"))





data <- SetIdent(macros, value = "Clustering.Round2")
pdf(paste0(outdir,"/QC2/Clustering.Macro8.pdf"), width=14, height=8)
DimPlot_scCustom(macros, reduction = "umap", pt.size = 0.5, label=FALSE) +
  scale_color_manual(values = nora.colors) 

dev.off()


## Subsetting macros for Velocity


macros <- macros_filtrado

macros$barcode <- colnames(macros)
macros$UMAP_1 <- macros@reductions$umap@cell.embeddings[,1]
macros$UMAP_2 <- macros@reductions$umap@cell.embeddings[,2]
write.csv(macros@meta.data, file = paste0(rsdir, "/objects/scVelo/metadata.csv"), 
          quote = FALSE, row.names = FALSE)


# write expression counts matrix
library(Matrix)
counts_matrix <- GetAssayData(macros, assay='SCT', slot='counts')
writeMM(counts_matrix, file=paste0(rsdir, '/objects/scVelo/counts.mtx'))

# write dimesnionality reduction matrix, in this example case pca matrix

write.csv(macros@reductions$pca@cell.embeddings, file = paste0(rsdir, "/objects/scVelo/pca.csv"), 
          quote = FALSE, row.names = FALSE)

# write gene names
write.table(
  data.frame('gene'=rownames(counts_matrix)),file=paste0(rsdir,'/objects/scVelo/gene_names.csv'),
  quote=F,row.names=F,col.names=F
)




nora.colors <- c(
  # --- MONOCITOS & TAMs (ROJOS) ---
  "Ly6cHi Monocytes"     = "#FF0000",  # Rojo inflamatorio puro
  "Ly6cLo Monocytes"     = "#FF6A6A",  # Rojo salmón transición
  "Early IFN|MHCII-TAMs"       = "#B22222",  # Rojo vino (pre-TAM IFN)
  "Trem1|Ptgs2|Plaur|Celc4e Mac" = "#EE7942",   # naranja vivo (inflam. activados)
  "Mrc1|C1QC2|Cbr2|Gas6 Mac"      = "#FFD92F",   # amarillo brillante (TAM residentes)
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac" = "#4DAF4A", # verde TAM reparadores
  "Npr2|Actn1 Mac"       = "#A6D854",   # verde lima (TAM estructurales)
  "Marco+ Mac"           = "#1E3A8A",   # azul intermedio (scavenger)
  "Mmp9|Ctsk Mac"        = "#00723F",   # verde botella (remodelado matriz)
  "IFN Mac"              = "#C080FF",   # morado claro (TAM interferón maduros)
  "Fn1|Vegfa Mac"        = "#FFA500",   # naranja angiogénico
  "MHCII|Siglec Mac"     = "#1E90FF",   # azul cobalto (antigen presenting)
  "MHCII|Ccl12 Mac"      = "#4682B4",   # azul acero

  # --- CÉLULAS LINFOIDES ---
  "Cd8 T cells"          = "#00BFC4",   # turquesa
  "Cd4 T cells"          = "#FF69B4",   # rosa fuerte
  "Cd8 Effector"         = "#A52A2A",   # rojo ladrillo
  "Tgd"                  = "#D9B3FF",   # lila pastel
  "NK"                   = "#AB82FF",   # violeta oscuro
  "Activated B cells"    = "#FF1493",   # fucsia intenso
  "B cells"              = "#DC143C",   # rojo intenso

  # --- INNATAS ---
  "Neutrophils"          = "#4876FF",   # azul fuerte
  "DCs"                  = "#87CEEB",   # azul claro
  "Mastocytes"           = "#3CB371"    # verde saturado
)




data$Cluster <- factor(data$Clustering.Round2, levels = names(nora.colors))

data <- SetIdent(data, value = "Clustering.Round2")
pdf(paste0(outdir,"/QC2/Clustering1.pdf"), width=12, height=8)
DimPlot_scCustom(data, reduction = "umap", label = TRUE, repel = TRUE, pt.size = 0.5) +
  scale_color_manual(values = nora.colors) +
  NoLegend()

dev.off()



data <- SetIdent(data, value = "Clustering.Round2")
pdf(paste0(outdir,"/QC2/Clustering2.pdf"), width=14, height=8)
DimPlot_scCustom(data, reduction = "umap", pt.size = 0.5, label=FALSE) +
  scale_color_manual(values = nora.colors) 

dev.off()



### Markers

 pdf(paste0(outdir,"/QC2/Umap_Markers.Pattern recognition receptors.pdf"), width=12, height=24)
FeaturePlot_scCustom(data, 
features = c("Marco", "Cd163", "Clec9a", "Tlr2", "Tlr4", "Cd93"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()




 pdf(paste0(outdir,"/QC2/Umap_Markers.ECM remodeling.pdf"), width=12, height=20)
FeaturePlot_scCustom(data, 
features = c("Mmp12", "Timp2", "Mmp8", "Hpse"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()




 pdf(paste0(outdir,"/QC2/Umap_Markers.LAMs.pdf"), width=12, height=24)
FeaturePlot_scCustom(data, 
features = c("Lgal3", "Cd36", "Fabp5", "Cd9", "Trem2"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()


 pdf(paste0(outdir,"/QC2/Umap_Markers.Immune suppression.pdf"), width=12, height=16)
FeaturePlot_scCustom(data, 
features = c("Nos2", "Arg1", "Vegfa"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()



 pdf(paste0(outdir,"/QC2/Umap_Markers.Fibrotic Macs.pdf"), width=12, height=24)
FeaturePlot_scCustom(data, 
features = c("Tgfb2", "Tgfb1", "Acta2", "Acta2", "Col1a1", "Timp1"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()


 pdf(paste0(outdir,"/QC2/Umap_Markers.TREM1 Macs.pdf"), width=12, height=18)
FeaturePlot_scCustom(data, 
features = c("Il1b", "Il1a", "Nlrp3", "Trem1"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()


 pdf(paste0(outdir,"/QC2/Umap_Markers.SAMs.pdf"), width=12, height=24)
FeaturePlot_scCustom(data, 
features = c("Spp1", "Lgals3", "Tnfsf12", "Pdgfb", "Vegfa"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()


 pdf(paste0(outdir,"/QC2/Umap_Markers.Resident Macs.pdf"), width=12, height=16)
FeaturePlot_scCustom(data, 
features = c("Cd163", "Timd4", "Lyve1", "Folr2"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()


 pdf(paste0(outdir,"/QC2/Umap_Markers.Ag-presenting Macs.pdf"), width=12, height=26)
FeaturePlot_scCustom(data, 
features = c("H2-Aa", "H2-Ab1" ,"H2-K1", "H2-Eb1", "Cd74", "Ciita", "B2m"), split.by="group") & theme(plot.title = element_text(size=10))
 dev.off()




## Rambo signature

rambo <- c(
  "Alox15","Tgfb2","Cav1","F5","Ccl24","Serpine1","Ltbp1","Cxcl13","Ankrd1","Garnl3",
  "Prg4","Serpinb2","Tgfbr3","Kank2","Tmem176b","Me1","Cd55","Plxdc2","Selp","Socs2",
  "Nt5e","Htra3","Cxcl2","Vmn2r26","Ccn1","Thbs1","Epha4","Ptgs2","Olr1","Slc7a11",
  "Id1","Flnb","Ccn4","Map1a","Cxcl5","Kcnq5","Kif26b","Akap12","Fam43a","Ppbp",
  "Il1b","Ofcc1","Tiam2","Csf1","Cxcl1","Timp3","Acpp","Plscr2","Olfml3","Aldh1a1"
)




# Asegúrate de que 'rambo' contiene genes presentes en tus datos
genes_present <- rambo[rambo %in% rownames(data)]

# Crear lista de genes (AddModuleScore requiere lista)
gene_list <- list(genes_present)

# Añadir el score al objeto Seurat
# name = "RAMBO" → creará columna 'RAMBO1' en data@meta.data
data <- AddModuleScore(
  object = data,
  features = gene_list,
  name = "RAMBO"
)

 pdf(paste0(outdir,"/QC2/Umap_Markers.Rambo.signature.pdf"), width=12, height=12)
FeaturePlot_scCustom(data, features = "RAMBO1") +
  ggtitle("RAMBO Module Score")
 dev.off()




 pdf(paste0(outdir,"/QC2/Umap_Markers.Rambo.signature.split.pdf"), width=24, height=12)
FeaturePlot_scCustom(data, features = "RAMBO1", split.by="group")
 dev.off()



 pdf(paste0(outdir,"/QC2/Umap_Markers.Rambo.signature.split.Il1b.pdf"), width=24, height=12)
FeaturePlot_scCustom(data, features = "Il1b", split.by="group")
 dev.off()



### Pseudobulk


data@meta.data$sample <- data@meta.data$tag
table(data@meta.data$sample)

library(Matrix)
library(limma)
library(edgeR)
library(Matrix.utils)
library(dplyr)

# --------------------------------------------
# 1️⃣ Pseudobulk por cluster usando counts crudos
# --------------------------------------------
pseudobulk_cluster <- function(seurat_obj, cluster_name, cluster_col="Cluster", sample_col="tag"){
  
  # counts crudos del assay RNA
  counts <- GetAssayData(seurat_obj, assay="RNA", slot="counts")
  meta   <- seurat_obj@meta.data
  
  # seleccionar células del cluster
  cells <- rownames(meta)[meta[[cluster_col]] == cluster_name]
  counts <- counts[, cells, drop=FALSE]
  
  # vector de samples
  samples <- as.character(meta[cells, sample_col])
  
  # sumar counts por sample
  pb_mat <- t(aggregate.Matrix(t(counts), groupings=samples, fun="sum"))
  
  # número de células por sample
  n_cells <- table(samples)
  n_cells <- n_cells[colnames(pb_mat)]  # asegurar orden
  
  # total counts por sample
  n_counts <- colSums(pb_mat)
  
  list(pb_mat=pb_mat, n_cells=n_cells, n_counts=n_counts)
}

# --------------------------------------------
# 2️⃣ Metadata de samples
# --------------------------------------------
make_sample_metadata <- function(sample_names){
  df <- data.frame(sample = sample_names)
  df$genotype <- ifelse(grepl("KO", sample_names), "KO", "WT")
  df$dsred    <- ifelse(grepl("DsRedP", sample_names), "DsRedP", "DsRedN")
  rownames(df) <- df$sample
  df
}

# --------------------------------------------
# 3️⃣ Limma-voom ajustando por nCells y nCounts
# --------------------------------------------
run_limma_cluster <- function(pb_mat, sample_meta, n_cells){
  
  # crear objeto DGEList
  dge <- DGEList(counts = pb_mat)
  keep <- rowSums(dge$counts) > 10
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  # asegurar factor de genotype
  sample_meta$genotype <- factor(sample_meta$genotype, levels = c("WT", "KO"))
  
  if (nlevels(sample_meta$genotype) < 2) {
    stop("Solo hay un nivel de genotype en este cluster → no se puede hacer contraste KO vs WT")
  }
  
  # añadir covariables: nCells y nCounts
  sample_meta$nCells <- as.numeric(n_cells[colnames(pb_mat)])

  
  # diseño
  design <- model.matrix(~ nCells + genotype, data = sample_meta)
  
  # voom
  v <- voom(dge, design, plot=FALSE)
  
  # ajuste
  fit <- lmFit(v, design)
  
  # contraste KO vs WT
  cont.matrix <- makeContrasts(KOvsWT = genotypeKO, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  # tabla de resultados
  tt <- topTable(fit2, coef="KOvsWT", number=Inf, sort.by="none")
  tt <- tt[order(tt$adj.P.Val), ]
  tt$gene <- rownames(tt)
  tt
}

# --------------------------------------------
# 4️⃣ Ejemplo de uso
# --------------------------------------------

# Ly6c|Ms4ac Monocytes

res <- pseudobulk_cluster(data, cluster_name="Ly6c|Ms4ac Monocytes")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)


p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Ly6c|Ms4ac Monocytes")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Ly6c|Ms4ac Monocytes.pdf"), width=16, height=12)
print(p)
dev.off()



write.table(res_limma, paste0(rsdir,"table.macros.Ly6c|Ms4ac Monocytes.tsv"), sep='\t')


# Trem1|Ptgs2|Plaur|Celc4e


res <- pseudobulk_cluster(data, cluster_name="Trem1|Ptgs2|Plaur|Celc4e Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples.Trem1|Ptgs2|Plaur|Celc4e Mac")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Trem1|Ptgs2|Plaur|Celc4e Mac.pdf"), width=16, height=12)
print(p)
dev.off()

write.table(res_limma, paste0(rsdir,"table.macros.Trem1|Ptgs2|Plaur|Celc4e Mac.tsv"), sep='\t')




# Mrc1|C1QC2|Cbr2|Gas6 Mac


res <- pseudobulk_cluster(data, cluster_name="Mrc1|C1QC2|Cbr2|Gas6 Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples.Mrc1|C1QC2|Cbr2|Gas6 Mac")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Mrc1|C1QC2|Cbr2|Gas6 Mac.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Mrc1|C1QC2|Cbr2|Gas6 Mac.tsv"), sep='\t')






# Arg1|Spp1|Mmp12|Mmp19|Il1a Mac


res <- pseudobulk_cluster(data, cluster_name="Arg1|Spp1|Mmp12|Mmp19|Il1a Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)




p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Arg1|Spp1|Mmp12|Mmp19|Il1a Mac")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Arg1|Spp1|Mmp12|Mmp19|Il1a Mac.pdf"), width=16, height=12)
print(p)
dev.off()

write.table(res_limma, paste0(rsdir,"table.macros.Arg1|Spp1|Mmp12|Mmp19|Il1a Mac.tsv"), sep='\t')




# Npr2|Actn1 Mac


res <- pseudobulk_cluster(data, cluster_name="Npr2|Actn1 Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Npr2|Actn1 Mac")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Npr2|Actn1 Mac.pdf"), width=16, height=12)
print(p)
dev.off()



write.table(res_limma, paste0(rsdir,"table.macros.Npr2|Actn1 Mac.tsv"), sep='\t')




# Mmp9|Ctsk Mac


res <- pseudobulk_cluster(data, cluster_name="Mmp9|Ctsk Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)


p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Mmp9|Ctsk Mac")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Mmp9|Ctsk Mac.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Mmp9|Ctsk Mac.tsv"), sep='\t')




# Ciita|Siglec Mac 


res <- pseudobulk_cluster(data, cluster_name="Ciita|Siglec Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Ciita|Siglec Mac")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Ciita|Siglec Mac.pdf"), width=16, height=12)
print(p)
dev.off()



write.table(res_limma, paste0(rsdir,"table.macros.Ciita|Siglec Mac.tsv"), sep='\t')




 # Ciita|Ccl12 Mac

res <- pseudobulk_cluster(data, cluster_name="Ciita|Ccl12 Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Ciita|Ccl12 Mac")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Ciita|Ccl12 Mac.pdf"), width=16, height=12)
print(p)
dev.off()




write.table(res_limma, paste0(rsdir,"table.macros.Ciita|Ccl12 Mac.tsv"), sep='\t')


# IFN Mac 


res <- pseudobulk_cluster(data, cluster_name="IFN Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. IFN Mac")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.IFN Mac.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.IFN Mac.tsv"), sep='\t')




# Slamf Mac 


res <- pseudobulk_cluster(data, cluster_name="Slamf Monocytes")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Slamf Monocytes")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Slamf Monocytes.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Slamf Monocytes.tsv"), sep='\t')






# Fn1 Mac 


res <- pseudobulk_cluster(data, cluster_name="Fn1|Vegfa Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Fn1|Vegfa Mac")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Fn1|Vegfa Mac.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Fn1|Vegfa Mac.tsv"), sep='\t')




# Neutrophils 


res <- pseudobulk_cluster(data, cluster_name="Neutrophils")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Neutrophils")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.Neutrophils.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Neutrophils.tsv"), sep='\t')







# All macros clusters 

library(Matrix)
library(limma)
library(edgeR)
library(Matrix.utils)
library(dplyr)

# 1️⃣ Pseudobulk agrupando varios clusters
pseudobulk_group <- function(seurat_obj, clusters, cluster_col="Cluster", sample_col="tag"){
  
  counts <- GetAssayData(seurat_obj, assay="RNA", slot="counts")
  meta   <- seurat_obj@meta.data
  
  # células de los clusters de interés
  cells <- rownames(meta)[meta[[cluster_col]] %in% clusters]
  counts <- counts[, cells, drop=FALSE]
  
  # vector de samples
  samples <- as.character(meta[cells, sample_col])
  
  # sumar counts por sample
  pb_mat <- t(aggregate.Matrix(t(counts), groupings=samples, fun="sum"))
  
  # número de células por sample
  n_cells <- table(samples)
  n_cells <- n_cells[colnames(pb_mat)]
  
  # total counts por sample
  n_counts <- colSums(pb_mat)
  
  list(pb_mat=pb_mat, n_cells=n_cells)
}

# 2️⃣ Metadata de samples
make_sample_metadata <- function(sample_names){
  df <- data.frame(sample = sample_names)
  df$genotype <- ifelse(grepl("KO", sample_names), "KO", "WT")
  df$dsred    <- ifelse(grepl("DsRedP", sample_names), "DsRedP", "DsRedN")
  rownames(df) <- df$sample
  df
}

# 3️⃣ Limma-voom
library(limma)



# 3️⃣ Limma-voom
run_limma <- function(pb_mat, sample_meta, n_cells){
  
  dge <- DGEList(counts = pb_mat)
  keep <- rowSums(dge$counts) > 10
  dge <- dge[keep, , keep.lib.sizes=FALSE]
  dge <- calcNormFactors(dge)
  
  sample_meta$genotype <- factor(sample_meta$genotype, levels = c("WT", "KO"))
  sample_meta$nCells <- as.numeric(n_cells[colnames(pb_mat)])
  
  # diseño → corregimos solo por nº de células
  design <- model.matrix(~ nCells + genotype, data = sample_meta)
  
  v <- voom(dge, design, plot=FALSE)
  fit <- lmFit(v, design)
  
  cont.matrix <- makeContrasts(KOvsWT = genotypeKO, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)
  
  tt <- topTable(fit2, coef="KOvsWT", number=Inf, sort.by="none")
  tt <- tt[order(tt$adj.P.Val), ]
  tt$gene <- rownames(tt)
  tt
}


# ---------------------------
# 🚀 Ejemplo: todos los clusters de macrófagos
# ---------------------------


# Lista de clusters
macro_clusters <- c(
  "Ly6c|Ms4ac Monocytes",
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac",
  "Ciita|Ccl12 Mac",
  "Ciita|Siglec Mac",
  "Fn1|Vegfa Mac",
  "IFN Mac",
  "Mmp9|Ctsk Mac",
  "Mrc1|C1QC2|Cbr2|Gas6 Mac",
  "Neutrophils",
  "Npr2|Actn1 Mac",
  "Slamf Monocytes",
  "Trem1|Ptgs2|Plaur|Celc4e Mac"
)

pb_out <- pseudobulk_group(data, clusters = macro_clusters, cluster_col="Cluster", sample_col="tag")
sample_meta <- make_sample_metadata(colnames(pb_out$pb_mat))

res_allmacs <- run_limma(pb_out$pb_mat, sample_meta, pb_out$n_cells)

head(res_allmacs)


res_allmacs$diffexpressed <- "NO"
res_allmacs$diffexpressed[res_allmacs$adj.P.Val < 0.05 & res_allmacs$logFC > 1] <- "Up"
res_allmacs$diffexpressed[res_allmacs$adj.P.Val < 0.05 & res_allmacs$logFC < -1] <- "Down"



# Top 20 genes
top_genes <- res_allmacs %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_allmacs$diffexpressed)


write.table(res_allmacs, paste0(rsdir,"table.macros.ALL Mac.tsv"), sep='\t')


Volcano2 <- function(data, legend, logFC_threshold = 0.5) {
  data$gene <- rownames(data)
  
  # Solo etiquetar genes diferencialmente expresados
  data$label <- ifelse(data$diffexpressed != "NO", data$gene, NA)
  
  plot <- ggplot(data, aes(x = logFC, y = -log10(adj.P.Val),
                           color = diffexpressed)) +
    geom_point() +
    geom_text_repel(
      data = subset(data, !is.na(label)),
      aes(label = label),
      size = 3,
      segment.colour = NA,
      force = 2,
      box.padding = 0.5,
      point.padding = 0.5
    ) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold),
               linetype = "dashed", color = "black") +
    theme_minimal() +
    scale_color_manual(values = c("Down" = "dodgerblue4",
                                  "NO"   = "dimgrey",
                                  "Up"   = "brown1")) +
    ggtitle(paste("DEG", legend)) +
    labs(x = "logFC", y = "-log10(adj.pval)")
  
  return(plot)
}




p <- Volcano2(data = res_allmacs, legend = "KO vs WT samples. All Macrophages clusters")


pdf(paste0(outdir,"/QC2/Volcano.KO vs WT Macrophages.pdf"), width=16, height=12)
print(p)
dev.off()






### Man plot for all 
cluster.Ly6c.Ms4ac <- read.table(paste0(rsdir,"table.macros.Ly6c|Ms4ac Monocytes.tsv"), sep='\t', header=T)
cluster.Trem1.Ptgs2.Plaur.Celc4e <- read.table(paste0(rsdir,"table.macros.Trem1|Ptgs2|Plaur|Celc4e Mac.tsv"), sep='\t', header=T)
cluster.Mrc1.C1QC2.Cbr2.Gas6 <- read.table(paste0(rsdir,"table.macros.Mrc1|C1QC2|Cbr2|Gas6 Mac.tsv"), sep='\t', header=T)
cluster.Arg1.Spp1.Mmp12.Mmp19.Il1a <- read.table(paste0(rsdir,"table.macros.Arg1|Spp1|Mmp12|Mmp19|Il1a Mac.tsv"), sep='\t', header=T)
cluster.Npr2.Actn1 <- read.table(paste0(rsdir,"table.macros.Npr2|Actn1 Mac.tsv"), sep='\t', header=T)
cluster.Slamf <- read.table(paste0(rsdir,"table.macros.Slamf Monocytes.tsv"), sep='\t', header=T)
cluster.Mmp9.Ctsk <- read.table(paste0(rsdir,"table.macros.Mmp9|Ctsk Mac.tsv"), sep='\t', header=T)
cluster.IFN <- read.table(paste0(rsdir,"table.macros.IFN Mac.tsv"), sep='\t', header=T)
cluster.Fn1.Vegfa <- read.table(paste0(rsdir,"table.macros.Fn1|Vegfa Mac.tsv"), sep='\t', header=T)
cluster.Ciita.Siglec <- read.table(paste0(rsdir,"table.macros.Ciita|Siglec Mac.tsv"), sep='\t', header=T)
cluster.Ciita.Ccl12 <- read.table(paste0(rsdir,"table.macros.Ciita|Ccl12 Mac.tsv"), sep='\t', header=T)
cluster.Neutrophils <- read.table(paste0(rsdir,"table.macros.Neutrophils.tsv"), sep='\t', header=T)




cluster.Ly6c.Ms4ac$cluster <- "Ly6c_Ms4ac"
cluster.Trem1.Ptgs2.Plaur.Celc4e$cluster <- "Trem1_Ptgs2_Plaur_Celc4e"
cluster.Mrc1.C1QC2.Cbr2.Gas6$cluster <- "Mrc1_C1QC2_Cbr2_Gas6"
cluster.Arg1.Spp1.Mmp12.Mmp19.Il1a$cluster <- "Arg1_Spp1_Mmp12_Mmp19_Il1a"
cluster.Npr2.Actn1$cluster <- "Npr2_Actn1"
cluster.Slamf$cluster <- "Slamf"
cluster.Mmp9.Ctsk$cluster <- "Mmp9_Ctsk"
cluster.IFN$cluster <- "IFN"
cluster.Fn1.Vegfa$cluster <- "Fn1_Vegfa"
cluster.Ciita.Siglec$cluster <- "Ciita_Siglec"
cluster.Ciita.Ccl12$cluster <- "Ciita_Ccl12"
cluster.Neutrophils$cluster <- "Neutrophils"

# Unir todo
all_clusters <- bind_rows(
  cluster.Ly6c.Ms4ac,
  cluster.Trem1.Ptgs2.Plaur.Celc4e,
  cluster.Mrc1.C1QC2.Cbr2.Gas6,
  cluster.Arg1.Spp1.Mmp12.Mmp19.Il1a,
  cluster.Npr2.Actn1,
  cluster.Slamf,
  cluster.Mmp9.Ctsk,
  cluster.IFN,
  cluster.Fn1.Vegfa,
  cluster.Ciita.Siglec,
  cluster.Ciita.Ccl12,
  cluster.Neutrophils
)





nora.colors2 <- c(
  "Ly6c|Ms4ac Monocytes"           = "#FF3B30",
  "Trem1|Ptgs2|Plaur|Celc4e Mac"   = "#EE7942",
  "Mrc1|C1QC2|Cbr2|Gas6 Mac"        = "#FFD92F",
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac" = "#4DAF4A",
  "Npr2|Actn1 Mac"                 = "#A6D854",
  "Neutrophils"                    = "#4876FF",
  "Slamf Monocytes"                = "#66B2FF",  # corregido
  "Mmp9|Ctsk Mac"                  = "#00723F",
  "IFN Mac"                        = "#C080FF",
  "Fn1|Vegfa Mac"                  = "#FFA500",
  "Ciita|Siglec Mac"               = "#1E90FF",
  "Ciita|Ccl12 Mac"                = "#4682B4"
)




# Reemplazar "_" por "|" y mapear a nora.colors
nora_names <- names(nora.colors2)
nora_clusters <- sapply(strsplit(nora_names, " "), `[`, 1)  # parte del cluster sin tipo de célula
cluster_map <- setNames(nora_names, gsub("\\|", "_", nora_clusters))

# Reemplazar los nombres en tu df
all_clusters$cluster <- cluster_map[all_clusters$cluster]

# Verificar que todos los clusters tengan mapeo
all_clusters$cluster[is.na(all_clusters$cluster)]  # Estos serían los que no se mapean



deg_shared <- all_clusters %>%
  filter(diffexpressed %in% c("Up", "Down")) %>%
  group_by(gene, diffexpressed) %>%
  summarise(n_clusters = n_distinct(cluster), .groups = "drop") %>%
  filter(n_clusters >= 3) %>%  # ≥3 clusters
  arrange(desc(n_clusters))

as.data.frame(deg_shared)



library(ggplot2)
library(dplyr)
library(ggrepel)

# Vector de genes a etiquetar
genes_to_label <- c(
  "Cxcl1","C1QC2","Hexim1","Lgals1","Marcksl1","Nfkbia","Oas2","Ahnak","Cxcl9","Egr2",
  "Ifit3","Ifit3b","Ifitm3","Il1a","Irf7","Ms4a4c","Ly6a","Stat2","Atp6v0a1","rp2",
  "C3","Cxcl3","Per2","Atp5a1","Socs3","Abca9","Ehd1","Oas1g","H2-T24","Rtp4","Cxcl2",
  "Timeless"
)

# Filtrar solo genes DE (Up/Down) que queremos etiquetar
label_points <- all_clusters %>%
  filter(gene %in% genes_to_label & diffexpressed %in% c("Up","Down"))

# Suavizar colores rojo y azul
color_up <- "#FF6666"   # tono más claro que rojo
color_down <- "#6699FF" # tono más claro que azul


cluster_levels <- levels(all_clusters$cluster)


cluster_boxes <- data.frame(
  cluster = cluster_levels,
  xmin = seq_along(cluster_levels) - 0.4,  # mantener ancho actual
  xmax = seq_along(cluster_levels) + 0.4,
  ymin = -0.3,   # más alto
  ymax = 0.3     # más alto
)
cluster_boxes$cluster <- factor(cluster_boxes$cluster, levels = cluster_levels)


# Plot
pdf(paste0(outdir, "/QC2/Shared_genes_metric_plot_rects_final_labels.pdf"), width = 22, height = 14)

ggplot(all_clusters, aes(x = cluster, y = metric)) +
  # recuadros centrados en 0
  geom_rect(data = cluster_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
            inherit.aes = FALSE, alpha = 1) +
  # todos los genes
  geom_jitter(width = 0.25, alpha = 0.3, color = "grey") +
  # genes Up/Down con color suavizado
  geom_point(data = filter(all_clusters, diffexpressed %in% c("Up","Down")),
             aes(x = cluster, y = metric, color = diffexpressed),
             position = position_jitter(width = 0.25),
             size = 2) +
  scale_color_manual(values = c("Up" = color_up, "Down" = color_down)) +
  scale_fill_manual(name = "Cluster", values = nora.colors2) +
  # etiquetas con ggrepel
  geom_text_repel(data = label_points,
                  aes(label = gene, x = cluster, y = metric),
                  size = 5,
                  box.padding = 0.4,
                  point.padding = 0.3,
                  segment.size = 0.3,
                  max.overlaps = Inf) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = "logFC * -log10(adj.P.Val)",
    title = "KO vs WT samples by cluster. Differentially expressed genes"
  )

dev.off()




library(dplyr)
library(tidyr)
library(pheatmap)

# --- 1. Top 80 genes ---
deg_shared_df <- as.data.frame(deg_shared)
deg_shared_ordered <- deg_shared_df[order(-deg_shared_df$n_clusters), ]
top_genes <- deg_shared_ordered$gene[1:80]

# --- 2. Preparar matriz genes x clusters numérica ---
heat_wide <- all_clusters %>%
  filter(gene %in% top_genes) %>%
  mutate(value = case_when(
    diffexpressed == "Up"   ~ 1,
    diffexpressed == "Down" ~ -1,
    TRUE                    ~ 0
  )) %>%
  select(gene, cluster, value) %>%
  pivot_wider(names_from = cluster, values_from = value, values_fill = 0)

heat_mat <- as.matrix(heat_wide[,-1])
rownames(heat_mat) <- heat_wide$gene

# --- 3. Generar PDF con heatmap ---
pdf(paste0(outdir, "/QC2/Top80_genes_heatmap_discrete_numeric.pdf"), width = 16, height = 20)

pheatmap(
  mat = heat_mat,
  color = colorRampPalette(c("#6699FF","white","#FF6666"))(3),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 8,
  fontsize_col = 10,
  main = "Top 80 shared genes",
  legend_breaks = c(-1, 0, 1),
  legend_labels = c("Down regulated","Not DE","Up regulated")
)

dev.off()













differential_allgenes_by_group_in_cluster <- function(
  seurat_obj,
  cluster_name,
  cluster_col = "Cluster",
  group_col = "group",
  ident1,
  ident2
) {
  # Filtrar solo las células del cluster deseado
  cluster_cells <- subset(seurat_obj, subset = !!as.name(cluster_col) == cluster_name)
  
  # Establecer identidad de agrupamiento por la columna group_col
  Idents(cluster_cells) <- group_col
  
  # Ejecutar FindMarkers sin umbrales para traer TODOS los genes
  markers <- FindMarkers(
    object = cluster_cells,
    ident.1 = ident1,
    ident.2 = ident2,
    logfc.threshold = 0,   # sin filtro de logFC
    min.pct = 0,           # sin filtro de pct
    return.thresh = 1.01,  # devuelve todos los genes
    only.pos = FALSE,       # también los genes down
    test.use = "MAST",
    slot= "data",
    latent.vars= "percent.mt"  # método estadístico
  )
  
  # Añadir meta-info
  markers$Cluster <- cluster_name
  markers$Comparison <- paste(ident1, "vs", ident2)
  markers$Gene <- rownames(markers)
  rownames(markers) <- NULL
  
  return(markers)
}




markers_cluster <- differential_allgenes_by_group_in_cluster(
  seurat_obj = data,
  cluster_name = "Neutrophils",
  cluster_col = "Cluster",
  group_col = "group",
  ident1 = "KO",
  ident2 = "WT"

)

head(markers_cluster, n=50)

write.table(markers_cluster, paste0(rsdir,"table.macros.Arg1|Mrc1 Mac.tsv"), sep='\t')



markers_cluster <- differential_allgenes_by_group_in_cluster(
  seurat_obj = data,
  cluster_name = "Arg1+ Mac",
  cluster_col = "Cluster",
  group_col = "group",
  ident1 = "KO",
  ident2 = "WT"

)

head(markers_cluster)

write.table(markers_cluster, paste0(rsdir,"table.macros.Arg1+ Mac.tsv"), sep='\t')



markers_cluster <- differential_allgenes_by_group_in_cluster(
  seurat_obj = data,
  cluster_name = "Mrc1+ Mac",
  cluster_col = "Cluster",
  group_col = "group",
  ident1 = "KO",
  ident2 = "WT"

)

head(markers_cluster)

write.table(markers_cluster, paste0(rsdir,"table.macros.Mrc1+ Mac.tsv"), sep='\t')



markers_cluster <- differential_allgenes_by_group_in_cluster(
  seurat_obj = data,
  cluster_name = "Mrc1+ Mac",
  cluster_col = "Cluster",
  group_col = "group",
  ident1 = "KO",
  ident2 = "WT"

)

head(markers_cluster)

write.table(markers_cluster, paste0(rsdir,"table.macros.Mrc1+ Mac.tsv"), sep='\t')




markers_cluster <- differential_allgenes_by_group_in_cluster(
  seurat_obj = data,
  cluster_name = "Mki67+ Mac",
  cluster_col = "Cluster",
  group_col = "group",
  ident1 = "KO",
  ident2 = "WT"

)

head(markers_cluster)

write.table(markers_cluster, paste0(rsdir,"table.macros.Mki67+ Mac.tsv"), sep='\t')




markers_cluster <- differential_allgenes_by_group_in_cluster(
  seurat_obj = data,
  cluster_name = "Monocytes Ly6c2Hi",
  cluster_col = "Cluster",
  group_col = "group",
  ident1 = "KO",
  ident2 = "WT"

)

head(markers_cluster)

write.table(markers_cluster, paste0(rsdir,"table.macros.Monocytes Ly6c2Hi.tsv"), sep='\t')




markers_cluster <- differential_allgenes_by_group_in_cluster(
  seurat_obj = data,
  cluster_name = "Gpnmb+ Mac",
  cluster_col = "Cluster",
  group_col = "group",
  ident1 = "KO",
  ident2 = "WT"

)

head(markers_cluster)

write.table(markers_cluster, paste0(rsdir,"table.macros.Gpnmb+ Mac.tsv"), sep='\t')




markers_cluster <- differential_allgenes_by_group_in_cluster(
  seurat_obj = data,
  cluster_name = "IFN Mac",
  cluster_col = "Cluster",
  group_col = "group",
  ident1 = "KO",
  ident2 = "WT"

)

head(markers_cluster)

write.table(markers_cluster, paste0(rsdir,"table.macros.Ciita|Mrc1|MHC-II Mac.tsv"), sep='\t')









# Highlight KO3-DsRedp cells in the UMAP plot
# Asegúrate de que la columna 'tag' esté presente en los metadatos

# Añade una nueva columna lógica para marcar las células KO3-DsRedp
data$highlight_KO3 <- ifelse(data$tag == "KO3-DsRedp", "KO3-DsRedp", "Other")

# Especificar los niveles del factor para controlar el orden y el color
data$highlight_KO3 <- factor(data$highlight_KO3, levels = c("Other", "KO3-DsRedp"))


# UMAP coloreado por esa nueva columna
png(paste0(outdir,"/QC2/Umap_Markers2.png"), width=1800, height=2400)
DimPlot(data, group.by = "highlight_KO3", cols = c("gray", "red"), split.by="group") +
  ggtitle("KO3-DsRedp Cells")
dev.off()




library(Seurat)
library(dplyr)

calcular_resumen <- function(obj, assay_name) {
  meta <- obj@meta.data %>% select(tag, percent.mt,
                                   nFeature_RNA, nCount_RNA,
                                   nFeature_SCT, nCount_SCT)
  
  if (assay_name == "RNA") {
    resumen <- meta %>%
      group_by(tag) %>%
      summarise(
        mean_genes = mean(nFeature_RNA),
        mean_counts = mean(nCount_RNA),
        mean_percent_mt = mean(percent.mt) 
      ) %>%
      mutate(assay = assay_name)
  } else if (assay_name == "SCT") {
    resumen <- meta %>%
      group_by(tag) %>%
      summarise(
        mean_genes = mean(nFeature_SCT),
        mean_counts = mean(nCount_SCT),
        mean_percent_mt = mean(percent.mt) 
      ) %>%
      mutate(assay = assay_name)
  } else {
    stop("Assay no reconocido")
  }
  
  return(resumen)
}

resumen_RNA <- calcular_resumen(data, "RNA")
resumen_SCT <- calcular_resumen(data, "SCT")

resumen_total <- bind_rows(resumen_RNA, resumen_SCT)

print(resumen_total)
write.csv(resumen_total, "resumen_por_muestra.csv", row.names = FALSE)


