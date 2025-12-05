# Seurat Enviroment

library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(celldex)
library(SingleR)
library(SingleCellExperiment)
library(SoupX)
library(future)
library(dplyr)
library(DoubletFinder)


# Reading the config file

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

options(future.globals.maxSize = 16*1024^3)  # 2 GB
plan(sequential)


 # SoupX


process_library(paste0(datadir, "/Cellranger/smR065"), output_dir = paste0(datadir, "/SoupX/2.0/smR065"))
process_library(paste0(datadir, "/Cellranger/SmR039a.9"), output_dir = paste0(datadir, "/SoupX/2.0/SmR039a.9"))
process_library(paste0(datadir, "/Cellranger/SmR039b.9"), output_dir = paste0(datadir, "/SoupX/2.0/SmR039b.9"))


#### Doublets


smpR065_DsRedP_KO3 <-readRDS(paste0(datadir,"SoupX/2.0/smR065/Seurat_smpR065_DsRedP-KO3.rds"))

smR039a_WT1_DsRedn <-readRDS(paste0(datadir,"SoupX/2.0/SmR039a.9/Seurat_WT1-DsRedn.rds"))
smR039a_WT1_DsRedp <-readRDS(paste0(datadir,"SoupX/2.0/SmR039a.9/Seurat_WT1-DsRedp.rds"))
smR039a_DsRedN_KO2 <-readRDS(paste0(datadir,"SoupX/2.0/SmR039a.9/Seurat_KO2-DsRedn.rds"))
smR039a_DsRedP_KO2 <-readRDS(paste0(datadir,"SoupX/2.0/SmR039a.9/Seurat_KO2-DsRedp.rds"))

smR039b_WT2_DsRedn <-readRDS(paste0(datadir,"SoupX/2.0/SmR039b.9/Seurat_WT2-DsRedn.rds"))
smR039b_WT2_DsRedp <-readRDS(paste0(datadir,"SoupX/2.0/SmR039b.9/Seurat_WT2-DsRedp.rds"))
smR039b_DsRedN_KO3 <-readRDS(paste0(datadir,"SoupX/2.0/SmR039b.9/Seurat_KO3-DsRedn.rds"))


seurat_list <- list(

  smpR065_DsRedP_KO3,
  smR039a_WT1_DsRedn,
  smR039a_WT1_DsRedp,
  smR039a_DsRedN_KO2,
  smR039a_DsRedP_KO2,
  smR039b_WT2_DsRedn,
  smR039b_WT2_DsRedp,
  smR039b_DsRedN_KO3
)

seurat_names <- c(
  "smpR065_DsRedP_KO3",
  "smR039a_WT1_DsRedn",
  "smR039a_WT1_DsRedp",
  "smR039a_DsRedN_KO2",
  "smR039a_DsRedP_KO2",
  "smR039b_WT2_DsRedn",
  "smR039b_WT2_DsRedp",
  "smR039b_DsRedN_KO3"
)

process_seurat_samples(seurat_list, seurat_names)


##### PCA


PCA.DsRedP_KO3<- PCA(smpR065_DsRedP_KO3)
PCA.DsRedP_KO3[1]
PCA.DsRedP_KO3[2]
#12

PCA.WT1_DsRedn <- PCA(smR039a_WT1_DsRedn)
PCA.WT1_DsRedn[1]
PCA.WT1_DsRedn[2]
#15

PCA.WT1_DsRedp <- PCA(smR039a_WT1_DsRedp)
PCA.WT1_DsRedp[1]
PCA.WT1_DsRedp[2]
#17

PCA.DsRedN_KO2 <- PCA(smR039a_DsRedN_KO2)
PCA.DsRedN_KO2[1]
PCA.DsRedN_KO2[2]
#10

PCA.DsRedP_KO2 <- PCA(smR039a_DsRedP_KO2)
PCA.DsRedP_KO2[1]
PCA.DsRedP_KO2[2]
#14

PCA.WT2_DsRedn <- PCA(smR039b_WT2_DsRedn)
PCA.WT2_DsRedn[1]
PCA.WT2_DsRedn[2]
#14

PCA.WT2_DsRedp <- PCA(smR039b_WT2_DsRedp)
PCA.WT2_DsRedp[1]
PCA.WT2_DsRedp[2]
#12

PCA.DsRedN_KO3<- PCA(smR039b_DsRedN_KO3)
PCA.DsRedN_KO3[1]
PCA.DsRedN_KO3[2]
#11



#UMAP


smpR065_DsRedP_KO3 <- umap(smpR065_DsRedP_KO3,"12")

smR039a_WT1_DsRedn <- umap(smR039a_WT1_DsRedn,"15")
smR039a_WT1_DsRedp <- umap(smR039a_WT1_DsRedp,"17")
smR039a_DsRedN_KO2 <- umap(smR039a_DsRedN_KO2,"10")
smR039a_DsRedP_KO2 <- umap(smR039a_DsRedP_KO2,"14")


smR039b_WT2_DsRedn <- umap(smR039b_WT2_DsRedn,"14")
smR039b_WT2_DsRedp <- umap(smR039b_WT2_DsRedp,"12")
smR039b_DsRedN_KO3 <- umap(smR039b_DsRedN_KO3,"11")




## Doublets



process_doublets_single_with_rate <- function(seurat_obj, nPCs, pN = 0.25, doublet_rate_base = 0.075) {


# Normalización y preparación si falta
if (!"data" %in% Layers(seurat_obj)) {
  seurat_obj <- NormalizeData(seurat_obj)
}
if (length(VariableFeatures(seurat_obj)) == 0) {
  seurat_obj <- FindVariableFeatures(seurat_obj)
}
if (!"scale.data" %in% Layers(seurat_obj)) {
  seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj))
}
if (!"pca" %in% names(seurat_obj@reductions)) {
  seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
}
if (!"seurat_clusters" %in% colnames(seurat_obj@meta.data)) {
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:nPCs)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
}

clusters <- seurat_obj$seurat_clusters

# ParamSweep y summarize
sweep.res.list <- paramSweep(seurat_obj, PCs = 1:nPCs, sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list, GT.calls = FALSE)

sweep.stats$pN     <- as.numeric(as.character(sweep.stats$pN))
sweep.stats$pK     <- as.numeric(as.character(sweep.stats$pK))
sweep.stats$BCreal <- as.numeric(as.character(sweep.stats$BCreal))

sweep.stats <- data.frame(
  pN = sweep.stats$pN,
  pK = sweep.stats$pK,
  BCreal = sweep.stats$BCreal,
  stringsAsFactors = FALSE
)



# Asegurar formato numérico
sweep.stats$pK <- as.numeric(as.character(sweep.stats$pK))
sweep.stats$BCreal <- as.numeric(as.character(sweep.stats$BCreal))

# Extraer pK óptimo con máximo BCreal
pk.max <- sweep.stats$pK[which.max(sweep.stats$BCreal)]
pk.max


cat("pK óptimo:", pk.max, "\n")

# Número esperado de dobletes
n_cells <- ncol(seurat_obj)
estimated_rate <- doublet_rate_base * (n_cells / 10000)  # estimación proporcional
homotypic.prop <- modelHomotypic(clusters)
nExp <- round(estimated_rate * n_cells * (1 - homotypic.prop))

cat("Número de células:", n_cells, "\n")
cat("Doublet rate estimado para esta muestra:", round(estimated_rate * 100, 2), "%\n")
cat("Número esperado de dobletes ajustado:", nExp, "\n")

seurat_obj$seurat_clusters <- as.factor(seurat_obj$seurat_clusters)

pk.max <- as.numeric(pk.max)[1]
pN     <- as.numeric(pN)[1]
nExp   <- as.integer(nExp)[1]


# Ejecutar DoubletFinder
seurat_obj <- doubletFinder(
  seurat_obj,
  PCs = 1:nPCs,
  pN = pN,
  pK = pk.max,
  nExp = nExp,
  reuse.pANN = NULL,
  sct = FALSE
)

  
  return(seurat_obj)
}





smpR065_DsRedP_KO3 <- process_doublets_single_with_rate(smpR065_DsRedP_KO3, nPCs = 12, pN = 0.25, doublet_rate_base = 0.075)

smR039a_WT1_DsRedn <- process_doublets_single_with_rate(smR039a_WT1_DsRedn, nPCs = 15, pN = 0.25, doublet_rate_base = 0.075)
smR039a_WT1_DsRedp <- process_doublets_single_with_rate(smR039a_WT1_DsRedp, nPCs = 17, pN = 0.25, doublet_rate_base = 0.075)
smR039a_DsRedN_KO2 <- process_doublets_single_with_rate(smR039a_DsRedN_KO2, nPCs = 10, pN = 0.25, doublet_rate_base = 0.075)
smR039a_DsRedP_KO2 <- process_doublets_single_with_rate(smR039a_DsRedP_KO2, nPCs = 14, pN = 0.25, doublet_rate_base = 0.075)

smR039b_WT2_DsRedn <- process_doublets_single_with_rate(smR039b_WT2_DsRedn, nPCs = 14, pN = 0.25, doublet_rate_base = 0.075)
smR039b_WT2_DsRedp <- process_doublets_single_with_rate(smR039b_WT2_DsRedp, nPCs = 12, pN = 0.25, doublet_rate_base = 0.075)
smR039b_DsRedN_KO3 <- process_doublets_single_with_rate(smR039b_DsRedN_KO3, nPCs = 11, pN = 0.25, doublet_rate_base = 0.075)





saveRDS(smpR065_DsRedP_KO3,paste0(datadir,"SoupX/2.0/smR065/Seurat_smpR065_DsRedP-KO3.2.rds"))

saveRDS(smR039a_WT1_DsRedn,paste0(datadir,"SoupX/2.0/SmR039a.9/Seurat_WT1-DsRedn.2.rds"))
saveRDS(smR039a_WT1_DsRedp,paste0(datadir,"SoupX/2.0/SmR039a.9/Seurat_WT1-DsRedp.2.rds"))
saveRDS(smR039a_DsRedN_KO2,paste0(datadir,"SoupX/2.0/SmR039a.9/Seurat_KO2-DsRedn.2.rds"))
saveRDS(smR039a_DsRedP_KO2,paste0(datadir,"SoupX/2.0/SmR039a.9/Seurat_KO2-DsRedp.2.rds"))


saveRDS(smR039b_WT2_DsRedn,paste0(datadir,"SoupX/2.0/SmR039b.9/Seurat_WT2-DsRedn.2.rds"))
saveRDS(smR039b_WT2_DsRedp,paste0(datadir,"SoupX/2.0/SmR039b.9/Seurat_WT2-DsRedp.2.rds"))
saveRDS(smR039b_DsRedN_KO3,paste0(datadir,"SoupX/2.0/SmR039b.9/Seurat_KO3-DsRedn.2.rds"))


# Plots



psmpR065_DsRedP_KO3 <- plot(smpR065_DsRedP_KO3,"DsRedP_KO3")
png(paste0(outdir,"/SoupX.Doublet/2.0/DsRedP_KO3.png"), width=800, height=600)
psmpR065_DsRedP_KO3
dev.off()


psmR039a_WT1_DsRedn <- plot(smR039a_WT1_DsRedn,"WT1_DsRedn")
png(paste0(outdir,"/SoupX.Doublet/2.0/WT1_DsRedn.png"), width=800, height=600)
psmR039a_WT1_DsRedn
dev.off()

psmR039a_WT1_DsRedp <- plot(smR039a_WT1_DsRedp,"WT1_DsRedp")
png(paste0(outdir,"/SoupX.Doublet/2.0/WT1_DsRedp.png"), width=800, height=600)
psmR039a_WT1_DsRedp
dev.off()

psmR039a_DsRedN_KO2 <- plot(smR039a_DsRedN_KO2,"DsRedN_KO2")
png(paste0(outdir,"/SoupX.Doublet/2.0/DsRedN_KO2.png"), width=800, height=600)
psmR039a_DsRedN_KO2
dev.off()

psmR039a_DsRedP_KO2 <- plot(smR039a_DsRedP_KO2,"DsRedP_KO2")
png(paste0(outdir,"/SoupX.Doublet/2.0/DsRedP_KO2.png"), width=800, height=600)
psmR039a_DsRedP_KO2
dev.off()


psmR039b_WT2_DsRedn <- plot(smR039b_WT2_DsRedn,"WT2_DsRedn")
png(paste0(outdir,"/SoupX.Doublet/WT2_DsRedn.png"), width=800, height=600)
psmR039b_WT2_DsRedn
dev.off()


psmR039b_WT2_DsRedp <- plot(smR039b_WT2_DsRedp,"WT2_DsRedp")
png(paste0(outdir,"/SoupX.Doublet/WT2_DsRedp.png"), width=800, height=600)
psmR039b_WT2_DsRedp
dev.off()

psmR039b_DsRedN_KO3 <- plot(smR039b_DsRedN_KO3,"DsRedN_KO3")
png(paste0(outdir,"/SoupX.Doublet/DsRedN_KO3.png"), width=800, height=600)
psmR039b_DsRedN_KO3
dev.off()



png(paste0(outdir,"/SoupX.Doublet/2.0/doublets.png"), width=2200, height=1200)
plot_grid(psmpR065_DsRedP_KO3,
psmR039a_WT1_DsRedn,
psmR039a_WT1_DsRedp,
psmR039a_DsRedN_KO2,
psmR039a_DsRedP_KO2,
psmR039b_WT2_DsRedn,
psmR039b_WT2_DsRedp,
psmR039b_DsRedN_KO3, ncol=4)
dev.off()