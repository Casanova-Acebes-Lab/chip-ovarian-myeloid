# Seurat Enviroment

library(Seurat)
library(dplyr)
library(cowplot)
library(scCustomize)
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


process_library(paste0(datadir, "/Cellranger/smR065"), output_dir = paste0(datadir, "/SoupX/smR065"))
process_library(paste0(datadir, "/Cellranger/SmR039a.9"), output_dir = paste0(datadir, "/SoupX/SmR039a.9"))
process_library(paste0(datadir, "/Cellranger/SmR039b.9"), output_dir = paste0(datadir, "/SoupX/SmR039b.9"))



#### Doublets


smpR065_DsRedN_KO2 <-readRDS(paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedN-KO2.rds"))
smpR065_DsRedN_KO3 <-readRDS(paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedN-KO3.rds"))
smpR065_DsRedP_KO2 <-readRDS(paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedP-KO2.rds"))
smpR065_DsRedP_KO3 <-readRDS(paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedP-KO3.rds"))

smR039a_WT1_DsRedn <-readRDS(paste0(datadir,"SoupX/SmR039a.9/Seurat_WT1-DsRedn.rds"))
smR039a_WT1_DsRedp <-readRDS(paste0(datadir,"SoupX/SmR039a.9/Seurat_WT1-DsRedp.rds"))

smR039b_WT2_DsRedn <-readRDS(paste0(datadir,"SoupX/SmR039b.9/Seurat_WT2-DsRedn.rds"))
smR039b_WT2_DsRedp <-readRDS(paste0(datadir,"SoupX/SmR039b.9/Seurat_WT2-DsRedp.rds"))


seurat_list <- list(
  smpR065_DsRedN_KO2,
  smpR065_DsRedN_KO3,
  smpR065_DsRedP_KO2,
  smpR065_DsRedP_KO3,
  smR039a_WT1_DsRedn,
  smR039a_WT1_DsRedp,
  smR039b_WT2_DsRedn,
  smR039b_WT2_DsRedp
)

seurat_names <- c(
  "smpR065_DsRedN_KO2",
  "smpR065_DsRedN_KO3",
  "smpR065_DsRedP_KO2",
  "smpR065_DsRedP_KO3",
  "smR039a_WT1_DsRedn",
  "smR039a_WT1_DsRedp",
  "smR039b_WT2_DsRedn",
  "smR039b_WT2_DsRedp"
)

process_seurat_samples(seurat_list, seurat_names)


##### PCA


PCA.DsRedN_KO2 <- PCA(smpR065_DsRedN_KO2)
PCA.DsRedN_KO2[1]
PCA.DsRedN_KO2[2]
#14

PCA.DsRedN_KO3 <- PCA(smpR065_DsRedN_KO3)
PCA.DsRedN_KO3[1]
PCA.DsRedN_KO3[2]
#19

PCA.DsRedP_KO2 <- PCA(smpR065_DsRedP_KO2)
PCA.DsRedP_KO2[1]
PCA.DsRedP_KO2[2]
#16

PCA.DsRedP_KO3 <- PCA(smpR065_DsRedP_KO3)
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

PCA.WT2_DsRedn <- PCA(smR039b_WT2_DsRedn)
PCA.WT2_DsRedn[1]
PCA.WT2_DsRedn[2]
#14

PCA.WT2_DsRedp <- PCA(smR039b_WT2_DsRedp)
PCA.WT2_DsRedp[1]
PCA.WT2_DsRedp[2]
#11



#UMAP

smpR065_DsRedN_KO2 <- umap(smpR065_DsRedN_KO2,"14")
smpR065_DsRedN_KO3 <- umap(smpR065_DsRedN_KO3,"19")
smpR065_DsRedP_KO2 <- umap(smpR065_DsRedP_KO2,"16")
smpR065_DsRedP_KO3 <- umap(smpR065_DsRedP_KO3,"12")

smR039a_WT1_DsRedn <- umap(smR039a_WT1_DsRedn,"15")
smR039a_WT1_DsRedp <- umap(smR039a_WT1_DsRedp,"17")

smR039b_WT2_DsRedn <- umap(smR039b_WT2_DsRedn,"14")
smR039b_WT2_DsRedp <- umap(smR039b_WT2_DsRedp,"11")





## Doublets



smpR065_DsRedN_KO2 <- pk(smpR065_DsRedN_KO2,"19")

pm1 <- plot(m1,"M1")
png(paste0(outdir,"/pk/m1.png"), width=800, height=600)
pm1
dev.off()

m2 <- pk(m2,"16")

pm2 <- plot(m2,"M2")
png(paste0(outdir,"/pk/m2.png"), width=800, height=600)
pm2
dev.off()

m3 <- pk(m3,"12")

pm3 <- plot(m3,"M3")
png(paste0(outdir,"/pk/m3.png"), width=800, height=600)
pm3
dev.off()

m4 <- pk(m4,"18")
pm4 <- plot(m4,"M4")
png(paste0(outdir,"/pk/m4.png"), width=800, height=600)
pm4
dev.off()


m10 <- pk(m10,"18")
pm10 <- plot(m10,"M10")
png(paste0(outdir,"/pk/m10.png"), width=800, height=600)
pm10
dev.off()

m12 <- pk(m12,"16")
pm12 <- plot(m12,"M12")
png(paste0(outdir,"/pk/m12.png"), width=800, height=600)
pm12
dev.off()

m13 <- pk(m13,"20")
pm13 <- plot(m13,"M13")
png(paste0(outdir,"/pk/m13.png"), width=800, height=600)
pm13
dev.off()

m14 <- pk(m14,"15")
pm14 <- plot(m14,"M14")
png(paste0(outdir,"/pk/m14.png"), width=800, height=600)
pm14
dev.off()

m15 <- pk(m15,"14")
pm15 <- plot(m15,"M15")
png(paste0(outdir,"/pk/m15.png"), width=800, height=600)
pm15
dev.off()

m16 <- pk(m16,"24")
pm16 <- plot(m16,"M16")
png(paste0(outdir,"/pk/m16.png"), width=800, height=600)
pm16
dev.off()

m17 <- pk(m17,"24")
pm17 <- plot(m17,"M17")
png(paste0(outdir,"/pk/m17.png"), width=800, height=600)
pm17
dev.off()

m18 <- pk(m18,"12")
pm18 <- plot(m18,"M18")
png(paste0(outdir,"/pk/m18.png"), width=800, height=600)
pm18
dev.off()

png(paste0(outdir,"/pk/all.png"), width=1800, height=1200)
plot_grid(pm1,pm2,pm3,pm4,pm10,pm12,pm13,pm14,pm15,pm16,pm17,pm18, ncol=4)
dev.off()




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

# Asegurar formato numérico
sweep.stats$pK <- as.numeric(as.character(sweep.stats$pK))
sweep.stats$BCreal <- as.numeric(as.character(sweep.stats$BCreal))

# Extraer pK óptimo con máximo BCreal
pk.max <- sweep.stats$pK[which.max(sweep.stats$BCreal)]
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
  reuse.pANN = FALSE,
  sct = FALSE
)

  
  return(seurat_obj)
}


smpR065_DsRedN_KO2 <- process_doublets_single_with_rate(smpR065_DsRedN_KO2, nPCs = 14, pN = 0.25, doublet_rate_base = 0.075)
smpR065_DsRedN_KO3 <- process_doublets_single_with_rate(smpR065_DsRedN_KO3, nPCs = 19, pN = 0.25, doublet_rate_base = 0.075)
smpR065_DsRedP_KO2 <- process_doublets_single_with_rate(smpR065_DsRedP_KO2, nPCs = 16, pN = 0.25, doublet_rate_base = 0.075)
smpR065_DsRedP_KO3 <- process_doublets_single_with_rate(smpR065_DsRedP_KO3, nPCs = 12, pN = 0.25, doublet_rate_base = 0.075)

smR039a_WT1_DsRedn <- process_doublets_single_with_rate(smR039a_WT1_DsRedn, nPCs = 15, pN = 0.25, doublet_rate_base = 0.075)
smR039a_WT1_DsRedp <- process_doublets_single_with_rate(smR039a_WT1_DsRedp, nPCs = 17, pN = 0.25, doublet_rate_base = 0.075)

smR039b_WT2_DsRedn <- process_doublets_single_with_rate(smR039b_WT2_DsRedn, nPCs = 14, pN = 0.25, doublet_rate_base = 0.075)
smR039b_WT2_DsRedp <- process_doublets_single_with_rate(smR039b_WT2_DsRedp, nPCs = 11, pN = 0.25, doublet_rate_base = 0.075)




saveRDS(smpR065_DsRedN_KO2,paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedN-KO2.2.rds"))
saveRDS(smpR065_DsRedN_KO3,paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedN-KO3.2.rds"))
saveRDS(smpR065_DsRedP_KO2,paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedP-KO2.2.rds"))
saveRDS(smpR065_DsRedP_KO3,paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedP-KO3.2.rds"))

saveRDS(smR039a_WT1_DsRedn,paste0(datadir,"SoupX/SmR039a.9/Seurat_WT1-DsRedn.2.rds"))
saveRDS(smR039a_WT1_DsRedp,paste0(datadir,"SoupX/SmR039a.9/Seurat_WT1-DsRedp.2.rds"))

saveRDS(smR039b_WT2_DsRedn,paste0(datadir,"SoupX/SmR039b.9/Seurat_WT2-DsRedn.2.rds"))
saveRDS(smR039b_WT2_DsRedp,paste0(datadir,"SoupX/SmR039b.9/Seurat_WT2-DsRedp.2.rds"))





smpR065_DsRedN_KO2 <- readRDS(paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedN-KO2.2.rds"))
smpR065_DsRedN_KO3 <- readRDS(paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedN-KO3.2.rds"))
smpR065_DsRedP_KO2 <- readRDS(paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedP-KO2.2.rds"))
smpR065_DsRedP_KO3 <- readRDS(paste0(datadir,"SoupX/smR065/Seurat_smpR065_DsRedP-KO3.2.rds"))

smR039a_WT1_DsRedn <- readRDS(paste0(datadir,"SoupX/SmR039a.9/Seurat_WT1-DsRedn.2.rds"))
smR039a_WT1_DsRedp <- readRDS(paste0(datadir,"SoupX/SmR039a.9/Seurat_WT1-DsRedp.2.rds"))

smR039b_WT2_DsRedn <- readRDS(paste0(datadir,"SoupX/SmR039b.9/Seurat_WT2-DsRedn.2.rds"))
smR039b_WT2_DsRedp <- readRDS(paste0(datadir,"SoupX/SmR039b.9/Seurat_WT2-DsRedp.2.rds"))

# Plots



psmpR065_DsRedN_KO2 <- plot(smpR065_DsRedN_KO2,"DsRedN_KO2")
png(paste0(outdir,"/SoupX.Doublet/DsRedN_KO2.png"), width=800, height=600)
psmpR065_DsRedN_KO2
dev.off()

psmpR065_DsRedN_KO3 <- plot(smpR065_DsRedN_KO3,"DsRedN_KO3")
png(paste0(outdir,"/SoupX.Doublet/DsRedN_KO3.png"), width=800, height=600)
psmpR065_DsRedN_KO3
dev.off()

psmpR065_DsRedP_KO2 <- plot(smpR065_DsRedP_KO2,"DsRedP_KO2")
png(paste0(outdir,"/SoupX.Doublet/DsRedP_KO2.png"), width=800, height=600)
psmpR065_DsRedP_KO2
dev.off()

psmpR065_DsRedP_KO3 <- plot(smpR065_DsRedP_KO3,"DsRedP_KO3")
png(paste0(outdir,"/SoupX.Doublet/DsRedP_KO3.png"), width=800, height=600)
psmpR065_DsRedP_KO3
dev.off()


psmR039a_WT1_DsRedn <- plot(smR039a_WT1_DsRedn,"WT1_DsRedn")
png(paste0(outdir,"/SoupX.Doublet/WT1_DsRedn.png"), width=800, height=600)
psmR039a_WT1_DsRedn
dev.off()

psmR039a_WT1_DsRedp <- plot(smR039a_WT1_DsRedp,"WT1_DsRedp")
png(paste0(outdir,"/SoupX.Doublet/WT1_DsRedp.png"), width=800, height=600)
psmR039a_WT1_DsRedp
dev.off()



psmR039b_WT2_DsRedn <- plot(smR039b_WT2_DsRedn,"WT2_DsRedn")
png(paste0(outdir,"/SoupX.Doublet/WT2_DsRedn.png"), width=800, height=600)
psmR039b_WT2_DsRedn
dev.off()


psmR039b_WT2_DsRedp <- plot(smR039b_WT2_DsRedp,"WT2_DsRedp")
png(paste0(outdir,"/SoupX.Doublet/WT2_DsRedp.png"), width=800, height=600)
psmR039b_WT2_DsRedp
dev.off()



png(paste0(outdir,"/SoupX.Doublet/doublets.png"), width=2200, height=1200)
plot_grid(psmpR065_DsRedN_KO2,
psmpR065_DsRedN_KO3,
psmpR065_DsRedP_KO2,
psmpR065_DsRedP_KO3,
psmR039a_WT1_DsRedn,
psmR039a_WT1_DsRedp,
psmR039b_WT2_DsRedn,
psmR039b_WT2_DsRedp, ncol=4)
dev.off()