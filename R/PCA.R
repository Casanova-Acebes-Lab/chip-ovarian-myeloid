

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





## Correlation of Technical and Biological variables with PCs

# Getting a Circadian expression score
circadian_genes <- c("Arntl", "Clock", "Per1", "Per2", "Cry1", "Cry2", "Nr1d1", "Nr1d2", "Dbp") 

# Make sure they are present in your object
circadian_genes <- intersect(circadian_genes, rownames(Circadian))

# Get expression matrix for those genes (RNA assay, normalized data)
expr_circadian <- GetAssayData(Circadian, assay = "RNA", slot = "data")[circadian_genes, ]

# Average circadian expression per cell (can also use module score)
circadian_score <- colMeans(as.matrix(expr_circadian))


Circadian <- AddModuleScore(Circadian, features = list(circadian_genes), name = "Circadian.Score")


# Heatmap of correlations with PCs and technical/biological variables
# Continue variablaes + discrete numerized

covariables <- data.frame(
  nCount_RNA   = Circadian$nCount_decontXcounts,
  nFeature_RNA = Circadian$nFeature_decontXcounts,
  percent.mt   = Circadian$percent.mt,
  circadian_score = circadian_score,
  Phase.ScoreDiff = Circadian$Phase.ScoreDiff,
  ZT = Circadian$ZT.Num  # añadimos ZT directamente
)

# As Numbers
covariables$Time <- as.numeric(factor(Circadian$Time))
covariables$Tag  <- as.numeric(factor(Circadian$Tag))

# Clustering.Simplified as dummy
dummies_clusters <- model.matrix(~ Clustering.Simplified - 1, data = Circadian@meta.data)

# Binding
covariables <- cbind(covariables, dummies_clusters)

# Check
str(covariables)

# PCA Embeddings
pcs <- Embeddings(Circadian, "pca")[, 1:20]

# Spearman correlations
cors_all <- sapply(1:ncol(pcs), function(i) {
  pc <- pcs[, i]
  apply(covariables, 2, function(var) cor(pc, var, method = "spearman"))
})

# PCs as rows
cors_all <- t(cors_all)
rownames(cors_all) <- paste0("PC", 1:ncol(pcs))

colnames(cors_all)

colnames(cors_all) <- c("nCount_RNA","nFeature_RNA","percent.mt","circadian_score",
"Phase.ScoreDiff","ZT","Library_Day/Night","Sample","B Cells","Myeloid",
"Other_cells", "T Cells")

# Heatmap 
png(paste0(outdir,"/Analysis2/QC/heatmap.technical.png"), width=1200, height=800)
pheatmap(
  cors_all,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  main = "Correlation of PCs with technical and biological covariates before correction",
  color = colorRampPalette(c("blue", "white", "red"))(100),
  fontsize_row = 8,
  fontsize_col = 8,
  display_numbers = TRUE,
  number_format = "%.2f",
  number_color = "black",
  angle_col = 45
)
dev.off()



## Regresion model over PC2


metadata <- Circadian@meta.data
metadata$PC1 <- Embeddings(Circadian, "pca")[,1]
metadata$PC2 <- Embeddings(Circadian, "pca")[,2]


lm_PC2 <- lm(PC2 ~ nCount_decontXcounts + nFeature_decontXcounts + percent.mt +
             Clustering.Simplified +
             Phase.ScoreDiff, 
             data = metadata)

summary(lm_PC2)


metadata$PC2_technical <- predict(lm_PC2, newdata = metadata, 
                                  type = "response", 
                                  terms = ~ nCount_RNA_decontXcount + nFeature_decontXcounts + percent.mt)


metadata$PC2_biological <- metadata$PC2 - metadata$PC2_technical


library(ggplot2)
png(paste0(outdir,"/Analysis2/QC/heatmap.technical.vs1.png"), width=1200, height=800)
ggplot(metadata, aes(x = PC1, y = PC2_biological, color = Clustering.Simplified)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(title = "PC1 vs PC2 biológica")
dev.off()

png(paste0(outdir,"/Analysis2/QC/heatmap.technical.vs2.png"), width=1200, height=800)
ggplot(metadata, aes(x = PC1, y = PC2_technical, color = nFeature_RNA)) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  labs(title = "PC1 vs PC2 técnica")
dev.off()


library(car)
vif(lm_PC2) # para ver multicolinealidad
# R² de modelo completo
summary(lm_PC2)$r.squared
# R² parcial: variabilidad explicada solo por técnicos
lm_tech <- lm(PC2 ~ nCount_RNA + nFeature_RNA, data = metadata)
summary(lm_tech)$r.squared
