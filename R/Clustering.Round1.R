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



# Read clusterized 1

data <- readRDS(paste0(rsdir,"/objects/seurat_QC2_SCT_Harmony_res08_validated.SingleR.clusterized.rds"))


# Removing cluster 13. Fibroblasts and low quality 

data <- subset(data, subset = SCT_snn_res.0.8 != "13")

# Adding Monos subclustering to the main data

Monos <- read.csv(paste0(outdir,"/Subclustering.Monocytes/subclustering_monocytes_annotations.csv"))


# Crear columna Macros en data@meta.data
data@meta.data$Monos <- NA

# Asignar los clusters de tam a las células correspondientes en data
data@meta.data$Monos[match(Monos$barcode, rownames(data@meta.data))] <- Monos$Subclustering.Mon

data@meta.data$Clustering.Round1 <- NA

data@meta.data$Clustering.Round1 <- ifelse(
  !is.na(data@meta.data$Monos),
  as.character(data@meta.data$Monos),
  as.character(data@meta.data$Cluster)
)


sum(Monos$barcode %in% rownames(data@meta.data))

head(data)


data@meta.data$Clustering.Round1[data@meta.data$Clustering.Round1 == "Early TAMs"] <- "Early IFN TAMs"
data@meta.data$Clustering.Round1[data@meta.data$Clustering.Round1 == "Ly6cHi Monocytes"] <- "Ly6c2Hi Monocytes"
data@meta.data$Clustering.Round1[data@meta.data$Clustering.Round1 == "Ly6cLo Monocytes"] <- "Ly6c2Lo Monocytes"
data@meta.data$Clustering.Round1[data@meta.data$Clustering.Round1 == "Monocytes Ly6c2Hi"] <- "Ly6c2Hi Monocytes"
data@meta.data$Clustering.Round1[data@meta.data$Clustering.Round1 == "Monocytes Ly6c2Lo"] <- "Ly6c2Lo Monocytes"





data <- SetIdent(data, value = "Clustering.Round1")
png(paste0(outdir,"/Clustering.Round1/Umap2.png"), width=1200, height=1200)
DimPlot(data, reduction = "umap", raster = FALSE, label=T, repel=T, label.size=4) + NoLegend()
dev.off()


png(paste0(outdir,"/Clustering.Round1/Umap.split2.png"), width=3800, height=2200)
DimPlot(data, reduction = "umap", raster = FALSE, label=T, repel=T, label.size=4, split.by = "tag",
ncol = 4, pt.size = 0.7) + NoLegend()
dev.off()




# Adding Mki67 subclustering to the main data

Mki67 <- read.csv(paste0(outdir,"/Subclustering.Mi67Mac/subclustering_Mki67.csv"))


# Crear columna Macros en data@meta.data
data@meta.data$Mki67 <- NA

# Asignar los clusters de tam a las células correspondientes en data
data@meta.data$Mki67[match(Mki67$barcode, rownames(data@meta.data))] <- Mki67$Subclustering.Mki67

data@meta.data$Clustering.Round2 <- NA

data@meta.data$Clustering.Round2 <- ifelse(
  !is.na(data@meta.data$Mki67),
  as.character(data@meta.data$Mki67),
  as.character(data@meta.data$Clustering.Round1)
)


sum(Mki67$barcode %in% rownames(data@meta.data))

head(data)


# Adding Tcells subclustering to the main data

Tcells <- read.csv(paste0(outdir,"/Tcells/subclustering_tcells_annotations.csv"))


# Crear columna Macros en data@meta.data
data@meta.data$Tcells <- NA

# Asignar los clusters de tam a las células correspondientes en data
data@meta.data$Tcells[match(Tcells$barcode, rownames(data@meta.data))] <- Tcells$Subclustering.T

data@meta.data$Clustering.Round3 <- NA

data@meta.data$Clustering.Round3 <- ifelse(
  !is.na(data@meta.data$Tcells),
  as.character(data@meta.data$Tcells),
  as.character(data@meta.data$Clustering.Round2)
)


sum(Tcells$barcode %in% rownames(data@meta.data))

head(data)




table(data$Clustering.Round3)




### Generic Clustering 

data@meta.data$Clustering.wide <- NA



data@meta.data$Clustering.wide <- NA
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Ly6c2Hi Monocytes"] <- "Monocytes"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Ly6c2Lo Monocytes"] <- "Monocytes"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "IFN Mac"] <- "IFN Mac"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Early IFN TAMs"] <- "Early IFN TAMs"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Trem1|Ptgs2|Plaur|F10 Mac"] <- "Macrophages"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "MHCII|Mgl2 Mac"] <- "Macrophages"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "MHCII|Siglec Mac"] <-  "Macrophages"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Nrp2|Emp1 Mac"] <- "Macrophages"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Mmp9|Ctsk Mac"] <- "Macrophages"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Gas6|Folr2 Mac"] <- "Macrophages"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Saa3 Mac"] <- "Macrophages"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Arg1|Spp1|Mmp12|Il1a Mac"] <- "Macrophages"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Mki67 IFN Mac"] <- "IFN Mac"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Nlrp3|Vegfa Mac"] <- "Macrophages"

data@meta.data$Clustering.wide[data$Clustering.Round3 =="Mki67|Cstk|Mmp9|S100a4 Mac"] <- "Macrophages"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Cd8 Exhausted"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Cd8 memory-like"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Cd8 Cytotoxic"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Cd8 Effector"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Cd4 Effector-Memory"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Cd4 Activated"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Treg"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Tregs activated"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Th17"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Tgd"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Proliferating Tcells"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Activated B cells"] <- "B Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "B cells"] <- "B Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Cd4 Naive"] <- "T Cells"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Neutrophils"] <- "Neutrophils"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "DCs"] <- "DCs"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "PDcs"] <- "DCs"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "NK"] <- "NK"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "Mastocytes"] <- "Mastocytes"
data@meta.data$Clustering.wide[data$Clustering.Round3 == "ILC"] <- "ILC"


  



data <- SetIdent(data, value = "Clustering.Round3")
png(paste0(outdir,"/Clustering.Round1/Umap3.png"), width=1200, height=1200)
DimPlot(data, reduction = "umap", raster = FALSE, label=T, repel=T, label.size=4) + NoLegend()
dev.off()


png(paste0(outdir,"/Clustering.Round1/Umap.split3.png"), width=3200, height=2200)
DimPlot(data, reduction = "umap", raster = FALSE, label=T, repel=T, label.size=4, split.by = "tag",
ncol = 4) + NoLegend()
dev.off()



data <- SetIdent(data, value = "Clustering.wide")
png(paste0(outdir,"/Clustering.Round1/Umap.wide.png"), width=1200, height=1200)
DimPlot(data, reduction = "umap", raster = FALSE, label=T, repel=T, label.size=4) + NoLegend()
dev.off()


png(paste0(outdir,"/Clustering.Round1/Umap.split.wide.png"), width=3200, height=2200)
DimPlot(data, reduction = "umap", raster = FALSE, label=T, repel=T, label.size=4, split.by = "tag",
ncol = 4) + NoLegend()
dev.off()





# Establecer identidades con Seurat >=4
Idents(data) <- "Clustering.Round3"

# Encontrar marcadores
markers <- FindAllMarkers(
  object =data,
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



saveRDS(data, paste0(rsdir,"/objects/data.Clustering.Round1.rds"))

data <- readRDS(paste0(rsdir,"/objects/data.Clustering.Round1.rds"))



####### Figures


cluster_order <- c(
  # MONOCITOS & TAMs
  "Ly6c2Hi Monocytes",
  "Ly6c2Lo Monocytes",
  "IFN Mac",
  "Early IFN TAMs",
  "Trem1|Ptgs2|Plaur|F10 Mac",
  "MHCII|Mgl2 Mac",
  "MHCII|Siglec Mac",
  "Nrp2|Emp1 Mac",
  "Mmp9|Ctsk Mac",
  "Gas6|Folr2 Mac",
  "Saa3 Mac",
  "Arg1|Spp1|Mmp12|Il1a Mac",
  "Mki67 IFN Mac",
  "Mki67|Cstk|Mmp9|S100a4 Mac",
  "Nlrp3|Vegfa Mac",

  # LINFOIDES
"Cd8 Exhausted",
  "Cd8 memory-like",
  "Cd8 Cytotoxic",
  "Cd8 Effector",

  "Cd4 Effector-Memory",
  "Cd4 Activated",
  "Treg",
  "Tregs activated",
  "Th17",
  "Tgd",
  "Proliferating Tcells",
  "Activated B cells",
  "B cells",
  "Cd4 Naive",
  
  # INNATAS
  "Neutrophils",
  "DCs",
  "PDcs",
  "NK",
  "Mastocytes",
  "ILC"
)

nora.colors <- c(

  # MONOCITOS & TAMs
  "Ly6c2Hi Monocytes"             = "#FF0000",
  "Ly6c2Lo Monocytes"             = "#FF6347",
  "IFN Mac"                        = "#C080FF",
  "Early IFN TAMs"                = "#FFA07A",
  "Trem1|Ptgs2|Plaur|F10 Mac"     = "#EE7942",
  "MHCII|Mgl2 Mac"                 = "#4682B4",
  "MHCII|Siglec Mac"               = "#1E90FF",
  "Nrp2|Emp1 Mac"                  = "#A6D854",
  "Mmp9|Ctsk Mac"                  = "#00723F",
  "Gas6|Folr2 Mac"                 = "#FFD700",
  "Saa3 Mac"                        = "#FFA500",
  "Arg1|Spp1|Mmp12|Il1a Mac"      = "#4DAF4A",
  "Mki67 IFN Mac"                 = "#504369",
  "Mki67|Cstk|Mmp9|S100a4 Mac"    = "#db44a9",
  "Nlrp3|Vegfa Mac"               = "#5f3121",

  # LINFOIDES
"Cd8 Exhausted"                  = "#1F77B4",  # azul fuerte
  "Cd8 memory-like"                = "#17BECF",  # celeste
  "Cd8 Cytotoxic"                  = "#2CA02C",  # verde
  "Cd8 Effector"                   = "#D62728",  # rojo intenso

  "Cd4 Effector-Memory"            = "#9467BD",  # morado
  "Cd4 Activated"                  = "#8C564B",  # marrón
  "Treg"                            = "#FF7F0E",  # naranja
  "Tregs activated"                 = "#BCBD22",  # verde oliva
  "Th17"                            = "#7F7F7F",  # gris
  "Tgd"                             = "#AEC7E8",  # celeste claro
  "Proliferating Tcells"            = "#F7B6D2",  # rosa pálido
  "Activated B cells"             = "#FF1493",
  "B cells"                       = "#DC143C",
  "Cd4 Naive"                       = "#FF69B4",

  # INNATAS
  "Neutrophils"                   = "#4876FF",
  "DCs"                           = "#87CEEB",
  "PDcs"                          = "#7FFFD4",
  "NK"                            = "#AB82FF",
  "Mastocytes"                    = "#3CB371",
  "ILC"                              = "#FFA500"
)




pdf(paste0(outdir,"/Clustering.Round1/Umap_Markers.pdf"), width=24, height=20)
FeaturePlot_scCustom(data, features= c("Ifit3","Ifit1", 
"Ifit2", "Ly6c2", "Gas6", "Ly6c2", "Ciita", "Siglece", "Cp", "Il1b", "Trem1", "Arg1", "Nlrp3"), reduction = "umap")
dev.off()


umap <- Embeddings(data, "umap")
cells_remove <- rownames(data@meta.data)[
  data$Clustering.Round2 == "Nlrp3|Vegfa Mac" &
  umap[,1] >= 0
]

data<- subset(data, cells = setdiff(colnames(data), cells_remove))



data <- SetIdent(data, value = "Clustering.Round3")
pdf(paste0(outdir,"/Clustering.Round1/Umap.clusterized.Round4.pdf"), width=12, height=12)
DimPlot(data, reduction = "umap", raster = FALSE, label=T, repel=T, label.size=5,
cols = nora.colors, pt.size = 0.6) + NoLegend()
dev.off()


data <- SetIdent(data, value = "Clustering.Round3")
pdf(paste0(outdir,"/Clustering.Round1/Umap.clusterized.Round4.legend.pdf"), width=18, height=12)
DimPlot(data, reduction = "umap", raster = FALSE, 
cols = nora.colors, pt.size = 0.6)
dev.off()



# Downsampling

set.seed(123)

n_cells <- 8000
cells_to_keep <- c()

for (group in unique(data$sample_id)) {
  # Subset de células por metadata
  group_cells <- rownames(data@meta.data)[data@meta.data$sample_id == group]
  
  if (length(group_cells) <= n_cells) {
    sampled_cells <- group_cells
  } else {
    sampled_cells <- sample(group_cells, n_cells)
  }
  
  cells_to_keep <- c(cells_to_keep, sampled_cells)
}

# Ahora sí crear el subset
data_downsampled <- subset(data, cells = cells_to_keep)

# Revisar número de células por grupo
table(data_downsampled$sample_id)


data_downsampled <- SetIdent(data_downsampled, value = "Clustering.Round3")
pdf(paste0(outdir,"/Clustering.Round1/Umap.split.Round4.Downsampled.8k.pdf"), width=32, height=14)
 DimPlot(data_downsampled, reduction = "umap", split.by="sample_id", 
 label = FALSE, raster = FALSE, cols = nora.colors, ncol = 4, pt.size = 0.5)

dev.off()



pdf(paste0(outdir,"/Clustering.Round1/Umap.split.Round4.Downsampled.16k.chimera.pdf"), width=25, height=18)
 DimPlot(data_downsampled, reduction = "umap", split.by="chimera", 
 label = FALSE, raster = FALSE, cols = nora.colors, ncol = 2, pt.size = 0.5)

dev.off()


data_downsampled <-SetIdent(data_downsampled, value = "group")
pdf(paste0(outdir,"/Clustering.Round1/Umap.split.Round4.Downsampled.64k.group.pdf"), width=12, height=12)
 DimPlot(data_downsampled, reduction = "umap", 
 label = FALSE, raster = FALSE, ncol = 2, pt.size = 0.5)

dev.off()


data_downsampled$type <- paste(data_downsampled$group, data_downsampled$DsRed)


data_downsampled <-SetIdent(data_downsampled, value = "Clustering.Round2")
pdf(paste0(outdir,"/Clustering.Round1/Umap.split.Round4.Downsampled.16k.type.pdf"), width=18, height=12)
 DimPlot(data_downsampled, reduction = "umap", split.by="type",
 label = FALSE, raster = FALSE, ncol = 2, pt.size = 0.5,
 cols= nora.colors)

dev.off()


### Quantification



# Proportions por tag
prop_df <- data_downsampled@meta.data %>%
  group_by(sample_id, Clustering.Round2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.Round2
prop_df$Clustering.Round2 <- factor(
  prop_df$Clustering.Round2,
  levels = names(nora.colors)
)



pdf(paste0(outdir,"/Clustering.Round1/quantification.sample_id.down.8k.pdf"), width=16, height=12)
ggplot(prop_df, aes(x = sample_id, y = prop, fill = Clustering.Round2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Sample ID",
    y = "Cell proportion (Downsampled to 8k cells)",
    fill = "Clustering"   # ← aquí cambias la leyenda
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# Proportions por tag
prop_df <- data_downsampled@meta.data %>%
  group_by(chimera, Clustering.Round2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(chimera) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.Round2
prop_df$Clustering.Round2 <- factor(
  prop_df$Clustering.Round2,
  levels = names(nora.colors)
)



pdf(paste0(outdir,"/Clustering.Round1/quantification.chimera.down.16k.pdf"), width=16, height=12)
ggplot(prop_df, aes(x = chimera, y = prop, fill = Clustering.Round2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Tag",
    y = "Cell proportion (Downsampled to 16k cells)",
    fill = "Clustering"   # ← aquí cambias la leyenda
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()






# Proportions por type
prop_df <- data_downsampled@meta.data %>%
  group_by(type, Clustering.Round2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(type) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.Round2
prop_df$Clustering.Round2 <- factor(
  prop_df$Clustering.Round2,
  levels = names(nora.colors)
)



pdf(paste0(outdir,"/Clustering.Round1/quantification.type.down.16k.pdf"), width=16, height=12)
ggplot(prop_df, aes(x = type, y = prop, fill = Clustering.Round2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Tag",
    y = "Cell proportion (Downsampled to 16k cells)",
    fill = "Clustering"   # ← aquí cambias la leyenda
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()


data_downsampled$type <- paste(data_downsampled$group, data_downsampled$DsRed)


library("scProportionTest")

prop_test <- sc_utils(data_downsampled)


prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round2",
	sample_1 = "WT DsRedP", sample_2 = "KO DsRedP",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round1/Proportion.test.KO-WT.DsredP.pdf"), width=10, height=8)
plot1 <- permutation_plot(prop_test) + ggtitle("Proportion test KO vs WT DsRedP samples")
plot1
dev.off()




prop_test <- sc_utils(data_downsampled)


prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round2",
	sample_1 = "WT DsRedN", sample_2 = "KO DsRedN",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round1/Proportion.test.KO-WT.DsredN.pdf"), width=10, height=8)
plot2 <- permutation_plot(prop_test) + ggtitle("Proportion test KO vs WT DsRedN samples")
plot2
dev.off()




png(paste0(outdir,"/Clustering.Round1/Umap_Markers.png"), width=1400, height=3000)
FeaturePlot_scCustom(data, features= c("Ifit3","Ifit1", 
"Ifit2", "Il1b", "Vegfa", "Oas2", "Stat3"), reduction = "umap",
split.by = "group", ncol=4)
dev.off()



pdf(paste0(outdir,"/Clustering.Round1/Umap_Markers.pdf"), width=14, height=40)
FeaturePlot_scCustom(data, features= c("Ifit3","Ifit1", 
"Ifit2", "Il1b", "Vegfa", "Oas2", "Arg1", "Trem1"), reduction = "umap",
split.by = "group")
dev.off()





pdf(paste0(outdir,"/Clustering.Round1/Umap_Markers2.pdf"), width=14, height=30)
FeaturePlot_scCustom(data, features= c("Ifnar1","Ifnar2", 
"Stat1", "Stat2", "Stat3", "Irf7", "Irf9"), reduction = "umap",
split.by = "group")
dev.off()




### Wide clustering 




nora.colors.wide <- c(

  # MONOCITOS & TAMs
  "Monocytes"             = "#FF0000",

  "IFN Mac"                        = "#C080FF",
  "Early IFN TAMs"                = "#FFA07A",


  "B Cells"                       = "#CD5C5C",
  "DCs"                          = "#98cee4",
  "Macrophages"                   = "#CD950C",
"T Cells"                       = "#6692b1",
  # INNATAS
  "Neutrophils"                   = "#4876FF",
  "NK"                            = "#8968CD",
  "Mastocytes"                    = "#3CB371",
  "ILC"                              = "#a57317"
)






data_downsampled <- SetIdent(data_downsampled, value = "Clustering.wide")
pdf(paste0(outdir,"/Clustering.Round1/Umap.split.Round4.Downsampled.8k.wide.pdf"), width=32, height=14)
 DimPlot(data_downsampled, reduction = "umap", split.by="sample_id", 
 label = FALSE, raster = FALSE, cols = nora.colors.wide, ncol = 4, pt.size = 0.5)

dev.off()



pdf(paste0(outdir,"/Clustering.Round1/Umap.split.Round4.Downsampled.16k.chimera.wide.pdf"), width=25, height=18)
 DimPlot(data_downsampled, reduction = "umap", split.by="chimera", 
 label = FALSE, raster = FALSE, cols = nora.colors.wide, ncol = 2, pt.size = 0.5)

dev.off()


### Quantification



# Proportions por tag
prop_df <- data_downsampled@meta.data %>%
  group_by(sample_id, Clustering.wide) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(sample_id) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.wide
prop_df$Clustering.wide <- factor(
  prop_df$Clustering.wide,
  levels = names(nora.colors.wide)
)



pdf(paste0(outdir,"/Clustering.Round1/quantification.sample_id.down.8k.wide.pdf"), width=16, height=12)
ggplot(prop_df, aes(x = sample_id, y = prop, fill = Clustering.wide)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors.wide) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Sample ID",
    y = "Cell proportion (Downsampled to 8k cells)",
    fill = "Clustering"   # ← aquí cambias la leyenda
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()



# Proportions por tag
prop_df <- data_downsampled@meta.data %>%
  group_by(chimera, Clustering.wide) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(chimera) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.wide
prop_df$Clustering.wide <- factor(
  prop_df$Clustering.wide,
  levels = names(nora.colors.wide)
)



pdf(paste0(outdir,"/Clustering.Round1/quantification.chimera.down.16k.wide.pdf"), width=16, height=12)
ggplot(prop_df, aes(x = chimera, y = prop, fill = Clustering.wide)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors.wide) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Tag",
    y = "Cell proportion (Downsampled to 16k cells)",
    fill = "Clustering"   # ← aquí cambias la leyenda
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()






library("scProportionTest")



data_downsampled@meta.data$type <- paste(data_downsampled$group, data_downsampled$DsRed, sep = "_")

data_downsampled$type <- factor(
 data_downsampled$type,
  levels = c(
    "WT_DsRedN",
    "WT_DsRedP",
    "KO_DsRedN",
    "KO_DsRedP"
  )
)

table(data_downsampled$type)


prop_test <- sc_utils(data_downsampled)


prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.wide",
	sample_1 = "WT", sample_2 = "KO",
	sample_identity = "group"
)


pdf(paste0(outdir,"/Clustering.Round1/Proportion.test.KO-WT.wide.pdf"), width=10, height=8)
plot1 <- permutation_plot(prop_test) + ggtitle("Proportion test KO vs WT Chimeras")
plot1
dev.off()




prop_test <- sc_utils(data_downsampled)


prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round2",
	sample_1 = "WT DsRedN", sample_2 = "KO DsRedN",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round1/Proportion.test.KO-WT.DsredN.pdf"), width=10, height=8)
plot2 <- permutation_plot(prop_test) + ggtitle("Proportion test KO vs WT DsRedN samples")
plot2
dev.off()


