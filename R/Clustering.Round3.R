
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



# Colors



nora.colors <- c(
  # --- MONOCITOS & TAMs (ROJOS) ---
  "Ly6cHi Monocytes"     = "#FF0000",  # Rojo inflamatorio puro
  "Ly6cLo Monocytes"     = "#FF6A6A",  # Rojo salmón transición
  "Early IFN|MHCII-TAMs"       = "#B22222",  # Rojo vino (pre-TAM IFN)
  "Trem1|Ptgs2|Plaur|Celc4e Mac" = "#EE7942",   # naranja vivo (inflam. activados)
  "Mrc1|C1qc|Cbr2|Gas6 Mac"      = "#FFD92F",   # amarillo brillante (TAM residentes)
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac" = "#4DAF4A", # verde TAM reparadores
  "Npr2|Actn1 Mac"       = "#A6D854",   # verde lima (TAM estructurales)
  "Marco+ Mac"           = "#1E3A8A",   # azul intermedio (scavenger)
  "Mmp9|Ctsk Mac"        = "#00723F",   # verde botella (remodelado matriz)
  "IFN Mac"              = "#C080FF",   # morado claro (TAM interferón maduros)
  "Fn1|Vegfa Mac"        = "#FFA500",   # naranja angiogénico
  "MHCII|Siglec Mac"     = "#1E90FF",   # azul cobalto (antigen presenting)
  "MHCII|Ccl12 Mac"      = "#4682B4",   # azul acero

  # --- CÉLULAS LINFOIDES ---
  "Cd4 Naive"          = "#00BFC4",   # turquesa
  "Cd4 Activated"          = "#FF69B4",   # rosa fuerte
  "Cd8 Effector"         = "#698B69",   # rojo ladrillo
  "Cd4 Th17"                  = "#D9B3FF",   # lila pastel
  "Cd4 Treg"                  = "#36648B",   # lila pastel
  "NK"                   = "#AB82FF",   
  "Mki67+ Tcell"         = "#43CD80",   
  "Activated B cells"    = "#FF1493",   # fucsia intenso
  "B cells"              = "#DC143C",   # rojo intenso

  # --- INNATAS ---
  "Neutrophils"          = "#4876FF",   # azul fuerte
  "DCs"                  = "#87CEEB",   # azul claro
  "Mastocytes"           = "#3CB371"    # verde saturado
)



### Adding Tcells subclustering

t.annotation <- read.csv(paste0(outdir,"/Tcells/subclustering_tcells_annotations.csv"))

# Asegura que el índice son los barcodes
rownames(t.annotation) <- t.annotation$barcode
t.annotation$barcode <- NULL

# Agregar metadata al objeto global
data <- AddMetaData(data, t.annotation)



# Crear la nueva columna como character para que permita nuevos valores
data$Clustering.Round3 <- as.character(data$Clustering.Round2)

# Reasignar según tus subclusters refinados
data$Clustering.Round3[data$Subclustering.T == "Cd4 Activated"] <- "Cd4 Activated"
data$Clustering.Round3[data$Subclustering.T == "Cd4 Naive"] <- "Cd4 Naive"
data$Clustering.Round3[data$Subclustering.T == "Cd4 Th17"] <- "Cd4 Th17"
data$Clustering.Round3[data$Subclustering.T == "Cd4 Treg"] <- "Cd4 Treg"
data$Clustering.Round3[data$Subclustering.T == "Cd8 Effector"] <- "Cd8 Effector"
data$Clustering.Round3[data$Subclustering.T == "Mki67+ Tcell"] <- "Mki67+ Tcell"
data$Clustering.Round3[data$Subclustering.T == "NK"] <- "NK"

# Convertimos de vuelta a factor, ya con todos los niveles presentes
data$Clustering.Round3 <- factor(data$Clustering.Round3)



data$type <- paste(data$group, data$DsRed)



#### Plots 

data <- SetIdent(data, value = "Clustering.Round3") 

pdf(paste0(outdir,"/Clustering.Round3/Umap.pdf"), width=18, height=8)
 DimPlot(data, reduction = "umap", label = FALSE, raster = FALSE, cols = nora.colors)

dev.off()





# Downsampling

set.seed(123)

n_cells <- 8000
cells_to_keep <- c()

for (group in unique(data$tag)) {
  # Subset de células por metadata
  group_cells <- rownames(data@meta.data)[data@meta.data$tag == group]
  
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
table(data_downsampled$tag)



pdf(paste0(outdir,"/Clustering.Round3/Umap.split.downsampled.pdf"), width=22, height=8)
 DimPlot(data, reduction = "umap", split.by="tag", 
 label = FALSE, raster = FALSE, cols = nora.colors, ncol = 4)

dev.off()




pdf(paste0(outdir,"/Clustering.Round3/Umap.split.downsampled2.8kcells.pdf"), width=18, height=10)
 DimPlot(data, reduction = "umap", split.by="type", 
 label = FALSE, raster = FALSE, cols = nora.colors, ncol = 2)

dev.off()




# Quantification


# Proportions por tag
prop_df <- data_downsampled@meta.data %>%
  group_by(type, Clustering.Round3) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(type) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.Round2
prop_df$Clustering.Round3 <- factor(
  prop_df$Clustering.Round3,
  levels = names(nora.colors)
)



pdf(paste0(outdir,"/Clustering.Round3/quantification.type.pdf"), width=24, height=12)
ggplot(prop_df, aes(x = type, y = prop, fill = Clustering.Round3)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Tag",
    y = "Cell proportion (Downsampled to 8k)",
    fill = "Clustering"   # ← aquí cambias la leyenda
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()






prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round3",
	sample_1 = "KO DsRedN", sample_2 = "KO DsRedP",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round3/Proportion.test.KO.pdf"), width=8, height=10)
plot1 <- permutation_plot(prop_test) + ggtitle("Proportion test.DsRedP vs DsRedN. KO samples")
plot1
dev.off()



prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round3",
	sample_1 = "WT DsRedN", sample_2 = "WT DsRedP",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round3/Proportion.test.WT.pdf"), width=8, height=10)
plot2 <- permutation_plot(prop_test) + ggtitle("Proportion test.DsRedP vs DsRedN. WT samples")
plot2
dev.off()




prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round3",
	sample_1 = "WT DsRedP", sample_2 = "KO DsRedP",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round3/Proportion.test.Positive.KOvsWT.pdf"), width=8, height=10)
plot3 <- permutation_plot(prop_test) + ggtitle("Proportion test.KO vs WT. DsRedP samples")
plot3
dev.off()



prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round3",
	sample_1 = "WT DsRedN", sample_2 = "KO DsRedN",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round3/Proportion.test.Negative.KOvsWT.pdf"), width=8, height=10)
plot4 <- permutation_plot(prop_test) + ggtitle("Proportion test.KO vs WT. DsRedN samples")
plot4
dev.off()



pdf(paste0(outdir,"/Clustering.Round3/Proportion.tests.all.clustR3.pdf"), width=30, height=8)
plot_grid(plot1, plot2, plot3, plot4, ncol=4)
dev.off() 

