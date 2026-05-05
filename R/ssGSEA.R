library(Seurat)
library(escape)
library(msigdbr)
library(dplyr)
library(GSVA)




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

data <- readRDS(paste0(rsdir,"objects/data.Clusterized.Round3.rds"))


Idents(data) <- "Clustering.Round3"

# Clusters of interest
macro_clusters <- c(
  "Ly6cHi Monocytes", "Ly6cLo Monocytes", "Early IFN|MHCII-TAMs",
    "Trem1|Ptgs2|Plaur|Celc4e Mac",
    "Mrc1|C1qc|Cbr2|Gas6 Mac",
    "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac",
    "Npr2|Actn1 Mac",
    "Mmp9|Ctsk Mac",
    "IFN Mac",
    "Fn1|Vegfa Mac",
    "MHCII|Siglec Mac",
    "MHCII|Ccl12 Mac"
)

data <- subset(data, subset = Clustering.Round3 %in% macro_clusters)

# Hallmark


msig_h <- msigdbr(
  species = "Mus musculus",
  category = "H"
)

genesets_hallmark <- split(
  msig_h$gene_symbol,
  msig_h$gs_name
)


DefaultAssay(data) <- "SCT"
library(BiocParallel)
register(MulticoreParam(workers = 8))


data_hallmark <- runEscape(
  data,
  gene.sets = genesets_hallmark,
  min.size = 10
)




pdf(paste0(outdir,"/ssGSEA/heatmap_KO_WT2.pdf"), width=10, height=18)

p <- heatmapEnrichment(
  data_hallmark, 
  group.by = "Clustering.Round3",
  gene.set.use = "all",
  assay = "escape"
)

# Si heatmapEnrichment devuelve un objeto ggplot, se puede modificar así:
p + theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )

dev.off()



# -------------------------------
# 1️⃣ Subset de KO y WT según metadata
# -------------------------------
meta_subset <- data_hallmark@meta.data

# Filas de las células KO y WT
cells_KO <- rownames(meta_subset)[meta_subset$group == "KO"]
cells_WT <- rownames(meta_subset)[meta_subset$group == "WT"]

# Subset del objeto Seurat
data_KO <- subset(data_hallmark, cells = cells_KO)
data_WT <- subset(data_hallmark, cells = cells_WT)

# -------------------------------
# 2️⃣ Heatmap KO
# -------------------------------
pdf(paste0(outdir,"/ssGSEA/heatmap_KO.pdf"), width=10, height=18)
heatmapEnrichment(
  data_KO,
  group.by = "Clustering.Round3",
  gene.set.use = "all",
  assay = "escape"
)
dev.off()

# -------------------------------
# 3️⃣ Heatmap WT
# -------------------------------
pdf(paste0(outdir,"/ssGSEA/heatmap_WT.pdf"), width=10, height=18)
heatmapEnrichment(
  data_WT,
  group.by = "Clustering.Round3",
  gene.set.use = "all",
  assay = "escape"
)
dev.off()




library(dplyr)
library(reshape2)
library(ggplot2)
library(Matrix)

# -------------------------------
# 1️⃣ Obtener matriz de scores ESCAPE
# -------------------------------
escape_mat <- GetAssayData(data_hallmark, assay = "escape", layer = "data")  # filas = pathways, columnas = células

# -------------------------------
# 2️⃣ Metadata correcta
# -------------------------------
meta <- data_hallmark@meta.data

# Asegurarnos de que coincidan las células
meta_subset <- meta[colnames(escape_mat), , drop = FALSE]

# Crear columna cluster + group si no existe
meta_subset$cluster_group <- paste(meta_subset$Clustering.Round3, meta_subset$group, sep = "_")

# -------------------------------
# 3️⃣ Convertir matriz a formato largo
# -------------------------------
escape_df_t <- as.data.frame(t(as.matrix(escape_mat)))  # filas = células, columnas = pathways
escape_df_t$Cluster <- meta_subset$Clustering.Round3
escape_df_t$Group <- meta_subset$group  # KO o WT

# Calcular promedio por cluster
escape_avg <- escape_df_t %>%
  group_by(Cluster, Group) %>%
  summarise(across(where(is.numeric), mean))

# -------------------------------
# 4️⃣ Preparar matriz para heatmap
# -------------------------------
escape_heatmap <- as.data.frame(t(escape_avg[, -c(1,2)]))  # filas = pathways, columnas = cluster
colnames(escape_heatmap) <- paste(escape_avg$Cluster, escape_avg$Group, sep = "_")

# -------------------------------
# 5️⃣ Convertir a formato largo para ggplot
# -------------------------------
escape_melt <- melt(as.matrix(escape_heatmap))
colnames(escape_melt) <- c("Pathway", "Cluster_Group", "Enrichment")

# Separar Cluster y Group
escape_melt <- escape_melt %>%
  mutate(
    Cluster = sub("_(KO|WT)$", "", Cluster_Group),
    Group = sub("^.*_", "", Cluster_Group)
  )

# Ordenar Pathways por promedio global
pathway_order <- rev(sort(tapply(escape_melt$Enrichment, escape_melt$Pathway, mean)))
escape_melt$Pathway <- factor(escape_melt$Pathway, levels = names(pathway_order))

# -------------------------------
# 6️⃣ Heatmap comparativo KO vs WT
# -------------------------------
pdf(paste0(outdir,"/ssGSEA/heatmap_KO_WT.pdf"), width = 18, height = 10)

ggplot(escape_melt, aes(x = Cluster, y = Pathway, fill = Enrichment)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("purple", "black", "yellow")) +  # misma escala para comparar
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  facet_wrap(~Group, nrow = 1, scales = "free_x")  # KO y WT lado a lado

dev.off()




library(dplyr)
library(reshape2)
library(ggplot2)
library(Matrix)
library(scales)  # para oob = squish

# -------------------------------
# 1️⃣ Obtener matriz ESCAPE
# -------------------------------
escape_mat <- GetAssayData(data_hallmark, assay = "escape", layer = "data")

# Metadata
meta <- data_hallmark@meta.data
meta_subset <- meta[colnames(escape_mat), , drop = FALSE]

# -------------------------------
# 2️⃣ Mediana por cluster y grupo
# -------------------------------
escape_df_t <- as.data.frame(t(as.matrix(escape_mat)))
escape_df_t$Cluster <- meta_subset$Clustering.Round3
escape_df_t$Group <- meta_subset$group

escape_avg <- escape_df_t %>%
  group_by(Cluster, Group) %>%
  summarise(across(where(is.numeric), median))  # ✅ usar mediana

# -------------------------------
# 3️⃣ Matriz combinada para clustering
# -------------------------------
colnames_comb <- paste(escape_avg$Cluster, escape_avg$Group, sep="_")
escape_matrix <- t(escape_avg[, -c(1,2)])
colnames(escape_matrix) <- colnames_comb

# -------------------------------
# 4️⃣ Clustering de pathways (filas)
# -------------------------------
row_clust <- hclust(dist(escape_matrix))
escape_matrix <- escape_matrix[row_clust$order, ]

# -------------------------------
# 5️⃣ Clustering de clusters (columnas) por cluster (ignorar KO/WT)
# -------------------------------
cluster_names <- unique(escape_avg$Cluster)
col_order <- sapply(cluster_names, function(cl) {
  grep(paste0("^", cl, "_"), colnames(escape_matrix))
}) |> unlist()
escape_matrix <- escape_matrix[, col_order]

# -------------------------------
# 6️⃣ Convertir a formato largo
# -------------------------------
escape_melt <- melt(escape_matrix)
colnames(escape_melt) <- c("Pathway", "Cluster_Group", "Enrichment")
escape_melt <- escape_melt %>%
  mutate(
    Cluster = sub("_(KO|WT)$", "", Cluster_Group),
    Group = sub("^.*_", "", Cluster_Group)
  )

escape_melt$Pathway <- factor(escape_melt$Pathway, levels = unique(escape_melt$Pathway))
escape_melt$Cluster <- factor(escape_melt$Cluster, levels = cluster_names)

# -------------------------------
# 7️⃣ Escala de colores global con saturación
# -------------------------------
min_val <- quantile(escape_melt$Enrichment, 0.05, na.rm = TRUE)
max_val <- quantile(escape_melt$Enrichment, 0.95, na.rm = TRUE)

# -------------------------------
# 8️⃣ Heatmap comparativo KO vs WT
# -------------------------------
pdf(paste0(outdir,"/ssGSEA/heatmap_KO_WT_median_green.pdf"), width = 18, height = 10)

ggplot(escape_melt, aes(x = Cluster, y = Pathway, fill = Enrichment)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("white", "lightgreen", "darkgreen"),
    limits = c(min_val, max_val),
    oob = squish  # saturar valores fuera del rango
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y = element_text(size = 8),
    panel.grid = element_blank()
  ) +
  facet_wrap(~Group, nrow = 1, scales = "free_x")  # KO y WT lado a lado

dev.off()

