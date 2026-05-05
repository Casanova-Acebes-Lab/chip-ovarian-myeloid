
library(Seurat)
library(dplyr)
library(scCustomize)
library(ggplot2)
library(RPresto)
library(ggrepel)
library(stringr)
library(DESeq2)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(clusterProfiler)
library(enrichplot)
library(SingleR)
library(cowplot)





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

data <- readRDS(paste0(rsdir,"objects/data.Clustering.Round1.rds"))




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
  "Stab1|Axl Mac",
  "Mki67 IFN Mac",
  "Mki67|Cstk|Mmp9|S100a4 Mac",
  "Nlrp3|Vegfa Mac",

  # LINFOIDES
  "Cd8 Effector",
  "Cd8 Cytotoxic",
  "Cd4 Naive",
  "Treg",
  "Activated B cells",
  "B cells",
  
  # INNATAS
  "Neutrophils",
  "DCs",
  "PDcs",
  "NK",
  "Mastocytes"
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
  "Stab1|Axl Mac"                  = "#8A2BE2",
  "Mki67 IFN Mac"                  = "#504369",
  "Mki67|Cstk|Mmp9|S100a4 Mac"    = "#db44a9",
  "Nlrp3|Vegfa Mac"                = "#5f3121",

  # T CELLS - NUEVOS COLORES
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

  # B CELLS
  "Activated B cells"               = "#FF1493",
  "B cells"                         = "#DC143C",

  # INNATAS
  "Neutrophils"                     = "#4876FF",
  "DCs"                             = "#87CEEB",
  "PDcs"                            = "#7FFFD4",
  "NK"                               = "#AB82FF",
  "ILC"                              = "#FFA500",
  "Mastocytes"                      = "#3CB371"
)

## Extracting T cells

data <- subset(
  data,
  subset = Clustering.Round2 %in% c(
    "Cd8 Effector",
    "Cd8 Cytotoxic",
    "Cd4 Naive",
    "Treg",
    "NK"
  )
)




# Subslustering


# removing contaminting myeloid cells
contam_cells <- WhichCells(data, expression = Msr1 > 0.5 | Fcgr4 > 0.5 | Clec4n > 0.5 | Hck > 0.5)
data <- subset(data, cells = setdiff(Cells(data), contam_cells))


# Identificación de genes variables ya está hecha por SCTransform
# PCA
data <- RunPCA(data, features = VariableFeatures(data))

# Vecindad y clustering
data <- FindNeighbors(data, dims = 1:30) %>%
         FindClusters(resolution = 0.6)

# UMAP
data <- RunUMAP(data, dims = 1:30)

png(paste0(outdir,"/Tcells/Umap.subclustering.png"), width=1200, height=800)
DimPlot(data, reduction = "umap", label = TRUE)
dev.off()



# Establecer identidades con Seurat >=4
Idents(data) <- "seurat_clusters"

# Encontrar marcadores
markers <- FindAllMarkers(
  object = data,
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
FeaturePlot_scCustom(data, features = c("Cd8a", "Cd4" , "Pdcd1", "Oas2", "Ifit2", "Ifit3"), split.by = "group")
dev.off()



png(paste0(outdir,"/Tcells/Umap.cd.group.png"), width=1600, height=3200)
FeaturePlot_scCustom(data, features = c("Cd8a", "Cd4" , "Fas", "Fasl", "Il17a", "Rorc"), split.by = "group")
dev.off()


png(paste0(outdir,"/Tcells/Umap.violin.markers.png"), width=1600, height=1600)
VlnPlot(
  object = data,
  features = c("Cd3e", "Cd4", "Cd8a", 
               "Msr1", "Fcgr4", "C3", "Hck", "Clec4n", "Csf1r"),
  group.by = "seurat_clusters",
  pt.size = 0
)
dev.off()



data@meta.data$Subcluster <- NA
data@meta.data$Subcluster[data$seurat_clusters == 0] <- "Cd8 Exhausted"
data@meta.data$Subcluster[data$seurat_clusters == 1] <- "Treg"
data@meta.data$Subcluster[data$seurat_clusters == 2] <- "Cd8 memory-like"
data@meta.data$Subcluster[data$seurat_clusters == 3] <- "Cd8 memory-like"
data@meta.data$Subcluster[data$seurat_clusters == 4] <- "Cd4 Effector-Memory"
data@meta.data$Subcluster[data$seurat_clusters == 5] <- "Cd4 Activated"
data@meta.data$Subcluster[data$seurat_clusters == 6] <- "Cd4 Activated"
data@meta.data$Subcluster[data$seurat_clusters == 7] <- "Cd8 Cytotoxic"
data@meta.data$Subcluster[data$seurat_clusters == 8] <- "NK"
data@meta.data$Subcluster[data$seurat_clusters == 9] <- "Cd8 Exhausted"
data@meta.data$Subcluster[data$seurat_clusters == 10] <- "NK"
data@meta.data$Subcluster[data$seurat_clusters == 11] <- "Cd8 Effector"
data@meta.data$Subcluster[data$seurat_clusters == 12] <- "Cd4 Effector-Memory"
data@meta.data$Subcluster[data$seurat_clusters == 13] <- "Proliferating Tcells"
data@meta.data$Subcluster[data$seurat_clusters == 14] <- "Tregs activated"
data@meta.data$Subcluster[data$seurat_clusters == 15] <- "Th17"
data@meta.data$Subcluster[data$seurat_clusters == 16] <- "Cd4 Effector-Memory"
data@meta.data$Subcluster[data$seurat_clusters == 17] <- "Tgd"
data@meta.data$Subcluster[data$seurat_clusters == 18] <- "ILC"





data <- SetIdent(data, value = "Subcluster")
pdf(paste0(outdir,"/Tcells/Umap.Clusterized.pdf"), width=12, height=8)
 DimPlot(data, reduction = "umap", label = TRUE)

dev.off()


data <- SetIdent(data, value = "Subcluster")
pdf(paste0(outdir,"/Tcells/Umap.Clusterized.group.pdf"), width=16, height=8)
 DimPlot(data, reduction = "umap", split.by = "group", label = FALSE,
 cols=nora.colors)

dev.off()


data <- SetIdent(data, value = "Subcluster")
pdf(paste0(outdir,"/Tcells/Umap.Clusterized.DsRed.pdf"), width=16, height=8)
 DimPlot(data, reduction = "umap", split.by = "DsRed", label = FALSE,cols=nora.colors)

dev.off()



data <- SetIdent(data, value = "Subcluster")
pdf(paste0(outdir,"/Tcells/Umap.Clusterized.tag.pdf"), width=24, height=16)
 DimPlot(data, reduction = "umap", split.by = "tag", label = TRUE, ncol=4)

dev.off()


## Saving annotation

t.annotation <- data.frame(
  barcode = colnames(data),
  Subclustering.T = data$Subcluster
)

write.csv(
  t.annotation,
  file = paste0(outdir,"/Tcells/subclustering_tcells_annotations.csv"),
  row.names = FALSE
)


saveRDS(data,paste0(rsdir,"tcells.clusterized.rds"))


data <- readRDS(paste0(rsdir,"tcells.clusterized.rds"))

## Quantification


data$type <- paste0(data$DsRed, data$group)

# Proportions por tag
prop_df <- data@meta.data %>%
  group_by(tag, Subcluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tag) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.Round2
prop_df$Subcluster <- factor(
  prop_df$Subcluster,
  levels = names(nora.colors)
)

# Stacked barplot
pdf(paste0(outdir,"/Tcells/quantification.tags.pdf"), width=24, height=12)
ggplot(prop_df, aes(x = tag, y = prop, fill = Subcluster)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(x = "Tag", y = "Cell proportion (Downsampled to 8k)", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




# Proportions por type
prop_df <- data@meta.data %>%
  group_by(type, Subcluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(type) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.Round2
prop_df$Subcluster <- factor(
  prop_df$Subcluster,
  levels = names(nora.colors)
)

prop_df$type <- factor(
  prop_df$type,
  levels = c("DsRedNKO", "DsRedPKO", "DsRedNWT", "DsRedPWT")
)

# Stacked barplot
pdf(paste0(outdir,"/Tcells/quantification.type.pdf"), width=12, height=12)
ggplot(prop_df, aes(x = type, y = prop, fill = Subcluster)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(x = "Tag", y = "Cell proportion (Downsampled to 8k)", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()




### DEG 



# Lista para guardar resultados por subcluster
results_list <- list()

subclusters <- unique(data@meta.data$Subcluster)

for (sc in subclusters) {
  
  # Filtrar solo las células de este subcluster
  cells_sc <- WhichCells(data, expression = Subcluster == sc)
  data_sc <- subset(data, cells = cells_sc)
  
  # Verificar que hay células de ambos grupos
  if(length(unique(data_sc$group)) < 2) next
  
  # MAST differential expression
  sc_markers <- FindMarkers(
    data_sc,
    ident.1 = "KO",
    ident.2 = "WT",
    group.by = "group",
    test.use = "MAST",
    slot = "counts",
    latent.vars = c("nCount_RNA", "percent.mt"),
    only.pos=TRUE,
    logfc.threshold = 0.3,
    min.pct = 0.1  # usando counts crudos
  )
  
  # Agregar columna del subcluster
  sc_markers$Subcluster <- sc
  
  # Guardar
  results_list[[sc]] <- sc_markers
}

# Unir todos los resultados en un data.frame
differential_results <- bind_rows(results_list, .id = "Cluster")


fas_genes <- differential_results[grep("^Fas$|^Fasl$", rownames(differential_results)), ]
fas_genes


pdf(paste0(outdir,"/Tcells/fas.pdf"), width=24, height=12)

VlnPlot(
  data,
  features = c("Fas", "Fasl"),
  group.by = "Subcluster",
  pt.size = 0,
  split.by= "group"
)
dev.off()