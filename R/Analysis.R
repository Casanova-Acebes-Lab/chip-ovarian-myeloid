
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




#######

library(ggplot2)
library(dplyr)
library(tidyr)

# 1️⃣ Tabla de conteo de células por cluster x muestra x genotipo
composition <- data@meta.data %>%
  group_by(Clustering.Round2, tag, group) %>%
  summarise(n_cells = n(), .groups = "drop")

head(composition)

# 2️⃣ Pivot para ver proporciones
composition_prop <- composition %>%
  group_by(Clustering.Round2) %>%
  mutate(prop = n_cells / sum(n_cells))

# 3️⃣ Plot de stacked bar por cluster
pdf(paste0(outdir,"/Analysis/cellbycluster.pdf"), width=24, height=12)
ggplot(composition_prop, aes(x=Clustering.Round2, y=prop, fill=tag)) +
  geom_bar(stat="identity") +
  facet_wrap(~group) +
  ylab("Proporción de células") +
  xlab("Cluster") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1)) +
  ggtitle("Composición de células por cluster y muestra")
dev.off()

####



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




data_downsampled <- SetIdent(data_downsampled, value = "Clustering.Round2")
pdf(paste0(outdir,"/Analysis/Umap.tag.Downsampled.8k.pdf"), width=24, height=12)
DimPlot(data_downsampled, reduction = "umap", split.by = "tag", ncol=4, raster = FALSE,
pt.size = 0.6) +  scale_color_manual(values = nora.colors)     
 
dev.off()




pdf(paste0(outdir,"/Analysis/Umap_expansion.downsampled.8k.pdf"), width=14, height=14)
DimPlot_scCustom(data_downsampled, reduction = "umap", group.by="group", 
ggplot_default_colors=T, raster = FALSE, pt.size=0.6) 
dev.off() 



pdf(paste0(outdir,"/Analysis/Umap_expansion2.downsampled8k.pdf"), width=24, height=14)
DimPlot_scCustom(data_downsampled, reduction = "umap", group.by="group", 
ggplot_default_colors=T, raster = FALSE, pt.size=0.3, split.by="DsRed") 
dev.off() 


# Reordenar niveles del factor para que WT quede primero
data_downsampled$group <- factor(data_downsampled$group, levels = c("WT", "KO"))

pdf(paste0(outdir,"/Analysis/Umap_expansion3.downsampled8k.pdf"), width=24, height=14)

DimPlot(
  data_downsampled,
  reduction = "umap",
  group.by = "DsRed",
  cols = c("DsRedN" = "grey", "DsRedP" = "red"),
  pt.size = 0.6,
  split.by = "group",
  raster = FALSE
)

dev.off()



interferon_genes_20 <- c(
  # Señalización y receptores
  "Ifna2", "Oas1a", "Oas2", "Oas3",
  
  # Factores reguladores de IFN
  "Irf7", "Irf9",
  
  # ISGs potentes
  "Isg15", "Mx1", "Ifit1", "Ifit3",
  "Rsad2", "Oas1a", "Oas2", "Oas3",
  
  # Antiviral / inmunoproteasoma
  "Bst2", "Zbp1", "Cxcl10",
  
  # Sensores antivirales
  "Ddx58",   # RIG-I
  "Ifih1"    # MDA5
)


pdf(paste0(outdir, "/Analysis/violin.IFN.pdf"), width = 36, height = 24)

VlnPlot(
  object = data_downsampled,
  features = interferon_genes_20,
  group.by = "Clustering.Round2",
  cols = nora.colors,
  pt.size = 0, split.by = "group"
)
dev.off()




## DEG


Idents(data) <- "Clustering.Round2"


markers <- FindAllMarkers(
  object = data,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

head(markers)


top30_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 30, with_ties = FALSE) %>%
  ungroup()


as.data.frame(top30_per_cluster)

write.table(markers, paste0(rsdir,"table.clusters.markers.filt.Clusterized.R2.tsv"), sep='\t')

write.table(markers, paste0(rsdir,"table.clusters.markers.Clusterized.R2.tsv"), sep='\t')


write.table(top30_per_cluster, paste0(rsdir,"table.clusters.markers.Clusterized.top30.R2.tsv"), sep='\t')



### Bubble plot for markers


clusters<- read.table(paste0(rsdir,"table.clusters.markers.Clusterized.top30.R2.tsv"), sep='\t', header=T)


# Epsilon to avoid zero multiplication
epsilon <- 1e-8

# New column 'score'
clusters <- clusters %>%
  mutate(score = avg_log2FC * abs(pct.1 - pct.2 + epsilon))


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
    "MHCII|Ccl12 Mac",
    "Neutrophils"
)


clusters_macro <- clusters %>%
  filter(cluster %in% macro_clusters)

clusters_otros <- clusters %>%
  filter(!cluster %in% macro_clusters)


macros <- subset(data, subset = Clustering.Round2 %in% macro_clusters)

others <- subset(data, subset = Clustering.Round2 %in% macro_clusters, invert = TRUE)





# For 6 genes per cluster
top_markers_subset <- clusters_macro %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = 6) %>%
  pull(gene)


top_markers_subset <- unique(top_markers_subset)


# For 8 genes per cluster
top_markers_subset2 <- clusters_otros  %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 8) %>%
  pull(gene)


cluster_order <- c(
  "Ly6cHi Monocytes", 
  "Ly6cLo Monocytes", 
  "Early IFN|MHCII-TAMs",
  "Trem1|Ptgs2|Plaur|Celc4e Mac",
  "Mrc1|C1qc|Cbr2|Gas6 Mac",
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac",
  "Npr2|Actn1 Mac",
  "Mmp9|Ctsk Mac",
  "IFN Mac",
  "Fn1|Vegfa Mac",
  "MHCII|Siglec Mac",
  "MHCII|Ccl12 Mac",
  "Neutrophils"
)




nora.colors2 <- c(
  "Ly6cHi Monocytes"     = "#FF0000",  # Rojo inflamatorio puro
  "Ly6cLo Monocytes"     = "#FF6A6A",  # Rojo salmón transición
  "Early IFN|MHCII-TAMs"       = "#B22222", 
  "Trem1|Ptgs2|Plaur|Celc4e Mac"   = "#EE7942",
  "Mrc1|C1qc|Cbr2|Gas6 Mac"        = "#FFD92F",
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac" = "#4DAF4A",
  "Npr2|Actn1 Mac"                 = "#A6D854",
  "Neutrophils"                    = "#4876FF",
  "Mmp9|Ctsk Mac"                  = "#00723F",
  "IFN Mac"                        = "#C080FF",
  "Fn1|Vegfa Mac"                  = "#FFA500",
  "MHCII|Siglec Mac"               = "#1E90FF",
  "MHCII|Ccl12 Mac"                = "#4682B4"
)


macros$Clustering.Round2 <- factor(macros$Clustering.Round2)
macros$Clustering.Round2 <- droplevels(macros$Clustering.Round2)
table(macros$Clustering.Round2)


macros$Clustering.Round2 <- factor(macros$Clustering.Round2, levels = cluster_order)
Idents(macros) <- "Clustering.Round2"


# --- Plots---

pdf(paste0(outdir, "/Analysis/bubble.macros.pdf"), width = 12, height = 14)
p <- Clustered_DotPlot(
  seurat_object = macros,
  features = top_markers_subset,
  colors_use_idents = nora.colors2
)

print(p[[1]])
dev.off()


pdf(paste0(outdir, "/Analysis/bubble.macros.pdf"), width = 12, height = 14)

p <- DotPlot(
  macros,
  features = unique(intersect(top_markers_subset, rownames(macros))),
  cols = c("lightgrey", "darkred")
) + 
  scale_color_manual(values = nora.colors2) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p)
dev.off()




# --- Plots ---
pdf(paste0(outdir, "/Analysis/bubble.macros.k5.pdf"), width=12, height=14)

p <- Clustered_DotPlot(seurat_object = macros, features = top_markers_subset,
colors_use_idents = nora.colors2)


print(p[[1]])
dev.off()



nora.colors3 <- nora.colors[!names(nora.colors) %in% names(nora.colors2)]


# 1️⃣ Identificar los clusters con >0 células
clusters_present <- names(table(others$Cluster))[table(others$Cluster) > 0]

# 2️⃣ Crear nora.colors3 con solo esos clusters
nora.colors3 <- nora.colors[clusters_present]

# 3️⃣ Revisar
nora.colors3




Idents(others) <- factor(Idents(others), levels = names(nora.colors3))

# --- Graficar y exportar ---
pdf(paste0(outdir, "/Analysis/bubble.other.pdf"), width=10, height=10)

p <- Clustered_DotPlot(seurat_object = others, features = top_markers_subset2,
colors_use_idents = nora.colors3)


print(p[[1]])
dev.off()


# --- Graficar y exportar ---
pdf(paste0(outdir, "/Analysis/bubble.other.k4.pdf"), width=12, height=14)

p <- Clustered_DotPlot(seurat_object = others, features = top_markers_subset2,
colors_use_idents = nora.colors3, k=4)


print(p[[1]])
dev.off()




# Filtrar nora.colors para que solo tenga los clusters presentes en 'macros'
nora.colors.filtered <- nora.colors[names(nora.colors) %in% levels(Idents(macros))]

# Asegurarnos de que los niveles del objeto Seurat estén en el mismo orden
Idents(macros) <- factor(Idents(macros), levels = names(nora.colors.filtered))


pdf(paste0(outdir, "/Analysis/bubble.macros.pdf"), width=12, height=14)

p <- Clustered_DotPlot(seurat_object = macros, features = top_markers_subset)

# Alinear colores con los niveles actuales del objeto
p[[1]] <- p[[1]] +
  scale_color_manual(values = nora.colors[levels(macros)]) +
  scale_fill_manual(values = nora.colors[levels(macros)])

print(p[[1]])
dev.off()




# 2️⃣ Seleccionar los 8 principales genes por cluster
top8_markers <- clusters %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 6) %>%
  ungroup()

# Eliminar duplicados de la lista de genes
top8_genes <- unique(top8_markers$gene)

# 3️⃣ Crear el dotplot sin duplicados
dotplot <- DotPlot(
  object = macros,
  features = top_markers_subset
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




pdf(paste0(outdir,"/Analysis/bubble1.pdf"), width=36, height=12)
print(dotplot)
dev.off()





# Quantification


# Proportions por tag
prop_df <- data_downsampled@meta.data %>%
  group_by(tag, Clustering.Round2) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(tag) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.Round2
prop_df$Clustering.Round2 <- factor(
  prop_df$Clustering.Round2,
  levels = names(nora.colors)
)

# Stacked barplot
pdf(paste0(outdir,"/Analysis/quantification.tags.pdf"), width=24, height=12)
ggplot(prop_df, aes(x = tag, y = prop, fill = Clustering.Round2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(x = "Tag", y = "Cell proportion (Downsampled to 8k)", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





# Crear la columna con 4 grupos combinados (robusto a orden distinto en 'tag')
data_downsampled$group4 <- case_when(
  str_detect(data_downsampled$tag, "WT")  & str_detect(data_downsampled$tag, "DsRedN") ~ "WT_DsRedN",
  str_detect(data_downsampled$tag, "WT")  & str_detect(data_downsampled$tag, "DsRedP") ~ "WT_DsRedP",
  str_detect(data_downsampled$tag, "KO")  & str_detect(data_downsampled$tag, "DsRedN") ~ "KO_DsRedN",
  str_detect(data_downsampled$tag, "KO")  & str_detect(data_downsampled$tag, "DsRedP") ~ "KO_DsRedP",
  TRUE ~ as.character(data_downsampled$tag)    # por si acaso alguna etiqueta rara
)



table(data_downsampled$group4)



prop_df <- data_downsampled@meta.data %>%
  group_by(group4, Clustering.Round2) %>%          # usa la columna Cluster que ya tenías
  summarise(n = n(), .groups = "drop") %>%
  group_by(group4) %>%
  mutate(prop = n / sum(n))

# Reorder levels 
prop_df$group4 <- factor(prop_df$group4,
                         levels = c("WT_DsRedN","WT_DsRedP","KO_DsRedN","KO_DsRedP"))


prop_df$Clustering.Round2 <- factor(
  prop_df$Clustering.Round2,
  levels = names(nora.colors)
)

pdf(paste0(outdir,"/Analysis/quantification.group4.pdf"), width = 16, height = 12)
ggplot(prop_df, aes(x = group4, y = prop, fill = Clustering.Round2)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(x = "Group", y = "Cell proportion (Downsampled to 8K)", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()





prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round2",
	sample_1 = "KO_DsRedN", sample_2 = "KO_DsRedP",
	sample_identity = "group4"
)


pdf(paste0(outdir,"/Analysis/Proportion.test.KO.pdf"), width=8, height=10)
plot1 <- permutation_plot(prop_test) + ggtitle("Proportion test.DsRedP vs DsRedN. KO samples")
plot1
dev.off()





prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round2",
	sample_1 = "WT_DsRedN", sample_2 = "WT_DsRedP",
	sample_identity = "group4"
)


pdf(paste0(outdir,"/Analysis/Proportion.test.WT.pdf"), width=8, height=10)
plot2 <- permutation_plot(prop_test) + ggtitle("Proportion test.DsRedP vs DsRedN. WT samples")
plot2
dev.off()




prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round2",
	sample_1 = "WT_DsRedP", sample_2 = "KO_DsRedP",
	sample_identity = "group4"
)


pdf(paste0(outdir,"/Analysis/Proportion.test.Positive.KOvsWT.pdf"), width=8, height=10)
plot3 <- permutation_plot(prop_test) + ggtitle("Proportion test.KO vs WT. DsRedP samples")
plot3
dev.off()



prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round2",
	sample_1 = "WT_DsRedN", sample_2 = "KO_DsRedN",
	sample_identity = "group4"
)


pdf(paste0(outdir,"/Analysis/Proportion.test.Negative.KOvsWT.pdf"), width=8, height=10)
plot4 <- permutation_plot(prop_test) + ggtitle("Proportion test.KO vs WT. DsRedN samples")
plot4
dev.off()



pdf(paste0(outdir,"/Analysis/Proportion.tests.all.pdf"), width=30, height=8)
plot_grid(plot1, plot2, plot3, plot4, ncol=4)
dev.off() 



# DEG for macros




# Change parameters to run the proper DEG analysis
Idents(data) <- "Clustering.Round2"

# Encontrar marcadores
markers <- FindAllMarkers(
  object = data,
  only.pos = FALSE,
  min.pct = 0.1,
  logfc.threshold = 0.25
)

head(markers)


top30_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 30, with_ties = FALSE) %>%
  ungroup()


as.data.frame(top30_per_cluster)

write.table(markers, paste0(rsdir,"table.clusters.markers.filt.Clusterized.R2.tsv"), sep='\t')

write.table(markers, paste0(rsdir,"table.clusters.markers.Clusterized.R2.tsv"), sep='\t')


write.table(top30_per_cluster, paste0(rsdir,"table.clusters.markers.Clusterized.top30.R2.tsv"), sep='\t')


clusters<- read.table(paste0(rsdir,"table.clusters.markers.Clusterized.R2.tsv"), sep='\t', header=T)




library(dplyr)
library(ggplot2)
library(ggrepel)

# --- Preparar los datos ---
markers_macro <- clusters %>%
  filter(cluster %in% macro_clusters) %>%        # filtrar clusters de interés
  mutate(
    diffexpressed = ifelse(avg_log2FC > 0, "Up", "Down"),
    cluster = factor(cluster, levels = macro_clusters)  # mantener orden
  )

# --- Seleccionar primeros 6 UP y 6 DOWN según orden de la tabla ---
top_genes <- markers_macro %>%
  group_by(cluster, diffexpressed) %>%
  slice_head(n = 10) %>%   # primeros 6 genes de cada tipo
  ungroup() %>%
  mutate(color = ifelse(diffexpressed == "Up", "red", "blue"))

# --- Crear dataframe para recuadros de color centrados en y = 0 ---
cluster_boxes <- data.frame(
  cluster = macro_clusters,
  xmin = seq_along(macro_clusters) - 0.4,
  xmax = seq_along(macro_clusters) + 0.4,
  ymin = -0.15,
  ymax = 0.15
)
cluster_boxes$cluster <- factor(cluster_boxes$cluster, levels = macro_clusters)
cluster_boxes$fill <- nora.colors[as.character(cluster_boxes$cluster)]



# --- Definir epsilon por si acaso hay NA ---
epsilon <- 1e-8

# --- Preparar los datos ---
markers_macro <- clusters %>%
  filter(cluster %in% macro_clusters) %>%        # filtrar clusters de interés
  mutate(
    # NUEVA MÉTRICA: logFC ponderado por especificidad (% diferencia)
    metric = avg_log2FC * abs(pct.1 - pct.2 + epsilon),
    diffexpressed = ifelse(avg_log2FC > 0, "Up", "Down"),
    cluster = factor(cluster, levels = macro_clusters)  # mantener orden
  )

# --- Seleccionar primeros 10 UP y 10 DOWN según orden original ---
top_genes <- markers_macro %>%
  group_by(cluster, diffexpressed) %>%
  slice_max(order_by = abs(metric), n = 10) %>%
  ungroup() %>%
  mutate(color = ifelse(diffexpressed == "Up", "red", "blue"))

# --- Crear dataframe para recuadros de color centrados en y = 0 ---
cluster_boxes <- data.frame(
  cluster = macro_clusters,
  xmin = seq_along(macro_clusters) - 0.4,
  xmax = seq_along(macro_clusters) + 0.4,
  ymin = -0.15,
  ymax = 0.15
)
cluster_boxes$cluster <- factor(cluster_boxes$cluster, levels = macro_clusters)
cluster_boxes$fill <- nora.colors[as.character(cluster_boxes$cluster)]

# --- Plot principal ---
yvar <- "metric"   # usamos la nueva métrica\

pdf(paste0(outdir, "/Analysis/Top_genes_metric_pctdiff_clean.pdf"), width = 22, height = 14)

ggplot(markers_macro, aes(x = cluster, y = !!sym(yvar))) +
  geom_rect(data = cluster_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
            inherit.aes = FALSE, alpha = 0.8) +
  geom_jitter(width = 0.25, alpha = 0.4, color = "grey") +
  geom_point(
  data = top_genes,
  aes(x = cluster, y = !!sym(yvar), color = color),
  size = 2,
  position = position_jitter(width = 0.25)  # <--- aquí está la clave
) +
  geom_text_repel(
    data = top_genes,
    aes(x = cluster, y = !!sym(yvar), label = gene),
    size = 3.5,
    max.overlaps = 60,
    force = 5,
    force_pull = 0.5,
    segment.size = 0.3,
    position = position_jitter(width = 0.25)  # dispersa labels horizontalmente
  ) +
  scale_color_identity(name = "Regulation", labels = c("Up", "Down"), guide = "legend") +
  scale_fill_manual(name = "Cluster",
                    values = nora.colors,
                    breaks = macro_clusters) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right"
  ) +
  labs(
    x = NULL,
    y = "Differential Expression Score",
    title = "Top differentially expressed genes per cluster (cluster vs. all others, weighted by expression specificity)"
  )

dev.off()






