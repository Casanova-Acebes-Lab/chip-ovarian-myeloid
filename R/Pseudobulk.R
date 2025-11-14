library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(ggrepel)
library(Matrix)
library(limma)
library(edgeR)
library(Matrix.utils)
library(dplyr)





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




### Pseudobulk


data@meta.data$sample <- data@meta.data$tag
table(data@meta.data$sample)



# --------------------------------------------
# Lets go
# --------------------------------------------




# Ly6cHi_Monocytes
Ly6cHi_Monocytes <- run_pseudobulk_analysis(data, "Ly6cHi Monocytes", outdir, rsdir)

Ly6cHi_Monocytes[[1]]
head(Ly6cHi_Monocytes[[2]], n=50)
Ly6cHi_Monocytes[[3]]
Ly6cHi_Monocytes[[4]]


# Ly6cLo_Monocytes
Ly6cLo_Monocytes <- run_pseudobulk_analysis(data, "Ly6cLo Monocytes", outdir, rsdir)

Ly6cLo_Monocytes[[1]]
head(Ly6cLo_Monocytes[[2]], n=50)
Ly6cLo_Monocytes[[3]]
Ly6cLo_Monocytes[[4]]



# Early IFN|MHCII-TAMs

Early_IFN_MHCII_TAMs <- run_pseudobulk_analysis(data, "Early IFN|MHCII-TAMs", outdir, rsdir)  
Early_IFN_MHCII_TAMs[[1]]       
head(Early_IFN_MHCII_TAMs[[2]], n=50)
Early_IFN_MHCII_TAMs[[3]]
Early_IFN_MHCII_TAMs[[4]]

# Arg1|Spp1|Mmp12|Mmp19|Il1a Mac
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac <- run_pseudobulk_analysis(data, "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac", outdir, rsdir)
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[1]]
head(Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[2]], n=50)
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[3]]
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[4]]

# Trem1|Ptgs2|Plaur|Celc4e Mac
Trem1_Ptgs2_Plaur_Celc4e_Mac <- run_pseudobulk_analysis(data, "Trem1|Ptgs2|Plaur|Celc4e Mac", outdir, rsdir)    
Trem1_Ptgs2_Plaur_Celc4e_Mac[[1]]
head(Trem1_Ptgs2_Plaur_Celc4e_Mac[[2]], n=50)
Trem1_Ptgs2_Plaur_Celc4e_Mac[[3]]
Trem1_Ptgs2_Plaur_Celc4e_Mac[[4]]

# MHCII|Ccl12 Mac
MHCII_Ccl12_Mac <- run_pseudobulk_analysis(data, "MHCII|Ccl12 Mac", outdir, rsdir)
MHCII_Ccl12_Mac[[1]]
head(MHCII_Ccl12_Mac[[2]], n=50)
MHCII_Ccl12_Mac[[3]]
MHCII_Ccl12_Mac[[4]]


# MHCII|Siglec Mac
MHCII_Siglec_Mac <- run_pseudobulk_analysis(data, "MHCII|Siglec Mac", outdir, rsdir)
MHCII_Siglec_Mac[[1]]
head(MHCII_Siglec_Mac[[2]], n=50)
MHCII_Siglec_Mac[[3]]
MHCII_Siglec_Mac[[4]]

# IFN Mac
IFN_Mac <- run_pseudobulk_analysis(data, "IFN Mac", outdir, rsdir)
IFN_Mac[[1]]
head(IFN_Mac[[2]], n=50)
IFN_Mac[[3]]
IFN_Mac[[4]]

# Mmp9|Ctsk Mac
Mmp9_Ctsk_Mac <- run_pseudobulk_analysis(data, "Mmp9|Ctsk Mac", outdir, rsdir)
Mmp9_Ctsk_Mac[[1]]
head(Mmp9_Ctsk_Mac[[2]], n=50)
Mmp9_Ctsk_Mac[[3]]
Mmp9_Ctsk_Mac[[4]]

# Mrc1|C1qc|Cbr2|Gas6 Mac
Mrc1_C1qc_Cbr2_Gas6_Mac <- run_pseudobulk_analysis(data, "Mrc1|C1qc|Cbr2|Gas6 Mac", outdir, rsdir)
Mrc1_C1qc_Cbr2_Gas6_Mac[[1]]
head(Mrc1_C1qc_Cbr2_Gas6_Mac[[2]], n=50)
Mrc1_C1qc_Cbr2_Gas6_Mac[[3]]
Mrc1_C1qc_Cbr2_Gas6_Mac[[4]]

# Npr2|Actn1 Mac
Npr2_Actn1_Mac <- run_pseudobulk_analysis(data, "Npr2|Actn1 Mac", outdir, rsdir)
Npr2_Actn1_Mac[[1]]
head(Npr2_Actn1_Mac[[2]], n=50)
Npr2_Actn1_Mac[[3]]
Npr2_Actn1_Mac[[4]]

# Fn1|Vegfa Mac
Fn1_Vegfa_Mac <- run_pseudobulk_analysis(data, "Fn1|Vegfa Mac", outdir, rsdir)
Fn1_Vegfa_Mac[[1]]
head(Fn1_Vegfa_Mac[[2]], n=50)
Fn1_Vegfa_Mac[[3]]
Fn1_Vegfa_Mac[[4]]

# Neutrophils
Neutrophils <- run_pseudobulk_analysis(data, "Neutrophils", outdir, rsdir)
Neutrophils[[1]]
head(Neutrophils[[2]], n=50)
Neutrophils[[3]]
Neutrophils[[4]]



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
# 🚀 For all macróphages clusters
# ---------------------------


# Cluster list
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

pb_out <- pseudobulk_group(data, clusters = macro_clusters, cluster_col="Clustering.Round2", sample_col="tag")
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


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.pdf"), width=16, height=12)
print(p)
dev.off()



### Plot for DEG in all clusters




Ly6cHi_Monocytes[[2]]$cluster <- "Ly6cHi Monocytes"
Ly6cLo_Monocytes[[2]]$cluster <- "Ly6cLo Monocytes"
Early_IFN_MHCII_TAMs[[2]]$cluster <- "Early_IFN_MHCII_TAMs"
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[2]]$cluster <- "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac"
Trem1_Ptgs2_Plaur_Celc4e_Mac[[2]]$cluster <- "Trem1|Ptgs2|Plaur|Celc4e Mac"
MHCII_Ccl12_Mac[[2]]$cluster <- "MHCII|Ccl12 Mac"
MHCII_Siglec_Mac[[2]]$cluster <- "MHCII|Siglec Mac"
IFN_Mac[[2]]$cluster <- "IFN Mac"
Mmp9_Ctsk_Mac[[2]]$cluster <- "Mmp9|Ctsk Mac"
Mrc1_C1qc_Cbr2_Gas6_Mac[[2]]$cluster <- "Mrc1|C1qc|Cbr2|Gas6 Mac"
Npr2_Actn1_Mac[[2]]$cluster <- "Npr2|Actn1 Mac"
Fn1_Vegfa_Mac[[2]]$cluster <- "Fn1|Vegfa Mac"
Neutrophils[[2]]$cluster <- "Neutrophils"


# Join all dataframes
all_clusters <- bind_rows(
Ly6cHi_Monocytes[[2]],
Ly6cLo_Monocytes[[2]],
Early_IFN_MHCII_TAMs[[2]],
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[2]],
Trem1_Ptgs2_Plaur_Celc4e_Mac[[2]],
MHCII_Ccl12_Mac[[2]],
MHCII_Siglec_Mac[[2]],
IFN_Mac[[2]],
Mmp9_Ctsk_Mac[[2]],
Mrc1_C1qc_Cbr2_Gas6_Mac[[2]],
Npr2_Actn1_Mac[[2]],
Fn1_Vegfa_Mac[[2]],
Neutrophils[[2]],
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
  "Cxcl1","C1qc","Hexim1","Lgals1","Marcksl1","Nfkbia","Oas2","Ahnak","Cxcl9","Egr2",
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
pdf(paste0(outdir, "/QC/Shared_genes_metric_plot_rects_final_labels.pdf"), width = 22, height = 14)

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
library(ggplot2)
library(ggrepel)

# --- Definir clusters de interés ---
macro_clusters <- unique(all_clusters$cluster)  # o tu lista de clusters seleccionados
n_top <- 10  # top N Up/Down por cluster

# --- Preparar los datos ---
all_clusters <- all_clusters %>%
  mutate(
    diffexpressed = ifelse(metric > 0, "Up", "Down"),
    cluster = factor(cluster, levels = macro_clusters)
  )

# --- Seleccionar top 10 Up y 10 Down por cluster ---
top_genes <- all_clusters %>%
  filter(diffexpressed %in% c("Up","Down")) %>%
  group_by(cluster, diffexpressed) %>%
  slice_max(abs(metric), n = n_top) %>%
  ungroup() %>%
  mutate(color = ifelse(diffexpressed == "Up", "#FF6666", "#6699FF"))

# --- Crear recuadros de color centrados en y = 0 ---
cluster_boxes <- data.frame(
  cluster = macro_clusters,
  xmin = seq_along(macro_clusters) - 0.4,
  xmax = seq_along(macro_clusters) + 0.4,
  ymin = -0.15,
  ymax = 0.15
)
cluster_boxes$cluster <- factor(cluster_boxes$cluster, levels = macro_clusters)
cluster_boxes$fill <- nora.colors2[as.character(cluster_boxes$cluster)]

# --- Plot ---
png(paste0(outdir, "/Pseudobulk/Top_genes_per_cluster.png"), width = 2200, height = 1400)

ggplot(all_clusters, aes(x = cluster, y = metric)) +
  geom_rect(data = cluster_boxes,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
            inherit.aes = FALSE, alpha = 0.8) +
  geom_jitter(width = 0.25, alpha = 0.3, color = "grey") +
  geom_point(data = top_genes,
             aes(x = cluster, y = metric, color = color),
             size = 3,
             position = position_jitter(width = 0.25)) +
  geom_text_repel(data = top_genes,
                  aes(x = cluster, y = metric, label = gene),
                  size = 4,
                  max.overlaps = 60,
                  force = 5,
                  force_pull = 0.5,
                  segment.size = 0.3,
                  position = position_jitter(width = 0.25)) +
  scale_color_identity(name = "Regulation", guide = "legend") +
  scale_fill_manual(name = "Cluster",
                    values = nora.colors2,
                    breaks = macro_clusters) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid = element_blank(),
        legend.position = "right") +
  labs(x = NULL,
       y = "Differential Expression Score",
       title = "Top 10 Up and 10 Down genes per cluster")

dev.off()

