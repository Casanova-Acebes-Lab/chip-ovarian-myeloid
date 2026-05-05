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

data <- readRDS(paste0(rsdir,"/objects/data.Clustering.Round1.rds"))




# -----------------------------
# 0. Rutas y objeto
# -----------------------------

# Ajustar si hace falta
# data <- readRDS(file.path(outdir, "seurat_QC2_noBadClusters_ClusteringRound3_final.rds"))

pb_dir <- file.path(outdir, "Pseudobulk_Macrophages_Monocytes_Round3")
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1. Comprobar columnas necesarias
# -----------------------------

required_cols <- c(
  "sample_id",
  "batch",
  "group",
  "DsRed",
  "mouse_id",
  "tet2_status",
  "cell_context",
  "Clustering.Round3",
  "Clustering.wide"
)

missing_cols <- setdiff(required_cols, colnames(data@meta.data))

if (length(missing_cols) > 0) {
  stop("Faltan columnas en metadata: ", paste(missing_cols, collapse = ", "))
}

# -----------------------------
# 2. Seleccionar macrófagos + monocitos
# -----------------------------

macmono_obj <- subset(
  data,
  subset = Clustering.wide %in% c("Macrophages", "IFN Mac", "Monocytes")
)

table(macmono_obj$Clustering.wide)
table(macmono_obj$Clustering.Round3)
table(macmono_obj$sample_id)



make_pseudobulk_for_cluster <- function(seu, cluster_name, assay = "RNA") {
  
  message("Processing cluster: ", cluster_name)
  
  obj_cl <- subset(
    seu,
    subset = Clustering.Round3 == cluster_name
  )
  
  counts <- GetAssayData(
    obj_cl,
    assay = assay,
    layer = "counts"
  )
  
  samples <- obj_cl$sample_id
  names(samples) <- colnames(obj_cl)
  samples <- samples[colnames(counts)]
  
  stopifnot(ncol(counts) == length(samples))
  stopifnot(all(names(samples) == colnames(counts)))
  stopifnot(!anyNA(samples))
  
  pb_counts <- Matrix::t(
    Matrix.utils::aggregate.Matrix(
      x = Matrix::t(counts),
      groupings = samples,
      fun = "sum"
    )
  )
  
  # SoupX puede dejar counts decimales.
  # Redondeamos después de sumar, no célula a célula.
  pb_counts <- round(pb_counts)
  pb_counts <- as.matrix(pb_counts)
  storage.mode(pb_counts) <- "integer"
  
  pb_meta <- obj_cl@meta.data %>%
    as.data.frame() %>%
    dplyr::distinct(
      sample_id,
      batch,
      group,
      DsRed,
      mouse_id,
      tet2_status,
      cell_context
    )
  
  rownames(pb_meta) <- pb_meta$sample_id
  pb_meta <- pb_meta[colnames(pb_counts), , drop = FALSE]
  
  pb_meta$group <- factor(pb_meta$group, levels = c("WT", "KO"))
  pb_meta$DsRed <- factor(pb_meta$DsRed, levels = c("DsRedN", "DsRedP"))
  pb_meta$batch <- factor(pb_meta$batch, levels = c("A", "B", "C"))
  
  n_cells <- table(obj_cl$sample_id)
  pb_meta$n_cells_cluster <- as.numeric(n_cells[rownames(pb_meta)])
  pb_meta$cluster <- cluster_name
  
  stopifnot(all(rownames(pb_meta) == colnames(pb_counts)))
  
  list(
    counts = pb_counts,
    meta = pb_meta
  )
}



macmono_clusters <- sort(unique(as.character(macmono_obj$Clustering.Round3)))

macmono_clusters

pb_list <- lapply(
  macmono_clusters,
  function(cl) {
    make_pseudobulk_for_cluster(
      seu = macmono_obj,
      cluster_name = cl,
      assay = "RNA"
    )
  }
)

names(pb_list) <- macmono_clusters



pb_dims <- data.frame(
  cluster = names(pb_list),
  n_genes = sapply(pb_list, function(x) nrow(x$counts)),
  n_samples = sapply(pb_list, function(x) ncol(x$counts)),
  min_cells = sapply(pb_list, function(x) min(x$meta$n_cells_cluster)),
  median_cells = sapply(pb_list, function(x) median(x$meta$n_cells_cluster)),
  max_cells = sapply(pb_list, function(x) max(x$meta$n_cells_cluster))
) %>%
  dplyr::arrange(min_cells)

as.data.frame(pb_dims)




run_deseq2_global_adjusted <- function(pb) {
  
  counts <- pb$counts
  meta <- pb$meta
  
  meta$group <- droplevels(factor(meta$group, levels = c("WT", "KO")))
  meta$DsRed <- droplevels(factor(meta$DsRed, levels = c("DsRedN", "DsRedP")))
  
  stopifnot(all(rownames(meta) == colnames(counts)))
  
  # Filtro de genes: al menos 10 counts en al menos 2 pseudobulks
  keep_genes <- rowSums(counts >= 10) >= 2
  counts <- counts[keep_genes, , drop = FALSE]
  
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = meta,
    design = ~ DsRed + group
  )
  
  dds <- DESeq2::DESeq(dds)
  
  res <- DESeq2::results(
    dds,
    contrast = c("group", "KO", "WT"),
    alpha = 0.05
  )
  
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  
  res <- res %>%
    dplyr::select(gene, everything()) %>%
    dplyr::arrange(padj)
  
  list(
    dds = dds,
    res = res,
    meta = meta
  )
}


clusters_de <- names(pb_list)
clusters_de

de_list <- lapply(
  clusters_de,
  function(cl) {
    message("Running DESeq2: ", cl)
    run_deseq2_global_adjusted(pb_list[[cl]])
  }
)

names(de_list) <- clusters_de




de_summary <- dplyr::bind_rows(
  lapply(names(de_list), function(cl) {
    
    res <- de_list[[cl]]$res
    meta <- de_list[[cl]]$meta
    
    data.frame(
      cluster = cl,
      n_samples = nrow(meta),
      min_cells = min(meta$n_cells_cluster),
      median_cells = median(meta$n_cells_cluster),
      max_cells = max(meta$n_cells_cluster),
      n_genes_tested = nrow(res),
      n_sig_005 = sum(res$padj < 0.05, na.rm = TRUE),
      n_sig_010 = sum(res$padj < 0.10, na.rm = TRUE),
      n_up_005 = sum(res$padj < 0.05 & res$log2FoldChange > 0, na.rm = TRUE),
      n_down_005 = sum(res$padj < 0.05 & res$log2FoldChange < 0, na.rm = TRUE),
      top_gene = res$gene[which.min(res$padj)],
      top_padj = min(res$padj, na.rm = TRUE),
      top_log2FC = res$log2FoldChange[which.min(res$padj)]
    )
  })
) %>%
  dplyr::arrange(desc(n_sig_005))

as.data.frame(de_summary)



top_genes_by_cluster <- dplyr::bind_rows(
  lapply(names(de_list), function(cl) {
    
    de_list[[cl]]$res %>%
      dplyr::mutate(cluster = cl) %>%
      dplyr::filter(!is.na(padj)) %>%
      dplyr::arrange(padj) %>%
      dplyr::slice_head(n = 20)
  })
) %>%
  dplyr::select(cluster, gene, log2FoldChange, lfcSE, stat, pvalue, padj)

as.data.frame(top_genes_by_cluster)





de_list[["IFN Mac"]]$res %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::arrange(padj) %>%
  head(60)



  de_list[["IFN Mac"]]$res %>%
  dplyr::filter(!is.na(padj), log2FoldChange > 0) %>%
  dplyr::arrange(padj) %>%
  head(60)



## Volcanos

volc_dir <- file.path(pb_dir, "VolcanoPlots")
dir.create(volc_dir, recursive = TRUE, showWarnings = FALSE)



library(ggplot2)
library(ggrepel)
library(dplyr)

plot_volcano_cluster <- function(
  res_df,
  cluster_name,
  padj_cutoff = 0.05,
  lfc_cutoff = 0.5,
  n_label = 30
) {
  
  df <- as.data.frame(res_df)
  df$gene <- rownames(df)
  
  # Limpiar
  df <- df %>%
    filter(!is.na(log2FoldChange), !is.na(padj)) %>%
    mutate(
      neglog10padj = -log10(padj),
      regulation = case_when(
        padj < padj_cutoff & log2FoldChange >= lfc_cutoff  ~ "Up in KO",
        padj < padj_cutoff & log2FoldChange <= -lfc_cutoff ~ "Up in WT",
        TRUE ~ "NS"
      )
    )
  
  # Genes a etiquetar: top por padj entre significativos
  label_df <- df %>%
    filter(regulation != "NS") %>%
    arrange(padj) %>%
    distinct(gene, .keep_all = TRUE) %>%
    slice_head(n = n_label)
  
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10padj)) +
    geom_point(aes(color = regulation), alpha = 0.7, size = 1.3) +
    scale_color_manual(
      values = c(
        "Up in KO" = "firebrick3",
        "Up in WT" = "steelblue3",
        "NS" = "grey75"
      )
    ) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
    geom_text_repel(
      data = label_df,
      aes(label = gene),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.35,
      point.padding = 0.2
    ) +
    labs(
      title = paste0("Volcano plot: ", cluster_name),
      subtitle = paste0("KO vs WT DEGs"),
      x = "log2FoldChange (KO / WT)",
      y = "-log10 adjusted p-value",
      color = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
  
  return(p)
}


names(de_list[[1]])



volcano_list <- list()

for (cl in names(de_list)) {
  
  cat("Plotting volcano:", cl, "\n")
  
  res_df <- de_list[[cl]]$res
  
  p <- plot_volcano_cluster(
    res_df = res_df,
    cluster_name = cl,
    padj_cutoff = 0.05,
    lfc_cutoff = 0.5,   # puedes cambiar a 1 si lo quieres más estricto
    n_label = 30
  )
  
  volcano_list[[cl]] <- p
  
  file_name <- gsub("[^[:alnum:]_]+", "_", cl)
  
  ggsave(
    filename = file.path(volc_dir, paste0("Volcano_", file_name, ".png")),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
}



# -----------------------------
# Guardar resultados DESeq2
# -----------------------------

res_dir <- file.path(pb_dir, "DESeq2_KO_vs_WT_adjusted_DsRed_all_clusters")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

saveRDS(
  de_list,
  file = file.path(res_dir, "de_list_KO_vs_WT_adjusted_DsRed_all_clusters.rds")
)

saveRDS(
  pb_list,
  file = file.path(res_dir, "pb_list_macmono_by_ClusteringRound3.rds")
)

write.csv(
  de_summary,
  file = file.path(res_dir, "DESeq2_KO_vs_WT_adjusted_DsRed_all_clusters_summary.csv"),
  row.names = FALSE
)


for (cl in names(de_list)) {
  
  safe_cl <- gsub("[^A-Za-z0-9_]+", "_", cl)
  
  write.csv(
    de_list[[cl]]$res,
    file = file.path(
      res_dir,
      paste0("DESeq2_KO_vs_WT_adjDsRed_", safe_cl, ".csv")
    ),
    row.names = FALSE
  )
}


list.files(res_dir)




#### Checking cytokines



cytokine_genes <- c(
  # IL1 family / inflamación
  "Il1b", "Il1a", "Il1rn", "Il18", "Il33",
  "Nlrp3", "Pycard", "Casp1",

  # Cytokines clásicas
  "Tnf", "Il6", "Il10", "Il12a", "Il12b", "Il23a",
  "Il27", "Il27ra", "Il4", "Il13", "Il15", "Il16",

  # CSFs
  "Csf1", "Csf2", "Csf3",

  # TGF / macrophage modulation
  "Tgfb1", "Tgfb2", "Tgfb3",

  # Chemokines inflamatorias
  "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl6", "Ccl7", "Ccl8", "Ccl9",
  "Ccl12", "Ccl17", "Ccl22", "Ccl24",

  # CXCLs
  "Cxcl1", "Cxcl2", "Cxcl3", "Cxcl5", "Cxcl9", "Cxcl10", "Cxcl11",
  "Cxcl12", "Cxcl13", "Cxcl16", "Cxcr3",

  # IFN-related cytokines
  "Ifnb1", "Ifng",

  # TNF family / related
  "Lta", "Ltb", "Tnfsf10", "Tnfsf13b", "Tnfsf14"
)




library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(patchwork)

###
### Heatmap de citoquinas con barra de color por cluster
###

cytokine_lfc_clean <- dplyr::bind_rows(
  lapply(names(de_list), function(cl) {
    
    res <- de_list[[cl]]$res
    
    if (!"gene" %in% colnames(res)) {
      res$gene <- rownames(res)
    }
    
    res %>%
      dplyr::mutate(
        gene = as.character(gene),
        cluster = cl
      ) %>%
      dplyr::filter(
        !is.na(gene),
        gene %in% cytokine_genes
      ) %>%
      dplyr::select(
        cluster,
        gene,
        log2FoldChange,
        stat,
        pvalue,
        padj
      )
  })
)

# Comprobar duplicados
cytokine_lfc_clean %>%
  dplyr::count(gene, cluster) %>%
  dplyr::filter(n > 1)

# Si hay duplicados, quedarnos con el resultado más significativo
cytokine_lfc_clean <- cytokine_lfc_clean %>%
  dplyr::group_by(gene, cluster) %>%
  dplyr::arrange(is.na(padj), padj, .by_group = TRUE) %>%
  dplyr::slice_head(n = 1) %>%
  dplyr::ungroup()

# Etiquetas de significancia
cytokine_lfc_clean <- cytokine_lfc_clean %>%
  dplyr::mutate(
    sig_label = dplyr::case_when(
      is.na(padj) ~ "",
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      padj < 0.10  ~ ".",
      TRUE ~ ""
    )
  )

# -----------------------------
# Orden manual de clusters
# -----------------------------

col_order_manual <- c(
  # Monocytes
  "Ly6c2Lo Monocytes",
  "Ly6c2Hi Monocytes",
  
  # IFN / IFN-related macrophages
  "Early IFN TAMs",
  "IFN Mac",
  "Mki67 IFN Mac",
  
  # MHCII macrophages
  "MHCII|Siglec Mac",
  "MHCII|Mgl2 Mac",
  "Saa3 Mac",
  
  # Resident / remodeling macrophages
  "Gas6|Folr2 Mac",
  "Nrp2|Emp1 Mac",
  "Nlrp3|Vegfa Mac",
  
  # Inflammatory / tissue-remodeling macrophages
  "Trem1|Ptgs2|Plaur|F10 Mac",
  "Arg1|Spp1|Mmp12|Il1a Mac",
  "Mki67|Cstk|Mmp9|S100a4 Mac",
  "Mmp9|Ctsk Mac"
)

# -----------------------------
# Colores para cada cluster
# -----------------------------

cluster_colors <- c(
  "Ly6c2Hi Monocytes"             = "#FF0000",
  "Ly6c2Lo Monocytes"             = "#FF6347",
  "Early IFN TAMs"                = "#FFA07A",
  "IFN Mac"                       = "#C080FF",
  "Mki67 IFN Mac"                 = "#504369",
  "MHCII|Siglec Mac"              = "#1E90FF",
  "MHCII|Mgl2 Mac"                = "#4682B4",
  "Saa3 Mac"                      = "#FFA500",
  "Gas6|Folr2 Mac"                = "#FFD700",
  "Nrp2|Emp1 Mac"                 = "#A6D854",
  "Nlrp3|Vegfa Mac"               = "#5F3121",
  "Trem1|Ptgs2|Plaur|F10 Mac"     = "#EE7942",
  "Arg1|Spp1|Mmp12|Il1a Mac"      = "#4DAF4A",
  "Mki67|Cstk|Mmp9|S100a4 Mac"    = "#DB44A9",
  "Mmp9|Ctsk Mac"                 = "#00723F"
)

# -----------------------------
# Matriz log2FC para clustering de genes
# -----------------------------

lfc_mat <- cytokine_lfc_clean %>%
  dplyr::select(gene, cluster, log2FoldChange) %>%
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = log2FoldChange,
    values_fn = mean
  ) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Mantener solo clusters presentes y en orden manual
col_order <- col_order_manual[col_order_manual %in% colnames(lfc_mat)]
lfc_mat <- lfc_mat[, col_order, drop = FALSE]

# Para clustering: NA a 0
lfc_mat_clust <- lfc_mat
lfc_mat_clust[is.na(lfc_mat_clust)] <- 0

# Quitar genes sin variabilidad
keep_rows <- apply(lfc_mat_clust, 1, sd, na.rm = TRUE) > 0
lfc_mat_clust <- lfc_mat_clust[keep_rows, , drop = FALSE]

# -----------------------------
# Clustering de genes
# Escalado por gen + euclidean + Ward.D2
# -----------------------------

lfc_mat_scaled <- t(scale(t(lfc_mat_clust)))
lfc_mat_scaled[is.na(lfc_mat_scaled)] <- 0

row_hc <- hclust(
  dist(lfc_mat_scaled, method = "euclidean"),
  method = "ward.D2"
)

row_order <- rownames(lfc_mat_clust)[row_hc$order]

# -----------------------------
# Tabla final para plot
# -----------------------------

cytokine_lfc_plot <- cytokine_lfc_clean %>%
  dplyr::filter(
    gene %in% row_order,
    cluster %in% col_order
  ) %>%
  dplyr::mutate(
    gene = factor(gene, levels = rev(row_order)),
    cluster = factor(cluster, levels = col_order)
  )

# -----------------------------
# Plot principal del heatmap
# -----------------------------

p_heatmap <- ggplot(
  cytokine_lfc_plot,
  aes(x = cluster, y = gene, fill = log2FoldChange)
) +
  geom_tile(color = NA) +
  geom_text(
    aes(label = sig_label),
    size = 3.5,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    na.value = "grey90",
    name = "log2FC\nKOvsWT"
  ) +
  scale_x_discrete(expand = c(0, 0), drop = FALSE) +
  scale_y_discrete(expand = c(0, 0)) +
  theme_bw() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_text(size = 9),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5.5, r = 5.5, b = 0, l = 5.5)
  ) +
  labs(
    title = "Cytokine / chemokine log2FC by macrophage-monocyte cluster",
    subtitle = "KO vs WT adjusted by DsRed | genes clustered by scaled log2FC pattern",
    y = "Gene"
  )

# -----------------------------
# Barra inferior de colores + nombres de cluster
# -----------------------------

cluster_annotation_df <- data.frame(
  cluster = factor(col_order, levels = col_order),
  y = 1
)

p_cluster_bar <- ggplot(
  cluster_annotation_df,
  aes(x = cluster, y = y, fill = cluster)
) +
  geom_tile(height = 0.35) +
  scale_fill_manual(values = cluster_colors, guide = "none") +
  scale_x_discrete(expand = c(0, 0), drop = FALSE) +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
    plot.margin = margin(t = 0, r = 5.5, b = 5.5, l = 5.5)
  )

# -----------------------------
# Combinar ambos plots
# -----------------------------

p_cytokine_lfc_clustered <- p_heatmap / p_cluster_bar +
  plot_layout(heights = c(12, 0.25)) +
  plot_annotation(
    caption = "DESeq2 adjusted p-value: *** padj < 0.001; ** padj < 0.01; * padj < 0.05; . padj < 0.10"
  ) &
  theme(
    plot.caption = element_text(
      hjust = 0,
      size = 9,
      face = "italic"
    )
  )

p_cytokine_lfc_clustered

# -----------------------------
# Guardar
# -----------------------------

ggsave(
  filename = file.path(
    cytokine_plot_dir,
    "Cytokine_Chemokine_log2FC_heatmap_genesWardD2_clustersManual_KO_vs_WT_adjDsRed_withClusterColorBar.pdf"
  ),
  plot = p_cytokine_lfc_clustered,
  width = 8,
  height = 12
)

ggsave(
  filename = file.path(
    cytokine_plot_dir,
    "Cytokine_Chemokine_log2FC_heatmap_genesWardD2_clustersManual_KO_vs_WT_adjDsRed_withClusterColorBar.png"
  ),
  plot = p_cytokine_lfc_clustered,
  width = 8,
  height = 12,
  dpi = 300
)
####



genes_tcell_recruitment <- c(
  "Cxcl9", "Cxcl10", "Cxcl11",
  "Ccl5",
  "Ifnb1", "Irf7", "Stat1", "Isg15",
  "Cxcl1", "Cxcl2", "Cxcl3",
  "Il1b", "Il18", "Nlrp3", "Casp1", "Pycard",
  "Il10", "Tnf", "Tnfsf10"
)

clusters_focus <- c(
  "IFN Mac",
  "Early IFN TAMs",
  "Mki67 IFN Mac",
  "Ly6c2Hi Monocytes",
  "MHCII|Siglec Mac"
)

clusters_focus <- clusters_focus[clusters_focus %in% names(de_list)]
clusters_focus



library(dplyr)
library(tidyr)
library(tibble)
library(DESeq2)
library(ggplot2)

extract_norm_counts_cluster <- function(de_obj, cluster_name, genes_use) {
  
  dds <- de_obj$dds
  
  norm_counts <- DESeq2::counts(dds, normalized = TRUE)
  
  genes_present <- genes_use[genes_use %in% rownames(norm_counts)]
  
  if (length(genes_present) == 0) {
    return(NULL)
  }
  
  meta <- as.data.frame(SummarizedExperiment::colData(dds))
  meta$sample_id <- rownames(meta)
  
  norm_df <- norm_counts[genes_present, , drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(
      cols = -gene,
      names_to = "sample_id",
      values_to = "norm_count"
    ) %>%
    dplyr::left_join(meta, by = "sample_id") %>%
    dplyr::mutate(
      cluster = cluster_name,
      log2_norm_count = log2(norm_count + 1)
    )
  
  norm_df
}



norm_cytokines <- dplyr::bind_rows(
  lapply(clusters_focus, function(cl) {
    extract_norm_counts_cluster(
      de_obj = de_list[[cl]],
      cluster_name = cl,
      genes_use = genes_tcell_recruitment
    )
  })
)

as.data.frame(norm_cytokines)



norm_plot_dir <- file.path(pb_dir, "Normalized_pseudobulk_cytokines")
dir.create(norm_plot_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  norm_cytokines,
  file = file.path(norm_plot_dir, "normalized_pseudobulk_cytokines_by_sample.csv"),
  row.names = FALSE
)



sample_order <- c(
  "WT1_DsRedN", "WT1_DsRedP",
  "WT2_DsRedN", "WT2_DsRedP",
  "KO2_DsRedN", "KO2_DsRedP",
  "KO3_DsRedN", "KO3_DsRedP"
)

sample_order <- sample_order[sample_order %in% norm_cytokines$sample_id]

norm_cytokines <- norm_cytokines %>%
  dplyr::mutate(
    sample_id = factor(sample_id, levels = sample_order),
    group = factor(group, levels = c("WT", "KO")),
    DsRed = factor(DsRed, levels = c("DsRedN", "DsRedP"))
  )



  genes_plot_core <- c(
  "Cxcl9", "Cxcl10", "Cxcl11",
  "Ccl5",
  "Ifnb1", "Irf7", "Stat1", "Isg15",
  "Il1b", "Il18", "Nlrp3"
)

genes_plot_core <- genes_plot_core[genes_plot_core %in% norm_cytokines$gene]

p_norm_core <- norm_cytokines %>%
  dplyr::filter(gene %in% genes_plot_core) %>%
  ggplot(
    aes(
      x = sample_id,
      y = log2_norm_count,
      color = group,
      shape = DsRed
    )
  ) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(mouse_id, gene)), alpha = 0.25) +
  facet_grid(gene ~ cluster, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Normalized pseudobulk expression of T-cell recruitment / inflammatory genes",
    subtitle = "DESeq2 normalized counts | each point = one pseudobulk sample",
    x = "",
    y = "log2(normalized count + 1)",
    color = "Group",
    shape = "DsRed"
  )

p_norm_core



ggsave(
  filename = file.path(norm_plot_dir, "Normalized_pseudobulk_core_cytokines_by_sample.pdf"),
  plot = p_norm_core,
  width = 16,
  height = 12
)

ggsave(
  filename = file.path(norm_plot_dir, "Normalized_pseudobulk_core_cytokines_by_sample.png"),
  plot = p_norm_core,
  width = 16,
  height = 12,
  dpi = 300
)




genes_cxcl_axis <- c(
  "Cxcl9", "Cxcl10", "Cxcl11",
  "Cxcl1", "Cxcl2", "Cxcl3",
  "Ccl5"
)

genes_cxcl_axis <- genes_cxcl_axis[genes_cxcl_axis %in% norm_cytokines$gene]

p_cxcl_axis <- norm_cytokines %>%
  dplyr::filter(gene %in% genes_cxcl_axis) %>%
  ggplot(
    aes(
      x = sample_id,
      y = log2_norm_count,
      color = group,
      shape = DsRed
    )
  ) +
  geom_point(size = 3) +
  geom_line(aes(group = interaction(mouse_id, gene)), alpha = 0.25) +
  facet_grid(gene ~ cluster, scales = "free_y") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text.x = element_text(face = "bold"),
    strip.text.y = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "CXCL / T-cell recruitment axis in pseudobulk",
    subtitle = "KO vs WT visual check using DESeq2 normalized counts",
    x = "",
    y = "log2(normalized count + 1)",
    color = "Group",
    shape = "DsRed"
  )

p_cxcl_axis



ggsave(
  filename = file.path(norm_plot_dir, "Normalized_pseudobulk_CXCL_axis_by_sample.pdf"),
  plot = p_cxcl_axis,
  width = 16,
  height = 9
)

ggsave(
  filename = file.path(norm_plot_dir, "Normalized_pseudobulk_CXCL_axis_by_sample.png"),
  plot = p_cxcl_axis,
  width = 16,
  height = 9,
  dpi = 300
)



res_dsred_ifn <- results(
  de_list[["IFN Mac"]]$dds,
  name = "DsRed_DsRedP_vs_DsRedN"
)

as.data.frame(res_dsred_ifn[c("Cxcl1", "Cxcl2", "Cxcl3", "Cxcl9", "Cxcl10"), ])



res_group_ifn <- results(
  de_list[["IFN Mac"]]$dds,
  name = "group_KO_vs_WT"
)

as.data.frame(res_group_ifn[c("Cxcl1", "Cxcl2", "Cxcl3", "Cxcl9", "Cxcl10"), ])





### Checking Angiogenesis-related genes


# -----------------------------
# 1. Genes de angiogénesis
# -----------------------------
angiogenesis_genes <- unique(c(
  # VEGF / VEGFR / neuropilins
  "Vegfa", "Vegfb", "Vegfc", "Vegfd",
  "Pgf",
  "Flt1", "Kdr", "Flt4",
  "Nrp1", "Nrp2",
  
  # ANGPT / Tie axis
  "Angpt1", "Angpt2", "Angptl4",
  "Tek", "Tie1",
  
  # PDGF / FGF
  "Pdgfa", "Pdgfb", "Pdgfc", "Pdgfd",
  "Pdgfra", "Pdgfrb",
  "Fgf1", "Fgf2", "Fgfr1",
  
  # Hypoxia / metabolic angiogenesis
  "Hif1a", "Epas1",
  "Adm", "Nos2",
  
  # Notch / aberrant sprouting
  "Dll4", "Jag1", "Jag2",
  "Notch1", "Notch2", "Notch4",
  "Hes1", "Hey1", "Hey2",
  "Apln", "Aplnr",
  "Esm1",
  
  # Growth factors / epithelial-endothelial interaction
  "Hbegf", "Areg", "Ereg",
  "Igf1",
  "Osm",
  
  # Matrix remodeling / proteases
  "Mmp2", "Mmp9", "Mmp12", "Mmp14", "Mmp19",
  "Timp1", "Timp2", "Timp3",
  "Ctsb", "Ctsk", "Ctsl",
  "Adam8", "Adamts1",
  
  # ECM / fibrosis / physical barrier
  "Fn1",
  "Col1a1", "Col1a2",
  "Col3a1",
  "Col4a1", "Col4a2",
  "Col18a1",
  "Lox", "Loxl2",
  "Postn",
  "Sparc",
  "Serpine1",
  "Plau", "Plaur",
  
  # Antiangiogenic / incomplete angiogenesis
  "Thbs1", "Thbs2",
  "Serpinf1",
  
  # TGF-beta / stromal vascular remodeling
  "Tgfb1", "Tgfb2", "Tgfb3",
  "Tgfbr1", "Tgfbr2",
  
  # Inflammatory angiogenesis
  "Il1b", "Il1a", "Il6", "Tnf",
  "Il10",
  "Cxcl12", "Ccl2",
  
  # TAM vascular/remodeling-associated
  "Spp1",
  "Mrc1", "Folr2", "Lyve1",
  "Stab1",
  "Lgals3",
  
  # Lymphatic / alternative vascular remodeling
  "Pdpn", "Prox1",
  
  # Endothelial activation / adhesion, useful as vascular interaction readout
  "Icam1", "Vcam1", "Sele", "Selp"
))



ecm_remodeling_genes <- c(
  # Fibrillar collagens / stromal ECM
  "Col1a1", "Col1a2",
  "Col3a1",
  "Col5a1", "Col5a2",
  "Col6a1", "Col6a2", "Col6a3",
  
  # Basement membrane / vascular ECM
  "Col4a1", "Col4a2", "Col4a3", "Col4a4",
  "Col18a1",
  "Lama2", "Lama4", "Lama5",
  "Lamb1", "Lamb2",
  "Lamc1",
  "Nid1", "Nid2",
  "Hspg2",
  
  # Glycoproteins / matricellular genes
  "Fn1",
  "Tnc",
  "Sparc",
  "Postn",
  "Thbs1", "Thbs2",
  "Vcan",
  "Dcn",
  "Lum",
  "Bgn",
  "Ccn1", "Ccn2", "Ccn3",
  
  # Crosslinking / stiffness
  "Lox", "Loxl1", "Loxl2", "Loxl3",
  
  # Proteases / matrix remodeling
  "Mmp2", "Mmp3", "Mmp8", "Mmp9", "Mmp12", "Mmp13", "Mmp14", "Mmp19",
  "Adam8", "Adam9", "Adam10", "Adam12", "Adam15", "Adam17",
  "Adamts1", "Adamts2", "Adamts4", "Adamts5", "Adamts9",
  
  # MMP inhibitors
  "Timp1", "Timp2", "Timp3",
  
  # Cathepsins, macrophage matrix degradation
  "Ctsb", "Ctsk", "Ctsl", "Ctss", "Ctsd",
  
  # Plasminogen / invasion / fibrinolysis
  "Plau", "Plaur",
  "Serpine1",
  "Serpinf1",
  
  # Integrins / adhesion to matrix
  "Itga4", "Itga5", "Itgav", "Itgax", "Itgam",
  "Itgb1", "Itgb2", "Itgb3",
  
  # Adhesion / endothelial interaction
  "Icam1", "Vcam1",
  "Sele", "Selp",
  
  # Macrophage remodeling-associated
  "Spp1",
  "Lgals1", "Lgals3",
  "Tgm2",
  "Mrc1",
  "Folr2",
  "Stab1",
  "Lyve1",
  "Pros1",
  "Gas6",
  "Dcn",
  "Lum",
  "Sparc",

  
  # TGF-beta / fibrotic remodeling
  "Tgfb1", "Tgfb2", "Tgfb3",
  "Tgfbr1", "Tgfbr2",
  "Smad2", "Smad3", "Smad4",
  
  # Hypoxia / vascular remodeling
  "Hif1a", "Epas1",
  "Vegfa", "Vegfc", "Pgf",
  "Flt1", "Flt4", "Nrp1", "Nrp2"
)



angiogenesis_genes <- unique(c(
  angiogenesis_genes,
  ecm_remodeling_genes
))


# -----------------------------
# 2. Extraer resultados por cluster
# -----------------------------

angiogenesis_lfc <- lapply(names(de_list), function(cl) {
  
  res_df <- as.data.frame(de_list[[cl]]$res)
  
  if (!"gene" %in% colnames(res_df)) {
    res_df$gene <- rownames(res_df)
  }
  
  res_df %>%
    dplyr::filter(gene %in% angiogenesis_genes) %>%
    dplyr::mutate(
      cluster = cl,
      sig_label = dplyr::case_when(
        !is.na(padj) & padj < 0.001 ~ "***",
        !is.na(padj) & padj < 0.01  ~ "**",
        !is.na(padj) & padj < 0.05  ~ "*",
        !is.na(padj) & padj < 0.10  ~ ".",
        TRUE ~ ""
      )
    ) %>%
    dplyr::select(
      cluster, gene,
      log2FoldChange, stat, pvalue, padj, sig_label
    )
})

angiogenesis_lfc <- dplyr::bind_rows(angiogenesis_lfc)

as.data.frame(angiogenesis_lfc)



# -----------------------------
# 3. Limpiar duplicados
# -----------------------------

angiogenesis_lfc_clean <- angiogenesis_lfc %>%
  dplyr::filter(!is.na(gene), gene != "") %>%
  dplyr::group_by(cluster, gene) %>%
  dplyr::summarise(
    log2FoldChange = log2FoldChange[which.max(abs(stat))],
    stat = stat[which.max(abs(stat))],
    pvalue = pvalue[which.max(abs(stat))],
    padj = padj[which.max(abs(stat))],
    sig_label = sig_label[which.max(abs(stat))],
    .groups = "drop"
  )

as.data.frame(angiogenesis_lfc_clean)



# -----------------------------
# 4. Construir matrices
# -----------------------------

lfc_mat <- angiogenesis_lfc_clean %>%
  dplyr::select(gene, cluster, log2FoldChange) %>%
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = log2FoldChange,
    values_fn = mean
  ) %>%
  dplyr::filter(!is.na(gene), gene != "") %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

sig_mat <- angiogenesis_lfc_clean %>%
  dplyr::select(gene, cluster, sig_label) %>%
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = sig_label,
    values_fn = ~ .x[1]
  ) %>%
  dplyr::filter(!is.na(gene), gene != "") %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# asegurar mismo orden
common_genes <- intersect(rownames(lfc_mat), rownames(sig_mat))
common_clusters <- intersect(colnames(lfc_mat), colnames(sig_mat))

lfc_mat <- lfc_mat[common_genes, common_clusters, drop = FALSE]
sig_mat <- sig_mat[common_genes, common_clusters, drop = FALSE]

# reemplazar NAs en LFC por 0 o dejar NA
lfc_mat[is.na(lfc_mat)] <- 0
sig_mat[is.na(sig_mat)] <- ""




library(pheatmap)

# -----------------------------
# 5. Heatmap
# -----------------------------

heat_dir <- file.path(pb_dir, "plots", "Angiogenesis")
dir.create(heat_dir, recursive = TRUE, showWarnings = FALSE)

# Misma paleta que citoquinas
cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

# Breaks simétricos centrados en 0
max_abs <- max(abs(lfc_mat), na.rm = TRUE)

# Opcional: capar para que outliers no apaguen el resto
max_abs <- min(max_abs, 3)

breaks <- seq(-max_abs, max_abs, length.out = length(cols) + 1)

pheatmap(
  mat = pmax(pmin(lfc_mat, max_abs), -max_abs),
  color = cols,
  breaks = breaks,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  display_numbers = sig_mat,
  number_color = "black",
  fontsize_number = 10,
  fontsize_row = 10,
  fontsize_col = 10,
  angle_col = 45,
  border_color = NA,
  main = "Angiogenesis / vascular remodeling genes\nlog2FC KO vs WT adjusted by DsRed",
  legend_breaks = c(-max_abs, 0, max_abs),
  legend_labels = c(
    paste0("-", max_abs, "\nHigher in WT"),
    "0",
    paste0("+", max_abs, "\nHigher in KO")
  ),
  filename = file.path(
    heat_dir,
    "Heatmap_Angiogenesis_log2FC_KO_vs_WT_adjDsRed_clustered_samePalette2.pdf"
  ),
  width = 10,
  height = 20
)


write.csv(
  angiogenesis_lfc_clean,
  file = file.path(heat_dir, "Angiogenesis_genes_log2FC_by_cluster.csv"),
  row.names = FALSE
)



# Clusters donde vimos el programa WT-high más claro
wt_focus_clusters <- c(
  "IFN Mac",
  "Early IFN TAMs",
  "Ly6c2Hi Monocytes"
)

# Resumen por gen
program_gene_summary <- angiogenesis_lfc_clean %>%
  dplyr::mutate(
    is_sig_005 = !is.na(padj) & padj < 0.05,
    is_sig_010 = !is.na(padj) & padj < 0.10
  ) %>%
  dplyr::group_by(gene) %>%
  dplyr::summarise(
    mean_lfc_all = mean(log2FoldChange, na.rm = TRUE),
    median_lfc_all = median(log2FoldChange, na.rm = TRUE),
    max_lfc_all = max(log2FoldChange, na.rm = TRUE),
    min_lfc_all = min(log2FoldChange, na.rm = TRUE),
    
    mean_lfc_wt_focus = mean(
      log2FoldChange[cluster %in% wt_focus_clusters],
      na.rm = TRUE
    ),
    min_lfc_wt_focus = min(
      log2FoldChange[cluster %in% wt_focus_clusters],
      na.rm = TRUE
    ),
    
    n_clusters_tested = dplyr::n(),
    
    n_KO_up_lfc05 = sum(log2FoldChange > 0.5, na.rm = TRUE),
    n_WT_up_lfc05 = sum(log2FoldChange < -0.5, na.rm = TRUE),
    
    n_KO_up_sig005 = sum(log2FoldChange > 0 & is_sig_005, na.rm = TRUE),
    n_WT_up_sig005 = sum(log2FoldChange < 0 & is_sig_005, na.rm = TRUE),
    
    n_KO_up_sig010 = sum(log2FoldChange > 0 & is_sig_010, na.rm = TRUE),
    n_WT_up_sig010 = sum(log2FoldChange < 0 & is_sig_010, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_lfc_wt_focus)



  wt_program_candidates <- program_gene_summary %>%
  dplyr::filter(
    mean_lfc_wt_focus < -0.5 |
      n_WT_up_sig005 >= 2 |
      n_WT_up_lfc05 >= 3
  ) %>%
  dplyr::arrange(mean_lfc_wt_focus)

as.data.frame(wt_program_candidates)



ko_program_candidates <- program_gene_summary %>%
  dplyr::filter(
    mean_lfc_all > 0.25 |
      n_KO_up_sig005 >= 2 |
      n_KO_up_lfc05 >= 3
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_KO_up_sig005),
    dplyr::desc(n_KO_up_lfc05),
    dplyr::desc(mean_lfc_all)
  )

as.data.frame(ko_program_candidates)


functional_angio_score_genes <- c(...)
imperfect_ecm_angio_score_genes <- c(...)





# -----------------------------
# Genes finales para heatmap/scores
# -----------------------------

wt_functional_angio_genes <- c(
  "Adm",
  "Hbegf",
  "Pgf",
  "Vegfa",
  "Flt1",
  "Nrp2",
  "Igf1",
  "Areg",
  "Nos2",
  "Hif1a",
  "Jag1",
  "Jag2",
  "Mmp9",
  "Mmp12",
  "Mmp19",
  "Ctsk",
  "Timp1",
  "Timp3",
  "Plau",
  "Plaur",
  "Vcam1",
  "Icam1",
  "Il6",
  "Il1a"
)

ko_imperfect_ecm_genes <- c(
  "Osm",
  "Tgfb3",
  "Tgfbr1",
  "Col4a1",
  "Col4a2",
  "Col5a1",
  "Col5a2",
  "Col6a1",
  "Col6a2",
  "Col6a3",
  "Col18a1",
  "Col3a1",
  "Tnc",
  "Lama4",
  "Lamb1",
  "Nid1",
  "Pdgfrb",
  "Flt4",
  "Kdr",
  "Aplnr",
  "Bgn",
  "Adamts2",
  "Cxcl12"
)

focused_angio_genes <- unique(c(
  wt_functional_angio_genes,
  ko_imperfect_ecm_genes
))



# -----------------------------
# Tabla enfocada
# -----------------------------

angio_focus_lfc <- angiogenesis_lfc_clean %>%
  dplyr::filter(gene %in% focused_angio_genes) %>%
  dplyr::mutate(
    program = dplyr::case_when(
      gene %in% wt_functional_angio_genes ~ "WT_high_functional_angio_remodeling",
      gene %in% ko_imperfect_ecm_genes ~ "KO_high_ECM_imperfect_remodeling",
      TRUE ~ "Other"
    )
  )

write.csv(
  angio_focus_lfc,
  file = file.path(
    heat_dir,
    "Focused_Angiogenesis_ECM_program_genes_log2FC_table.csv"
  ),
  row.names = FALSE
)



# -----------------------------
# Matriz log2FC
# -----------------------------

lfc_mat_focus <- angio_focus_lfc %>%
  dplyr::select(gene, cluster, log2FoldChange) %>%
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = log2FoldChange,
    values_fn = mean
  ) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

sig_mat_focus <- angio_focus_lfc %>%
  dplyr::select(gene, cluster, sig_label) %>%
  tidyr::pivot_wider(
    names_from = cluster,
    values_from = sig_label,
    values_fn = ~ .x[1]
  ) %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# Mismo orden
common_genes <- intersect(rownames(lfc_mat_focus), rownames(sig_mat_focus))
common_clusters <- intersect(colnames(lfc_mat_focus), colnames(sig_mat_focus))

lfc_mat_focus <- lfc_mat_focus[common_genes, common_clusters, drop = FALSE]
sig_mat_focus <- sig_mat_focus[common_genes, common_clusters, drop = FALSE]

lfc_mat_focus[is.na(lfc_mat_focus)] <- 0
sig_mat_focus[is.na(sig_mat_focus)] <- ""


gene_order_focus <- c(
  wt_functional_angio_genes,
  ko_imperfect_ecm_genes
)

gene_order_focus <- gene_order_focus[gene_order_focus %in% rownames(lfc_mat_focus)]

lfc_mat_focus <- lfc_mat_focus[gene_order_focus, , drop = FALSE]
sig_mat_focus <- sig_mat_focus[gene_order_focus, , drop = FALSE]

row_annot <- data.frame(
  Program = ifelse(
    rownames(lfc_mat_focus) %in% wt_functional_angio_genes,
    "WT-high functional angiogenesis/remodeling",
    "KO-high ECM/imperfect remodeling"
  )
)

rownames(row_annot) <- rownames(lfc_mat_focus)



library(pheatmap)

cols <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)

max_abs <- max(abs(lfc_mat_focus), na.rm = TRUE)
max_abs <- min(max_abs, 3)

breaks <- seq(-max_abs, max_abs, length.out = length(cols) + 1)

mat_plot <- pmax(pmin(lfc_mat_focus, max_abs), -max_abs)

pheatmap(
  mat = mat_plot,
  color = cols,
  breaks = breaks,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  clustering_method = "complete",
  annotation_row = row_annot,
  gaps_row = length(wt_functional_angio_genes[
    wt_functional_angio_genes %in% rownames(mat_plot)
  ]),
  display_numbers = sig_mat_focus,
  number_color = "black",
  fontsize_number = 10,
  fontsize_row = 9,
  fontsize_col = 10,
  angle_col = 45,
  border_color = NA,
  main = "Focused angiogenesis / ECM remodeling programs\nlog2FC KO vs WT adjusted by DsRed",
  legend_breaks = c(-max_abs, 0, max_abs),
  legend_labels = c(
    paste0("-", max_abs, "\nHigher in WT"),
    "0",
    paste0("+", max_abs, "\nHigher in KO")
  ),
  filename = file.path(
    heat_dir,
    "Heatmap_Focused_Angiogenesis_ECM_programs_log2FC_KO_vs_WT_adjDsRed.pdf"
  ),
  width = 11,
  height = 10
)






# -----------------------------
# Heatmap transpuesto:
# genes en eje X, clusters en eje Y
# -----------------------------
annotation_colors <- list(
  Program = c(
    "WT-high functional angiogenesis/remodeling" = "#00CD66",  # verde/teal suave
    "KO-high ECM/imperfect remodeling" = "#CD4F39"             # coral suave
  )
)



# -----------------------------
# Heatmap transpuesto:
# genes en eje X, clusters en eje Y
# -----------------------------

mat_plot_t <- t(mat_plot)
sig_mat_focus_t <- t(sig_mat_focus)

# Anotación de columnas = genes
col_annot <- data.frame(
  Program = ifelse(
    colnames(mat_plot_t) %in% wt_functional_angio_genes,
    "WT-high functional angiogenesis/remodeling",
    "KO-high ECM/imperfect remodeling"
  )
)

rownames(col_annot) <- colnames(mat_plot_t)

# Colores bonitos para los dos programas
annotation_colors <- list(
  Program = c(
    "WT-high functional angiogenesis/remodeling" = "#66C2A5",
    "KO-high ECM/imperfect remodeling" = "#FC8D62"
  )
)

# Gap entre programas
gap_col <- length(wt_functional_angio_genes[
  wt_functional_angio_genes %in% colnames(mat_plot_t)
])

# Etiquetas de la leyenda
max_abs_label <- round(max_abs, 1)

pheatmap(
  mat = mat_plot_t,
  color = cols,
  breaks = breaks,
  scale = "none",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  clustering_method = "complete",
  annotation_col = col_annot,
  annotation_colors = annotation_colors,
  gaps_col = gap_col,
  display_numbers = sig_mat_focus_t,
  number_color = "black",
  fontsize_number = 9,
  fontsize_row = 10,
  fontsize_col = 9,
  angle_col = 45,
  border_color = NA,
  main = "Focused angiogenesis / ECM remodeling programs\nKO vs WT adjusted by DsRed",
  legend_breaks = c(-max_abs, 0, max_abs),
  legend_labels = c(
    paste0("log2FC\n-", max_abs_label, "\nHigher in WT"),
    "0",
    paste0("+", max_abs_label, "\nHigher in KO")
  ),
  filename = file.path(
    heat_dir,
    "Heatmap_Focused_Angiogenesis_ECM_programs_log2FC_KO_vs_WT_adjDsRed_transposed.pdf"
  ),
  width = 22,
  height = 5
)

###



# -----------------------------
# ORA de las dos firmas angiogénicas/ECM
# -----------------------------

library(dplyr)
library(tidyr)
library(msigdbr)
library(clusterProfiler)
library(ggplot2)
library(forcats)

ora_dir <- file.path(pb_dir, "ORA_Focused_Angiogenesis_ECM_signatures")
dir.create(ora_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ora_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(ora_dir, "plots"), recursive = TRUE, showWarnings = FALSE)



signatures_angio <- list(
  WT_high_functional_angio_remodeling = wt_functional_angio_genes,
  KO_high_ECM_imperfect_remodeling = ko_imperfect_ecm_genes
)

signatures_angio <- lapply(signatures_angio, unique)

signatures_angio



universe_genes <- unique(unlist(lapply(de_list, function(x) {
  
  res <- as.data.frame(x$res)
  
  if (!"gene" %in% colnames(res)) {
    res$gene <- rownames(res)
  }
  
  res %>%
    dplyr::filter(!is.na(stat)) %>%
    dplyr::pull(gene)
})))

length(universe_genes)
head(universe_genes)



msig_hallmark <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "MH"
)

msig_go_bp <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "M5",
  subcollection = "GO:BP"
)

msig_go_mf <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "M5",
  subcollection = "GO:MF"
)

msig_go_cc <- msigdbr(
  db_species = "MM",
  species = "Mus musculus",
  collection = "M5",
  subcollection = "GO:CC"
)

msig_go <- bind_rows(
  msig_go_bp,
  msig_go_mf,
  msig_go_cc
)



make_term2gene <- function(msig_df) {
  
  gene_col <- dplyr::case_when(
    "gene_symbol" %in% colnames(msig_df) ~ "gene_symbol",
    "db_gene_symbol" %in% colnames(msig_df) ~ "db_gene_symbol",
    TRUE ~ NA_character_
  )
  
  if (is.na(gene_col)) {
    stop("No encuentro columna de símbolo génico en msigdbr.")
  }
  
  msig_df %>%
    dplyr::select(term = gs_name, gene = all_of(gene_col)) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::distinct(term, gene)
}

term2gene_hallmark <- make_term2gene(msig_hallmark)
term2gene_go <- make_term2gene(msig_go)




run_ora_signature <- function(
    genes,
    signature_name,
    term2gene,
    collection_name,
    universe_genes,
    minGSSize = 3,
    maxGSSize = 500
) {
  
  genes <- unique(genes)
  term_genes <- unique(term2gene$gene)
  
  universe_use <- intersect(universe_genes, term_genes)
  genes_use <- intersect(genes, universe_use)
  
  message(
    "Running ORA: ", signature_name,
    " | ", collection_name,
    " | genes used: ", length(genes_use)
  )
  
  if (length(genes_use) < 3) {
    return(data.frame())
  }
  
  enr <- clusterProfiler::enricher(
    gene = genes_use,
    universe = universe_use,
    TERM2GENE = term2gene,
    pAdjustMethod = "BH",
    pvalueCutoff = 1,
    qvalueCutoff = 1,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize
  )
  
  enr_df <- as.data.frame(enr)
  
  if (nrow(enr_df) == 0) {
    return(data.frame())
  }
  
  enr_df %>%
    dplyr::mutate(
      signature = signature_name,
      collection = collection_name,
      input_genes_total = length(genes),
      input_genes_used = length(genes_use),
      input_genes_used_list = paste(genes_use, collapse = ",")
    ) %>%
    dplyr::select(
      signature,
      collection,
      ID,
      Description,
      GeneRatio,
      BgRatio,
      pvalue,
      p.adjust,
      qvalue,
      Count,
      geneID,
      input_genes_total,
      input_genes_used,
      input_genes_used_list
    ) %>%
    dplyr::arrange(p.adjust, pvalue)
}



ora_results <- list()

for (sig_name in names(signatures_angio)) {
  
  sig_genes <- signatures_angio[[sig_name]]
  
  ora_results[[paste(sig_name, "Hallmark", sep = "__")]] <- run_ora_signature(
    genes = sig_genes,
    signature_name = sig_name,
    term2gene = term2gene_hallmark,
    collection_name = "Hallmark",
    universe_genes = universe_genes
  )
  
  ora_results[[paste(sig_name, "GO", sep = "__")]] <- run_ora_signature(
    genes = sig_genes,
    signature_name = sig_name,
    term2gene = term2gene_go,
    collection_name = "GO",
    universe_genes = universe_genes
  )
}

ora_all <- bind_rows(ora_results)

as.data.frame(ora_all)



write.csv(
  ora_all,
  file = file.path(
    ora_dir,
    "tables",
    "ORA_focused_angio_ECM_signatures_Hallmark_GO_all.csv"
  ),
  row.names = FALSE
)



ora_sig <- ora_all %>%
  dplyr::filter(!is.na(p.adjust), p.adjust < 0.05) %>%
  dplyr::arrange(signature, collection, p.adjust)

write.csv(
  ora_sig,
  file = file.path(
    ora_dir,
    "tables",
    "ORA_focused_angio_ECM_signatures_Hallmark_GO_significant_padj005.csv"
  ),
  row.names = FALSE
)

as.data.frame(ora_sig)



ora_top20 <- ora_all %>%
  dplyr::filter(!is.na(p.adjust)) %>%
  dplyr::group_by(signature, collection) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::ungroup()

write.csv(
  ora_top20,
  file = file.path(
    ora_dir,
    "tables",
    "ORA_focused_angio_ECM_signatures_top20_by_signature_collection.csv"
  ),
  row.names = FALSE
)

as.data.frame(ora_top20)




plot_ora_signature <- function(ora_df, sig_name, top_n = 15) {
  
  plot_df <- ora_df %>%
    dplyr::filter(
      signature == sig_name,
      !is.na(p.adjust)
    ) %>%
    dplyr::group_by(collection) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      minus_log10_padj = -log10(p.adjust),
      Description = stringr::str_replace_all(Description, "_", " "),
      Description = stringr::str_to_sentence(Description),
      Description = forcats::fct_reorder(Description, minus_log10_padj)
    )
  
  ggplot(
    plot_df,
    aes(
      x = minus_log10_padj,
      y = Description,
      size = Count,
      fill = collection
    )
  ) +
    geom_point(shape = 21, color = "black", alpha = 0.85) +
    facet_wrap(~ collection, scales = "free_y") +
    theme_bw() +
    theme(
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.text.y = element_text(size = 8)
    ) +
    labs(
      title = paste0("ORA: ", sig_name),
      subtitle = "Hallmark and GO | universe = genes tested in pseudobulk DESeq2",
      x = "-log10 adjusted p-value",
      y = "",
      size = "Genes",
      fill = "Collection"
    )
}




for (sig_name in names(signatures_angio)) {
  
  p <- plot_ora_signature(
    ora_df = ora_all,
    sig_name = sig_name,
    top_n = 15
  )
  
  print(p)
  
  ggsave(
    filename = file.path(
      ora_dir,
      "plots",
      paste0("ORA_", sig_name, "_Hallmark_GO_top15.pdf")
    ),
    plot = p,
    width = 10,
    height = 8
  )
  
  ggsave(
    filename = file.path(
      ora_dir,
      "plots",
      paste0("ORA_", sig_name, "_Hallmark_GO_top15.png")
    ),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300
  )
}



library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(stringr)

# -----------------------------
# 1. Preparar tabla ORA comparativa SOLO GO
# -----------------------------

parse_ratio <- function(x) {
  sapply(strsplit(as.character(x), "/"), function(y) {
    as.numeric(y[1]) / as.numeric(y[2])
  })
}

ora_compare_go <- ora_all %>%
  dplyr::filter(
    collection == "GO",
    !is.na(p.adjust)
  ) %>%
  dplyr::mutate(
    gene_ratio_num = parse_ratio(GeneRatio),
    bg_ratio_num = parse_ratio(BgRatio),
    fold_enrichment = gene_ratio_num / bg_ratio_num,
    minus_log10_padj = -log10(pmax(p.adjust, 1e-300)),
    pathway_clean = Description %>%
      stringr::str_replace_all("_", " ") %>%
      stringr::str_replace("^GOBP ", "") %>%
      stringr::str_replace("^GOMF ", "") %>%
      stringr::str_replace("^GOCC ", "") %>%
      stringr::str_to_sentence()
  )

  # -----------------------------
# 2. Seleccionar top GO pathways de cada programa
# -----------------------------

top_terms_go <- ora_compare_go %>%
  dplyr::group_by(signature) %>%
  dplyr::arrange(p.adjust, .by_group = TRUE) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::ungroup() %>%
  dplyr::pull(ID) %>%
  unique()

ora_compare_go_plot <- ora_compare_go %>%
  dplyr::filter(ID %in% top_terms_go) %>%
  dplyr::mutate(
    signature = factor(
      signature,
      levels = c(
        "WT_high_functional_angio_remodeling",
        "KO_high_ECM_imperfect_remodeling"
      ),
      labels = c(
        "WT-high functional\nangiogenesis/remodeling",
        "KO-high ECM / imperfect\nvascular remodeling"
      )
    ),
    pathway_clean = forcats::fct_reorder(
      pathway_clean,
      minus_log10_padj,
      .fun = max
    )
  )

  # -----------------------------
# 3. Dotplot comparativo SOLO GO
# -----------------------------

p_ora_compare_go <- ggplot(
  ora_compare_go_plot,
  aes(
    x = signature,
    y = pathway_clean,
    size = fold_enrichment,
    color = minus_log10_padj
  )
) +
  geom_point(alpha = 0.95) +
  scale_color_gradient(
    low = "#FEE0D2",
    high = "#99000D"
  ) +
  scale_size_continuous(
    range = c(2.5, 9)
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
    axis.text.y = element_text(size = 9),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold"),
    legend.title = element_text(face = "bold")
  ) +
  labs(
    title = "GO enrichment of focused angiogenesis / ECM programs",
    subtitle = "ORA using GO terms | universe = genes tested in pseudobulk Analysis
                        All macrophages clusters combined", 
    x = "",
    y = "",
    size = "Fold enrichment",
    color = "-log10 adj. p"
  )

p_ora_compare_go


p_ora_compare_go <- p_ora_compare_go +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6)
  )

p_ora_compare_go


ggsave(
  filename = file.path(
    ora_dir,
    "plots",
    "GO_ORA_comparison_WTfunctional_vs_KO_ECM_programs_noInternalGrid.pdf"
  ),
  plot = p_ora_compare_go,
  width = 10,
  height = 9
)

ggsave(
  filename = file.path(
    ora_dir,
    "plots",
    "GO_ORA_comparison_WTfunctional_vs_KO_ECM_programs_noInternalGrid.png"
  ),
  plot = p_ora_compare_go,
  width = 10,
  height = 9,
  dpi = 300
)