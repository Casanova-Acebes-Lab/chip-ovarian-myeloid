library(Seurat)
library(dplyr)
library(cowplot)
library(tidyr)
library(fgsea)
library(msigdbr)
library(ggplot2)






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



## Reading data 


de_list <- readRDS(
  file.path(
    res_dir,
    "de_list_KO_vs_WT_adjusted_DsRed_all_clusters.rds"
  )
)

de_summary <- read.csv(
  file.path(
    res_dir,
    "DESeq2_KO_vs_WT_adjusted_DsRed_all_clusters_summary.csv"
  )
)


### GSEA analysis



# ============================================================
# GSEA / fgsea
# KO vs WT adjusted by DsRed
# Macrophages + Monocytes, cluster by cluster
# ============================================================

# -----------------------------
# 0. Rutas
# -----------------------------


pb_dir <- file.path(outdir, "Pseudobulk_Macrophages_Monocytes_Round3")

res_dir <- file.path(
  pb_dir,
  "DESeq2_KO_vs_WT_adjusted_DsRed_all_clusters"
)

gsea_dir <- file.path(
  pb_dir,
  "GSEA_KO_vs_WT_adjusted_DsRed_all_clusters"
)

dir.create(gsea_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(gsea_dir, "tables"), recursive = TRUE, showWarnings = FALSE)



# -----------------------------
# 2. Gene sets MSigDB mouse-native
# Hallmark + GO only
# -----------------------------

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

# Checks
msigdbr_collections(db_species = "MM")

nrow(msig_hallmark)
nrow(msig_go)

head(msig_hallmark)
head(msig_go)




# -----------------------------
# 3. Convertir MSigDB a pathways para fgsea
# -----------------------------

make_pathways <- function(msig_df) {
  
  gene_col <- dplyr::case_when(
    "gene_symbol" %in% colnames(msig_df) ~ "gene_symbol",
    "db_gene_symbol" %in% colnames(msig_df) ~ "db_gene_symbol",
    TRUE ~ NA_character_
  )
  
  if (is.na(gene_col)) {
    stop("No encuentro columna de símbolo génico en msigdbr.")
  }
  
  msig_df %>%
    dplyr::select(gs_name, gene = all_of(gene_col)) %>%
    dplyr::filter(!is.na(gene), gene != "") %>%
    dplyr::distinct(gs_name, gene) %>%
    split(x = .$gene, f = .$gs_name)
}

pathways_hallmark <- make_pathways(msig_hallmark)
pathways_go <- make_pathways(msig_go)

length(pathways_hallmark)
length(pathways_go)



# -----------------------------
# 4. Ranking desde DESeq2
# -----------------------------

make_rank_from_deseq <- function(res_df) {
  
  res_df <- as.data.frame(res_df)
  
  if (!"gene" %in% colnames(res_df)) {
    res_df$gene <- rownames(res_df)
  }
  
  ranks_df <- res_df %>%
    dplyr::filter(
      !is.na(gene),
      !is.na(stat)
    ) %>%
    dplyr::select(gene, stat) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      stat = stat[which.max(abs(stat))],
      .groups = "drop"
    )
  
  ranks <- ranks_df$stat
  names(ranks) <- ranks_df$gene
  
  sort(ranks, decreasing = TRUE)
}

# Check rápido
ranks_test <- make_rank_from_deseq(de_list[[1]]$res)

head(ranks_test)
tail(ranks_test)
length(ranks_test)



# -----------------------------
# 5. Función fgsea
# -----------------------------

run_fgsea_one <- function(res_df, pathways, cluster_name, collection_name) {
  
  ranks <- make_rank_from_deseq(res_df)
  
  fg <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranks,
    minSize = 10,
    maxSize = 500
  )
  
  fg <- as.data.frame(fg) %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(
      cluster = cluster_name,
      collection = collection_name,
      contrast = "KO_vs_WT_adjDsRed",
      direction = dplyr::case_when(
        NES > 0 ~ "Enriched_in_KO",
        NES < 0 ~ "Enriched_in_WT",
        TRUE ~ "Neutral"
      ),
      leadingEdge = sapply(leadingEdge, paste, collapse = ",")
    ) %>%
    dplyr::select(
      contrast,
      cluster,
      collection,
      pathway,
      direction,
      NES,
      pval,
      padj,
      size,
      leadingEdge
    )
  
  fg
}



# -----------------------------
# 6. Correr GSEA por cluster
# Hallmark + GO only
# -----------------------------

gsea_results <- list()

for (cl in names(de_list)) {
  
  message("Running GSEA: ", cl)
  
  res_df <- de_list[[cl]]$res
  
  gsea_results[[paste(cl, "Hallmark", sep = "__")]] <- run_fgsea_one(
    res_df = res_df,
    pathways = pathways_hallmark,
    cluster_name = cl,
    collection_name = "Hallmark"
  )
  
  gsea_results[[paste(cl, "GO", sep = "__")]] <- run_fgsea_one(
    res_df = res_df,
    pathways = pathways_go,
    cluster_name = cl,
    collection_name = "GO"
  )
}

gsea_all <- bind_rows(gsea_results)



# -----------------------------
# 7. Guardar resultados completos
# -----------------------------

write.csv(
  gsea_all,
  file = file.path(
    gsea_dir,
    "tables",
    "GSEA_all_clusters_Hallmark_GO.csv"
  ),
  row.names = FALSE
)

saveRDS(
  gsea_all,
  file = file.path(
    gsea_dir,
    "GSEA_all_clusters_Hallmark_GO.rds"
  )
)



# -----------------------------
# 8. Significativos
# -----------------------------

gsea_sig_005 <- gsea_all %>%
  dplyr::filter(!is.na(padj), padj < 0.05) %>%
  dplyr::arrange(collection, cluster, padj)

gsea_sig_010 <- gsea_all %>%
  dplyr::filter(!is.na(padj), padj < 0.10) %>%
  dplyr::arrange(collection, cluster, padj)

write.csv(
  gsea_sig_005,
  file = file.path(
    gsea_dir,
    "tables",
    "GSEA_significant_padj005_Hallmark_GO.csv"
  ),
  row.names = FALSE
)

write.csv(
  gsea_sig_010,
  file = file.path(
    gsea_dir,
    "tables",
    "GSEA_significant_padj010_Hallmark_GO.csv"
  ),
  row.names = FALSE
)


# -----------------------------
# 9. Top 20 por cluster y colección
# -----------------------------

gsea_top20_by_cluster_collection <- gsea_all %>%
  dplyr::filter(!is.na(padj)) %>%
  dplyr::group_by(cluster, collection) %>%
  dplyr::arrange(padj, .by_group = TRUE) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::ungroup()

write.csv(
  gsea_top20_by_cluster_collection,
  file = file.path(
    gsea_dir,
    "tables",
    "GSEA_top20_by_cluster_collection_Hallmark_GO.csv"
  ),
  row.names = FALSE
)

as.data.frame(gsea_top20_by_cluster_collection)



# -----------------------------
# 10. Resumen
# -----------------------------

gsea_summary <- gsea_all %>%
  dplyr::group_by(cluster, collection) %>%
  dplyr::summarise(
    n_pathways_tested = dplyr::n(),
    n_sig_005 = sum(padj < 0.05, na.rm = TRUE),
    n_sig_010 = sum(padj < 0.10, na.rm = TRUE),
    n_KO_005 = sum(padj < 0.05 & NES > 0, na.rm = TRUE),
    n_WT_005 = sum(padj < 0.05 & NES < 0, na.rm = TRUE),
    top_pathway = pathway[which.min(padj)],
    top_NES = NES[which.min(padj)],
    top_padj = min(padj, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(cluster, collection)

write.csv(
  gsea_summary,
  file = file.path(
    gsea_dir,
    "tables",
    "GSEA_summary_by_cluster_collection_Hallmark_GO.csv"
  ),
  row.names = FALSE
)

as.data.frame(gsea_summary)





as.data.frame(gsea_summary)
as.data.frame(gsea_sig_005)
as.data.frame(gsea_top20_by_cluster_collection)



# Plotting things


plot_dir <- file.path(gsea_dir, "plots_top8_Hallmark_GO")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

library(dplyr)
library(ggplot2)

clean_gsea_name <- function(x) {
  x %>%
    gsub("^HALLMARK_", "", .) %>%
    gsub("^GOBP_", "", .) %>%
    gsub("^GOMF_", "", .) %>%
    gsub("^GOCC_", "", .) %>%
    gsub("_", " ", .)
}

sig_stars <- function(padj) {
  dplyr::case_when(
    is.na(padj) ~ "",
    padj < 0.001 ~ "***",
    padj < 0.01  ~ "**",
    padj < 0.05  ~ "*",
    padj < 0.10  ~ ".",
    TRUE ~ ""
  )
}

make_top_table <- function(gsea_df, cluster_name, collection_name, n_top = 8) {
  
  gsea_df %>%
    dplyr::filter(
      cluster == cluster_name,
      collection == collection_name,
      !is.na(padj)
    ) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = n_top) %>%
    dplyr::mutate(
      pathway_clean = clean_gsea_name(pathway),
      direction = dplyr::case_when(
        NES > 0 ~ "Enriched in KO",
        NES < 0 ~ "Enriched in WT",
        TRUE ~ "Neutral"
      ),
      sig_label = sig_stars(padj)
    )
}



plot_gsea_top_cluster <- function(gsea_df, cluster_name, n_top = 20) {
  
  top_h <- make_top_table(
    gsea_df = gsea_df,
    cluster_name = cluster_name,
    collection_name = "Hallmark",
    n_top = n_top
  )
  
  top_go <- make_top_table(
    gsea_df = gsea_df,
    cluster_name = cluster_name,
    collection_name = "GO",
    n_top = n_top
  )
  
  top_df <- dplyr::bind_rows(top_h, top_go)
  
  if (nrow(top_df) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = paste("No GSEA results for", cluster_name)) +
        theme_void()
    )
  }
  
  # ordenar dentro de cada colección
  top_df <- top_df %>%
    dplyr::group_by(collection) %>%
    dplyr::arrange(NES, .by_group = TRUE) %>%
    dplyr::mutate(
      pathway_plot = paste(collection, pathway_clean, sep = "___")
    ) %>%
    dplyr::ungroup()
  
  levs <- top_df %>%
    dplyr::arrange(collection, NES) %>%
    dplyr::pull(pathway_plot) %>%
    unique()
  
  top_df$pathway_plot <- factor(top_df$pathway_plot, levels = levs)
  
  max_abs_nes <- max(abs(top_df$NES), na.rm = TRUE)
  x_pad <- max(0.4, max_abs_nes * 0.18)
  
  top_df <- top_df %>%
    dplyr::mutate(
      star_x = ifelse(NES >= 0, NES + x_pad * 0.25, NES - x_pad * 0.25),
      star_hjust = ifelse(NES >= 0, 0, 1)
    )
  
  p <- ggplot(top_df, aes(x = NES, y = pathway_plot, fill = direction)) +
    geom_col(width = 0.75) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey40") +
    geom_text(
      aes(x = star_x, label = sig_label),
      hjust = top_df$star_hjust,
      size = 5,
      fontface = "bold"
    ) +
    facet_grid(collection ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) gsub("^.*___", "", x)) +
    scale_fill_manual(
      values = c(
        "Enriched in KO" = "#D55E00",
        "Enriched in WT" = "#0072B2",
        "Neutral" = "grey70"
      )
    ) +
    coord_cartesian(
      xlim = c(-max_abs_nes - x_pad, max_abs_nes + x_pad)
    ) +
    labs(
      title = paste0(cluster_name, " | Top ", n_top, " Hallmark + Top ", n_top, " GO"),
      subtitle = "NES > 0: enriched in KO | NES < 0: enriched in WT",
      x = "Normalized Enrichment Score (NES)",
      y = NULL,
      fill = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text.y = element_text(face = "bold"),
      axis.text.y = element_text(size = 9),
      legend.position = "bottom"
    )
  
  return(p)
}



clusters_gsea <- unique(gsea_all$cluster)

gsea_plot_list <- list()

for (cl in clusters_gsea) {
  
  message("Plotting GSEA: ", cl)
  
  p <- plot_gsea_top_cluster(
    gsea_df = gsea_all,
    cluster_name = cl,
    n_top = 60
  )
  
  gsea_plot_list[[cl]] <- p
  
  safe_cl <- gsub("[^A-Za-z0-9_\\-]", "_", cl)
  
  ggsave(
    filename = file.path(plot_dir, paste0(safe_cl, "_GSEA_top8_Hallmark_GO.png")),
    plot = p,
    width = 14,
    height = 22,
    dpi = 300
  )
  
  ggsave(
    filename = file.path(plot_dir, paste0(safe_cl, "_GSEA_top8_Hallmark_GO.pdf")),
    plot = p,
    width = 14,
    height = 22
  )
}