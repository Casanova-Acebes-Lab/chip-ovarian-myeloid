
PCA <- function(data){
## Selecting PCAs

# First metric 
#The point where the principal components only contribute 5% of standard deviation and the principal components cumulatively contribute 90% of the standard deviation.

# Determine percent of variation associated with each PC
pct <- data[["pca"]]@stdev / sum(data[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1

# Second metric
# The point where the percent change in variation between the consequtive PCs is less than 0.1%.

# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# last point where change of % of variation is more than 0.1%.

co2
list <- list(co1,co2)
return(list)

}



umap <- function(data,pc){

data <- FindNeighbors(data, dims = 1:paste0(pc))
data <- FindClusters(data, resolution = 0.5)

data <- RunUMAP(data, dims = 1:paste0(pc))

return(data)

}



# Función para procesar una librería
process_library <- function(base_dir, output_dir) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  sample_dirs <- list.dirs(base_dir, recursive = FALSE)
  results_list <- list()
  
  for (sample_dir in sample_dirs) {
    sample_name <- basename(sample_dir)
    cat("\nProcesando muestra:", sample_name, "\n")
    
    # Leer matrices raw y filtered
    raw_mat <- Read10X(file.path(sample_dir, "raw_feature_bc_matrix"))
    filt_mat <- Read10X(file.path(sample_dir, "filtered_feature_bc_matrix"))
    
    # Mantener solo genes comunes
    common_genes <- intersect(rownames(raw_mat), rownames(filt_mat))
    raw_mat <- raw_mat[common_genes, ]
    filt_mat <- filt_mat[common_genes, ]
    
    # Crear SoupChannel
    sc <- SoupChannel(tod = raw_mat, toc = filt_mat, channelName = sample_name)
    
    # Clustering temporal con Seurat
    seurat_tmp <- CreateSeuratObject(counts = filt_mat)
    seurat_tmp <- SCTransform(seurat_tmp, verbose = FALSE)
    seurat_tmp <- RunPCA(seurat_tmp, verbose = FALSE)
    seurat_tmp <- FindNeighbors(seurat_tmp, dims = 1:20)
    seurat_tmp <- FindClusters(seurat_tmp, resolution = 0.5)
    sc <- setClusters(sc, clusters = Idents(seurat_tmp))
    
    # Estimar contaminación
    sc <- autoEstCont(sc)
    
    # Mostrar rho estimado
    rho_value <- sc$metaData$rho
    cat("Tasa de contaminación estimada (rho) para", sample_name, ":", rho_value, "\n")
    
    # Corregir counts
    corrected_matrix <- adjustCounts(sc)
    
    # Crear Seurat object final
    seurat_obj <- CreateSeuratObject(counts = corrected_matrix)
    
    # Guardar Seurat object en disco
    saveRDS(seurat_obj, file = file.path(output_dir, paste0("Seurat_", sample_name, ".rds")))
    
    # Guardar rho para CSV
    results_list[[sample_name]] <- rho_value
  }
  
  # Crear CSV con las tasas de contaminación
  rho_df <- data.frame(
    sample = names(results_list),
    rho = unlist(results_list)
  )
  write.csv(rho_df, file = file.path(output_dir, paste0("rho_", basename(base_dir), ".csv")), row.names = FALSE)
  
  cat("\nProcesamiento de librería", basename(base_dir), "completado.\n")
}




process_seurat_samples <- function(seurat_list, seurat_names) {
  
  if(length(seurat_list) != length(seurat_names)) {
    stop("La longitud de seurat_list y seurat_names debe ser la misma.")
  }
  
  for (i in seq_along(seurat_list)) {
    cat("Procesando muestra:", seurat_names[i], "\n")
    
    # Normalizar
    seurat_list[[i]] <- NormalizeData(seurat_list[[i]], verbose = FALSE)
    
    # Encontrar genes variables
    seurat_list[[i]] <- FindVariableFeatures(seurat_list[[i]], selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    
    # Escalar usando solo genes variables
    seurat_list[[i]] <- ScaleData(seurat_list[[i]], features = VariableFeatures(seurat_list[[i]]), verbose = FALSE)
    
    # PCA
    seurat_list[[i]] <- RunPCA(seurat_list[[i]], features = VariableFeatures(seurat_list[[i]]), verbose = FALSE)
    
    # Asignar al entorno global
    assign(seurat_names[i], seurat_list[[i]], envir = .GlobalEnv)
  }
  
  cat("Procesamiento completado para todas las muestras.\n")
}






plot <- function(data,sample){

data$multi <- data@meta.data[7]

percent <- table(data$multi)[1] * 100 / (table(data$multi)[1] + table(data$multi)[2])
title <- paste0("Doublets ", table(data$multi)[1],
" Singlets ", table(data$multi)[2], 
" Percent of Doublets ", round(percent, digits = 4))

p <- DimPlot(data, reduction = "umap", group.by = "multi") +
ggtitle(sample,title)

return(p)

}




read.data <- function(path, sample, group, DsRed, library, batch, chimera){
data<- readRDS(paste0(datadir,path))
data$tag <- paste0(sample)
data$group <- group
data$DsRed <- DsRed
data$orig.ident <- library
data$batch <- batch
data$chimera <- chimera
data <- RenameCells(data, add.cell.id = sample)  

 if (ncol(data@meta.data) >= 6) {
    colnames(data@meta.data)[6] <- "pANN.doublet"
  }
  if (ncol(data@meta.data) >= 7) {
    colnames(data@meta.data)[7] <- "DF.class"
  }
  return(data)
}





GSEA <- function(data) {
 
  # Añadir el símbolo como columna
  data$mgi_symbol <- rownames(data)
  
  # Mapear símbolos a Entrez IDs usando org.Mm.eg.db (más confiable que biomaRt)
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = rownames(data),
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  # Añadir la columna con Entrez al dataframe
  data$entrezgene_id <- entrez_ids[rownames(data)]
  
  # Filtrar filas sin EntrezID
  data <- data[!is.na(data$entrezgene_id), ]
  
  # Quitar genes duplicados (si algún símbolo mapeó a la misma ID)
  data <- data[!duplicated(data$entrezgene_id), ]
  
  # Reasignar rownames por EntrezID
  rownames(data) <- data$entrezgene_id
  
  # Reemplazar p_val == 0 por un valor muy pequeño para evitar -Inf
  epsilon <- 1e-300
  data$adj.P.Val[data$adj.P.Val == 0] <- epsilon
  
  # Calcular la métrica de ranking
  data$metric <- data$logFC * -log10(data$adj.P.Val + 1e-8)

  
  # Ordenar por la métrica de ranking
  data <- data[order(data$metric, decreasing = TRUE), ]
  
  # Crear el vector nombrado para GSEA
  gene_metric <- data$metric
  names(gene_metric) <- rownames(data)
  
  # Ejecutar GSEA (ontología biológica por defecto)
  gseGO_result <- gseGO(geneList = gene_metric,
                        ont = "BP",
                        OrgDb = org.Mm.eg.db,
                        minGSSize = 150,
                        maxGSSize = 500,
                        eps = 1e-20,
                        nPermSimple = 10000,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05,
                        verbose = FALSE)
  
  # Simplificar GO (eliminar términos redundantes)
  gseGO_result <- simplify(gseGO_result,
                           cutoff = 0.9,
                           by = "p.adjust",
                           select_fun = min)
  
  # Hacer resultados legibles
  gseGO_result <- setReadable(gseGO_result,
                              OrgDb = org.Mm.eg.db,
                              keyType = "ENTREZID")
  
  return(gseGO_result)
}







GSEA2 <- function(data) {
 

  # Añadir el símbolo como columna
  data$mgi_symbol <- rownames(data)
  
  # Mapear símbolos a Entrez IDs usando org.Mm.eg.db (más confiable que biomaRt)
  entrez_ids <- mapIds(org.Mm.eg.db,
                       keys = rownames(data),
                       column = "ENTREZID",
                       keytype = "SYMBOL",
                       multiVals = "first")
  
  # Añadir la columna con Entrez al dataframe
  data$entrezgene_id <- entrez_ids[rownames(data)]
  
  # Filtrar filas sin EntrezID
  data <- data[!is.na(data$entrezgene_id), ]
  
  # Quitar genes duplicados (si algún símbolo mapeó a la misma ID)
  data <- data[!duplicated(data$entrezgene_id), ]
  
  # Reasignar rownames por EntrezID
  rownames(data) <- data$entrezgene_id
  
  # Reemplazar p_val == 0 por un valor muy pequeño para evitar -Inf
  epsilon <- 1e-300
  data$adj.P.Val[data$p_val_adj == 0] <- epsilon
  
  # Calcular la métrica de ranking
  data$metric <- sign(data$avg_log2FC) 
  
  # Ordenar por la métrica de ranking
  data <- data[order(data$metric, decreasing = TRUE), ]
  
  # Crear el vector nombrado para GSEA
  gene_metric <- data$metric
  names(gene_metric) <- rownames(data)
  
  # Ejecutar GSEA (ontología biológica por defecto)
  gseGO_result <- gseGO(geneList = gene_metric,
                        ont = "BP",
                        OrgDb = org.Mm.eg.db,
                        minGSSize = 50,
                        maxGSSize = 500,
                        eps = 1e-20,
                        nPermSimple = 10000,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.005,
                        verbose = FALSE)
  
  # Simplificar GO (eliminar términos redundantes)
  gseGO_result <- simplify(gseGO_result,
                           cutoff = 0.9,
                           by = "p.adjust",
                           select_fun = min)
  
  # Hacer resultados legibles
  gseGO_result <- setReadable(gseGO_result,
                              OrgDb = org.Mm.eg.db,
                              keyType = "ENTREZID")
  
  return(gseGO_result)
}





plot_gsea <- function(gsea_df, title = "GSEA Plot") {
  gsea_df <- gsea_df %>%
    mutate(
      NES_sign = ifelse(NES > 0, "Positive", "Negative"),
      NES_sign = factor(NES_sign, levels = c("Positive","Negative"))
    ) %>%
    dplyr::arrange(NES_sign, dplyr::desc(NES)) %>%   # usar dplyr:: explícito
    mutate(Description = factor(Description, levels = rev(Description))) 

  ggplot(gsea_df, aes(x = Description, y = NES, fill = NES_sign)) +
    geom_col() +
    scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue")) +
    coord_flip() +
    labs(
      x = "Pathway",
      y = "Normalized Enrichment Score",
      title = title
    ) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      legend.position = "none"
    )
}





# --------------------------------------------
# 1️⃣ Pseudobulk con raw counts
# --------------------------------------------
pseudobulk_cluster <- function(seurat_obj, cluster_name, cluster_col="Clustering.Round2", sample_col="tag"){
  # Matriz de counts raw
  
  counts <- GetAssayData(seurat_obj, assay="RNA", slot="counts")
  meta   <- seurat_obj@meta.data
  
  # Seleccionar células del cluster
  cells <- rownames(meta)[meta[[cluster_col]] == cluster_name]
  counts <- counts[, cells, drop=FALSE]
  
  # Vector de samples
  samples <- as.character(meta[cells, sample_col])
  
  # Matriz pseudobulk
  pb_mat <- t(aggregate.Matrix(t(counts), groupings=samples, fun="sum"))
  
  # Número de células por sample
  n_cells <- table(samples)
  n_cells <- n_cells[colnames(pb_mat)]
  
  # Total counts por sample
  n_counts <- colSums(pb_mat)
  
  # Media de counts por célula
  avg_counts_per_cell <- n_counts / n_cells


  # Número de genes detectados por célula (>=1 UMI)
n_genes_cell <- colSums(counts > 0)

# Promedio por muestra
n_genes_per_sample <- tapply(n_genes_cell, samples, mean)
  
  list(pb_mat=pb_mat, n_cells=n_cells, n_counts=n_counts, avg_counts_per_cell=avg_counts_per_cell,
       n_genes_per_sample=n_genes_per_sample)
}


# --------------------------------------------
# 2️⃣ Metadata de samples
# --------------------------------------------


make_sample_metadata <- function(sample_names, n_cells, n_counts, avg_counts_per_cell,
                                 n_genes_per_sample){
  depthPerGene <- n_counts / n_genes_per_sample  # <-- aquí calculamos correctamente
  
  df <- data.frame(
    sample = sample_names,
    genotype = ifelse(grepl("KO", sample_names), "KO", "WT"),
    dsred    = ifelse(grepl("DsRedP", sample_names), "DsRedP", "DsRedN"),
    nCells   = n_cells,
    nCounts  = n_counts,
    avgCountsPerCell = avg_counts_per_cell,
    nGenesPerCell = n_genes_per_sample,
    depthPerGene = depthPerGene
  )
  
  rownames(df) <- df$sample
  df$genotype <- factor(df$genotype, levels=c("WT","KO"))
  df
}


# --------------------------------------------
# 3️⃣ Limma-voom DE por cluster
# --------------------------------------------
run_limma_cluster <- function(pb_mat, sample_meta){
  
  # Crear objeto DGE
  dge <- DGEList(counts = pb_mat)

  # Filtrado por CPM: al menos 1 CPM en >=3 muestras
  cpm_counts <- cpm(dge)
  keep <- rowSums(cpm_counts > 1) >= 3
  dge <- dge[keep, , keep.lib.sizes = FALSE]

  # Normalización de factores de escala (TMM)
  dge <- calcNormFactors(dge)

  # Diseño: covariable técnica + genotipo
  design <- model.matrix(~ genotype + depthPerGene, data = sample_meta)

  v <- voom(dge, design, plot=FALSE)
  expr_normalized <- v$E  # matriz lista para GSVA

  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  # Contraste KO vs WT
  cont.matrix <- makeContrasts(KOvsWT = genotypeKO, levels = design)
  fit2 <- contrasts.fit(fit, cont.matrix)
  fit2 <- eBayes(fit2)

  # Tabla DE
  tt <- topTable(fit2, coef="KOvsWT", number=Inf, sort.by="none")
  tt <- tt[order(tt$adj.P.Val), ]
  tt$gene <- rownames(tt)

  # Devolver lista
  list(
    tt = tt,
    expr_normalized = expr_normalized
  )
}





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



# --------
# --------------------------------------------
# 5️⃣ Pipeline completo por cluster
# --------------------------------------------
run_pseudobulk_analysis <- function(data, cluster_name, outdir, rsdir){


  # 1️⃣ Pseudobulk
  res <- pseudobulk_cluster(data, cluster_name)
  pb_mat <- res$pb_mat
  
  # 2️⃣ Metadata

  sample_meta <- make_sample_metadata(
    colnames(pb_mat),
    n_cells = res$n_cells,
    n_counts = res$n_counts,
    avg_counts_per_cell = res$avg_counts_per_cell,
    n_genes_per_sample = res$n_genes_per_sample 
  )

colnames(sample_meta)[grepl("nCells.Freq", colnames(sample_meta))] <- "nCells"
colnames(sample_meta)[grepl("avgCountsPerCell.Freq", colnames(sample_meta))] <- "avgCountsPerCell"

# Mantener solo las columnas necesarias y limpiar nombres
sample_meta <- sample_meta[, c("sample", "genotype", "depthPerGene", "nCells")]
rownames(sample_meta) <- sample_meta$sample

# 3️⃣ Limma-voom DE
res_limma <- run_limma_cluster(pb_mat, sample_meta)
expr_normalized <- res_limma$expr_normalized  # matriz lista para GSVA
res_limma <- res_limma$tt                     # tabla DE para los siguientes pasos

  
  # Clasificación DEGs
  res_limma$diffexpressed <- "NO"
  res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
  res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"
  
  # Top genes
  top_genes <- res_limma %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 50) %>%
    as.data.frame()
  
  # Guardar tabla DE
  write.table(
    res_limma,
    file = paste0(rsdir, "/table.macros.", cluster_name, ".tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # Volcano plot
  p <- Volcano2(res_limma, legend = paste0("KO vs WT samples. ", cluster_name))
  pdf(paste0(outdir, "/Pseudobulk/Volcano.KO_vs_WT_", cluster_name, ".pdf"), width=16, height=12)
  print(p)
  dev.off()
  

  # Return
  list(
    cluster = cluster_name,
    res_limma = res_limma,
    top_genes = top_genes,
    n_cells = res$n_cells,
    n_counts = res$n_counts,
    avg_counts_per_cell = res$avg_counts_per_cell,
    plot = p,
    matrix = pb_mat,
    expr_normalized = expr_normalized,  # matriz lista para GSVA
    sample_meta = sample_meta
)
}


pseudobulk_limma <- function(
  seurat_obj,
  cluster_name,
  cluster_col = "Clustering.Round2",
  sample_col = "tag",
  genotype_col = "group",
  chimera_col = "chimera",
  min_cpm = 1,
  min_samples = 3
){
  library(Matrix.utils)
  library(edgeR)
  library(limma)
  
  # 1️⃣ Extraer counts y metadata
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  meta <- seurat_obj@meta.data
  
  # 2️⃣ Seleccionar células del cluster
  cells <- rownames(meta)[meta[[cluster_col]] == cluster_name]
  if(length(cells) == 0) stop("No se encontraron células para el cluster")
  
  counts <- counts[, cells, drop = FALSE]
  meta <- meta[cells, , drop = FALSE]
  meta <- meta[colnames(counts), , drop = FALSE]
  
  # 3️⃣ Crear pseudobulk
  sample_vector <- as.character(meta[[sample_col]])
  pb_mat <- t(aggregate.Matrix(t(counts), groupings = sample_vector, fun = "sum"))
  
  # 4️⃣ Metadata por sample
  # Contar células por pseudobulk
  n_cells <- sapply(colnames(pb_mat), function(s) sum(meta[[sample_col]] == s))
  # Calcular counts per cell
  total_counts <- colSums(pb_mat)
  counts_per_cell <- total_counts / n_cells
  
  sample_meta <- data.frame(
    sample = colnames(pb_mat),
    genotype = factor(ifelse(grepl("KO", colnames(pb_mat)), "KO", "WT"), levels = c("WT","KO")),
    chimera = factor(meta[[chimera_col]][match(colnames(pb_mat), sample_vector)]),
    n_cells = n_cells,
    total_counts = total_counts,
    counts_per_cell = counts_per_cell
  )
  rownames(sample_meta) <- sample_meta$sample
  
  # 5️⃣ Filtrado genes
  dge <- DGEList(counts = pb_mat)
  keep <- rowSums(cpm(dge) > min_cpm) >= min_samples
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  # 6️⃣ voom
  v <- voom(dge, plot = FALSE)
  
  # 7️⃣ Diseño limma usando counts_per_cell como covariable continua
  design <- model.matrix(~ genotype + counts_per_cell, data = sample_meta)
  
  # 8️⃣ Ajuste limma
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # 9️⃣ Contraste KO vs WT
  contrast_matrix <- makeContrasts(KOvsWT = genotypeKO, levels = design)
  fit2 <- contrasts.fit(fit, contrast_matrix)
  fit2 <- eBayes(fit2)
  tt <- topTable(fit2, coef = "KOvsWT", number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)

  # Clasificación DEGs
  tt$diffexpressed <- "NO"
  tt$diffexpressed[tt$adj.P.Val < 0.05 & tt$logFC > 0.5] <- "Up"
  tt$diffexpressed[tt$adj.P.Val < 0.05 & tt$logFC < -0.5] <- "Down"
  
  
  # 🔟 Devolver resultados
  return(list(
    pb_mat = pb_mat,
    sample_meta = sample_meta,
    voom_expr = v$E,
    limma_fit = fit2,
    tt = tt
  ))
}


pseudobulk_limma_avg_per_cell_with_umis <- function(
  seurat_obj,
  cluster_name,
  cluster_col = "Clustering.Round2",
  sample_col = "tag",
  genotype_col = "group",
  min_cpm = 1,
  min_samples = 3
){
  library(Matrix.utils)
  library(edgeR)
  library(limma)
  
  # 1️⃣ Extraer counts y metadata
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  meta <- seurat_obj@meta.data
  
  # 2️⃣ Seleccionar células del cluster
  cells <- rownames(meta)[meta[[cluster_col]] == cluster_name]
  if(length(cells) == 0) stop("No se encontraron células para el cluster")
  
  counts <- counts[, cells, drop = FALSE]
  meta <- meta[cells, , drop = FALSE]
  meta <- meta[colnames(counts), , drop = FALSE]
  
  # 3️⃣ Crear pseudobulk promedio por célula
  sample_vector <- as.character(meta[[sample_col]])
  n_cells <- sapply(unique(sample_vector), function(s) sum(sample_vector == s))
  
  pb_mat_sum <- t(aggregate.Matrix(t(counts), groupings = sample_vector, fun = "sum"))
  pb_mat_avg <- sweep(pb_mat_sum, 2, n_cells[colnames(pb_mat_sum)], FUN = "/")
  
  # 4️⃣ Calcular UMIs per cell
  total_counts <- colSums(pb_mat_sum)
  umis_per_cell <- total_counts / n_cells[colnames(pb_mat_sum)]
  
  # 5️⃣ Metadata por sample
  sample_meta <- data.frame(
    sample = colnames(pb_mat_avg),
    genotype = factor(meta[[genotype_col]][match(colnames(pb_mat_avg), sample_vector)],
                      levels = c("WT","KO")),
    n_cells = n_cells[colnames(pb_mat_avg)],
    umis_per_cell = umis_per_cell
  )
  rownames(sample_meta) <- sample_meta$sample
  
  # 6️⃣ Filtrado genes
  dge <- DGEList(counts = pb_mat_avg)
  keep <- rowSums(cpm(dge) > min_cpm) >= min_samples
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  # 7️⃣ voom
  v <- voom(dge, plot = FALSE)
  
  # 8️⃣ Diseño limma con genotipo
  design <- model.matrix(~ genotype, data = sample_meta)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # 9️⃣ Contraste KO vs WT
  coef_name <- "genotypeKO"
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  
  # Clasificación DEGs
  tt$diffexpressed <- "NO"
  tt$diffexpressed[tt$adj.P.Val < 0.05 & tt$logFC > 0.5]  <- "Up"
  tt$diffexpressed[tt$adj.P.Val < 0.05 & tt$logFC < -0.5] <- "Down"
  
  # 🔟 Matriz corregida para PCA/MDS usando UMIs per cell como covariable
  corrected_expr <- removeBatchEffect(v$E, covariates = sample_meta$umis_per_cell)
  
  return(list(
    pb_mat_avg    = pb_mat_avg,
    sample_meta   = sample_meta,
    voom_expr     = v$E,
    corrected_expr = corrected_expr,
    limma_fit     = fit,
    tt            = tt
  ))
}








pseudobulk_limma_SVA <- function(
  seurat_obj,
  cluster_name,
  cluster_col = "Clustering.Round2",
  sample_col = "tag",
  genotype_col = "group",
  min_cpm = 1,
  min_samples = 3,
  n_sv = 2
){
  library(Matrix.utils)
  library(edgeR)
  library(limma)
  library(sva)
  
  ### 1) Extraer counts y metadata
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  meta <- seurat_obj@meta.data
  
  ### 2) Seleccionar células del cluster
  cells <- rownames(meta)[meta[[cluster_col]] == cluster_name]
  if(length(cells) == 0) stop("No se encontraron células para el cluster especificado.")
  
  counts <- counts[, cells, drop = FALSE]
  meta <- meta[cells, , drop = FALSE]
  meta <- meta[colnames(counts), , drop = FALSE]
  
  ### 3) Crear pseudobulk
  sample_vector <- as.character(meta[[sample_col]])
  pb_mat <- t(aggregate.Matrix(t(counts), groupings = sample_vector, fun = "sum"))
  
  ### 4) Metadata por muestra
  sample_meta <- data.frame(
    sample    = colnames(pb_mat),
    genotype  = factor(meta[[genotype_col]][match(colnames(pb_mat), sample_vector)])
  )
  sample_meta$genotype <- factor(sample_meta$genotype, levels = c("WT", "KO"))
levels(sample_meta$genotype)
# [1] "WT" "KO"

  rownames(sample_meta) <- sample_meta$sample
  
  ### 5) Filtrado
  dge <- DGEList(counts = pb_mat)
  keep <- rowSums(cpm(dge) > min_cpm) >= min_samples
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  ### 6) Modelos para SVA
  mod  <- model.matrix(~ genotype, data = sample_meta)
  mod0 <- model.matrix(~ 1, data = sample_meta)
  
  ### 7) Estimar factores latentes
  svobj <- svaseq(as.matrix(dge$counts), mod = mod, mod0 = mod0, n.sv = n_sv)
  
  for(i in 1:n_sv){
    sample_meta[[paste0("SV", i)]] <- svobj$sv[, i]
  }
  
  ### 8) Diseño final
  design <- model.matrix(
    as.formula(
      paste("~ genotype +", paste(paste0("SV",1:n_sv), collapse=" + "))
    ),
    data = sample_meta
  )
  
  ### 9) Voom + limma
  v <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  ### 10) Contraste
  coef_name <- grep("genotype", colnames(design), value = TRUE)[1]
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  
  ### Clasificación DEGs
  tt$diffexpressed <- "NO"
  tt$diffexpressed[tt$adj.P.Val < 0.05 & tt$logFC >  0.5] <- "Up"
  tt$diffexpressed[tt$adj.P.Val < 0.05 & tt$logFC < -0.5] <- "Down"
  
  ### 11) MATRIZ CORREGIDA PARA MDS/PCA
  corrected_expr <- removeBatchEffect(
    v$E,
    covariates = svobj$sv
  )
  
  ### 12) Return
  return(list(
    pb_mat        = pb_mat,
    sample_meta   = sample_meta,
    dge           = dge,
    svobj         = svobj,
    design        = design,
    voom_expr     = v$E,
    corrected_expr = corrected_expr,   # <<--- NUEVA SALIDA
    fit           = fit,
    tt            = tt
  ))
}




pseudobulk_limma_RUVr <- function(
  seurat_obj,
  cluster_name,
  cluster_col = "Clustering.Round2",
  sample_col = "tag",
  genotype_col = "group",
  min_cpm = 1,
  min_samples = 3,
  k = 2  # número de factores latentes
){
  library(Matrix.utils)
  library(edgeR)
  library(limma)
  library(RUVSeq)
  
  ### 1) Extraer counts y metadata
  counts <- GetAssayData(seurat_obj, assay = "RNA", layer = "counts")
  meta <- seurat_obj@meta.data
  
  ### 2) Seleccionar células del cluster
  cells <- rownames(meta)[meta[[cluster_col]] == cluster_name]
  if(length(cells) == 0) stop("No se encontraron células para el cluster especificado.")
  
  counts <- counts[, cells, drop = FALSE]
  meta <- meta[cells, , drop = FALSE]
  meta <- meta[colnames(counts), , drop = FALSE]
  
  ### 3) Crear pseudobulk
  sample_vector <- as.character(meta[[sample_col]])
  pb_mat <- t(aggregate.Matrix(t(counts), groupings = sample_vector, fun = "sum"))
  
  ### 4) Metadata de muestra
  sample_meta <- data.frame(
    sample    = colnames(pb_mat),
    genotype  = factor(meta[[genotype_col]][match(colnames(pb_mat), sample_vector)],
                       levels = c("WT", "KO"))
  )
  rownames(sample_meta) <- sample_meta$sample
  
  ### 5) Filtrado EdgeR
  dge <- DGEList(counts = pb_mat)
  keep <- rowSums(cpm(dge) > min_cpm) >= min_samples
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  dge <- calcNormFactors(dge)
  
  ### 6) Modelo de interés
  mod <- model.matrix(~ genotype, data = sample_meta)
  
  ### 7) Identificar genes “invariantes” para RUVr
  # Usamos todos los genes con baja variación como negativos
  counts_for_ruv <- as.matrix(dge$counts)
  counts_for_ruv <- counts_for_ruv[keep, ]
  
  # Ajuste inicial voom
  v <- voom(dge, design = mod, plot = FALSE)
  
  # Residuals
  fit_init <- lmFit(v, mod)
  resids <- v$E - v$E %*% coef(fit_init) %*% t(mod)  # residuals matrix
  
  # RUVr: genera k factores latentes
  set <- newSeqExpressionSet(counts_for_ruv,
                             phenoData = data.frame(sample = colnames(counts_for_ruv), row.names = colnames(counts_for_ruv)))
  
  ruv <- RUVr(set, cIdx = rownames(counts_for_ruv), k = k, residuals = resids)
  
  # Añadir factores latentes al metadata
  for(i in 1:k){
    sample_meta[[paste0("W_",i)]] <- pData(ruv)[, paste0("W_",i)]
  }
  
  ### 8) Nuevo diseño con factores latentes
  design <- model.matrix(
    as.formula(
      paste("~ genotype +", paste(paste0("W_",1:k), collapse=" + "))
    ),
    data = sample_meta
  )
  
  ### 9) voom + limma con RUVr
  v2 <- voom(dge, design, plot = FALSE)
  fit <- lmFit(v2, design)
  fit <- eBayes(fit)
  
  ### 10) Contraste KO vs WT
  coef_name <- "genotypeKO"
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  tt$gene <- rownames(tt)
  
  ### Clasificación DEGs
  tt$diffexpressed <- "NO"
  tt$diffexpressed[tt$adj.P.Val < 0.05 & tt$logFC > 0.5]  <- "Up"
  tt$diffexpressed[tt$adj.P.Val < 0.05 & tt$logFC < -0.5] <- "Down"
  
  ### 11) Matriz corregida para MDS/PCA
  corrected_expr <- removeBatchEffect(v2$E,
                                      covariates = as.matrix(pData(ruv)[, paste0("W_",1:k)]))
  
  ### 12) Return
  return(list(
    pb_mat        = pb_mat,
    sample_meta   = sample_meta,
    dge           = dge,
    ruv           = ruv,
    design        = design,
    voom_expr     = v2$E,
    corrected_expr = corrected_expr,
    fit           = fit,
    tt            = tt
  ))
}
