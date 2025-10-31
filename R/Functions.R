
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




read.data <- function(path, sample, group, DsRed, library, batch){
data<- readRDS(paste0(datadir,path))
data$tag <- paste0(sample)
data$group <- group
data$DsRed <- DsRed
data$orig.ident <- library
data$batch <- batch
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
  data$metric <- sign(data$logFC) * -log10(data$P.Value)
  
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
                           cutoff = 0.7,
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
