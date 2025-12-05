library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(limma)
library(edgeR)
library(GSVA)
library(pheatmap)
library(msigdbr)
library(RColorBrewer)
library(matrixStats)
library(grid)






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



colnames(read.table(paste0(rsdir,"matrix.macros.Ly6cHi Monocytes.tsv"),
                           header = TRUE,
                           row.names = 1,
                           sep = "\t",
                           check.names = FALSE))

                           colnames(read.table(paste0(rsdir,"matrix.macros.Ly6cLo Monocytes.tsv"),
                           header = TRUE,
                           row.names = 1,
                           sep = "\t",
                           check.names = FALSE))
colnames(read.table(paste0(rsdir,"matrix.macros.Mrc1 C1qc Cbr2 Gas6 Mac.tsv"),
                           header = TRUE,
                           row.names = 1,
                           sep = "\t",
                           check.names = FALSE))



# lista de clusters con nombres "limpios"
cluster_files <- c(
  "Ly6cHi_Monocytes"               = "matrix.macros.Ly6cHi Monocytes.tsv",
  "Ly6cLo_Monocytes"               = "matrix.macros.Ly6cLo Monocytes.tsv",
  "Early_IFN_MHCII_TAMs"           = "matrix.macros.Early IFN MHCII TAMs.tsv",
  "Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac" = "matrix.macros.Arg1 Spp1 Mmp12 Mmp19 Il1a Mac.tsv",
  "Trem1_Ptgs2_Plaur_Celc4e_Mac"   = "matrix.macros.Trem1 Ptgs2 Plaur Celc4e Mac.tsv",
  "MHCII_Ccl12_Mac"                = "matrix.macros.MHCII Ccl12 Mac.tsv",
  "MHCII_Siglec_Mac"               = "matrix.macros.MHCII Siglec Mac.tsv",
  "IFN_Mac"                        = "matrix.macros.IFN Mac.tsv",
  "Mmp9_Ctsk_Mac"                  = "matrix.macros.Mmp9 Ctsk Mac.tsv",
  "Mrc1_C1qc_Cbr2_Gas6_Mac"        = "matrix.macros.Mrc1 C1qc Cbr2 Gas6 Mac.tsv",
  "Npr2_Actn1_Mac"                 = "matrix.macros.Npr2 Actn1 Mac.tsv",
  "Fn1_Vegfa_Mac"                  = "matrix.macros.Fn1 Vegfa Mac.tsv",
  "Neutrophils"                    = "matrix.macros.Neutrophils.tsv"
)

# leer todas las matrices
mat_list <- lapply(cluster_files, function(f) {
  m <- read.table(paste0(rsdir, f), header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  as.matrix(m)
})



expr_cluster_mat <- do.call(cbind, lapply(mat_list, function(m) {
  rowMeans(m)   # media log2CPM por cluster
}))

colnames(expr_cluster_mat) <- names(cluster_files)


dim(expr_cluster_mat)   # genes × 13 clusters





############################################################
# GSVA Top 50 Pathways por cluster (WT vs KO) - cluster a cluster
############################################################



# -------------------------------
# Parámetros
# -------------------------------
min_genes <- 30
max_genes <- 500
top_n <- 100
exclude_keywords <- c(
  "neuro", "axon", "brain", "neuronal", "synapse",
  "muscle", "myocyte", "cardiac", "heart",
  "eye", "retina", "lens", "spinal", "fear", "light"
)
rdbu_pal <- colorRampPalette(c("blue","white","red"))(100)

# -------------------------------
# Cluster a procesar (modificar según necesidad)
# -------------------------------
cluster_name <- "MHCII_Siglec_Mac"
cluster_file <- "matrix.macros.MHCII Siglec Mac.tsv"

# -------------------------------
# Leer matriz
# -------------------------------
m <- read.table(file.path(rsdir, cluster_file),
                header=TRUE, row.names=1, sep="\t", check.names=FALSE)
m[is.na(m)] <- 0

# -------------------------------
# Gene sets C5:BP
# -------------------------------
msig_c5 <- msigdbr(species="Mus musculus", category="C5", subcategory="BP")
gene_sets <- split(msig_c5$gene_symbol, msig_c5$gs_name)

# -------------------------------
# Normalizar solo este cluster con voom
# -------------------------------
y <- DGEList(counts = m)
y <- calcNormFactors(y)
voom_mat <- voom(y)$E

# -------------------------------
# Filtrar gene sets según genes presentes
# -------------------------------
gs_filtered <- lapply(gene_sets, function(x) intersect(x, rownames(voom_mat)))
gs_filtered <- gs_filtered[sapply(gs_filtered, length) >= min_genes &
                           sapply(gs_filtered, length) <= max_genes]

# -------------------------------
# GSVA
# -------------------------------
params <- gsvaParam(exprData = voom_mat, geneSets = gs_filtered, kcdf = "Gaussian")
gs <- gsva(params, verbose=TRUE)

# -------------------------------
# WT vs KO
# -------------------------------
WT_cols <- grep("WT", colnames(gs))
KO_cols <- grep("KO", colnames(gs))
if(length(WT_cols)==0 | length(KO_cols)==0){
  stop("No se encontraron columnas WT o KO en este cluster")
}

gs_cluster <- cbind(
  WT = rowMeans(gs[, WT_cols, drop=FALSE]),
  KO = rowMeans(gs[, KO_cols, drop=FALSE])
)

# -------------------------------
# Filtrar pathways irrelevantes
# -------------------------------
gs_cluster <- gs_cluster[
  !grepl(paste(exclude_keywords, collapse="|"), rownames(gs_cluster), ignore.case=TRUE),
]

# -------------------------------
# Seleccionar top 50 por promedio KO
# -------------------------------
KO_sorted <- sort(gs_cluster[, "KO"], decreasing=TRUE)
KO_sorted <- KO_sorted[!duplicated(names(KO_sorted))]
top_pathways <- names(KO_sorted)[1:min(top_n, length(KO_sorted))]
gs_top <- gs_cluster[top_pathways, , drop=FALSE]

# -------------------------------
# Escalado row-wise (z-score)
# -------------------------------
gs_scaled <- t(scale(t(gs_top)))

# -------------------------------
# Guardar heatmap PDF
# -------------------------------
pdf(file.path(outdir, paste0("/GSVA/GSVA_heatmap_", cluster_name, "_WTvsKO2.pdf")), width=10, height=22)

annotation_col <- data.frame(Condition=c("WT","KO"))
rownames(annotation_col) <- colnames(gs_scaled)
ann_colors <- list(Condition=c(WT="#B0B0B0", KO="#404040"))

pheatmap(
  gs_scaled,
  cluster_rows=TRUE,
  cluster_cols=FALSE,
  annotation_col=annotation_col,
  annotation_colors=ann_colors,
  color=rdbu_pal,
  main=paste0("GSVA - ", cluster_name, " (WT vs KO)"),
  fontsize_row=6,
  border_color=NA
)

dev.off()








plot_gsva_cluster <- function(cl, cluster_file) {

  message("Procesando cluster: ", cl)

  # -------------------------------
  # Parámetros fijos
  # -------------------------------
  min_genes <- 30
  max_genes <- 500
  top_n <- 40
  exclude_keywords <- c(
    "neuro", "axon", "brain", "neuronal", "synapse",
    "muscle", "myocyte", "cardiac", "heart",
    "eye", "retina", "lens", "spinal", "fear", "light"
  )
  rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))

  # -------------------------------
  # Leer matriz
  # -------------------------------
  m <- read.table(
    paste0(rsdir, cluster_file),
    header = TRUE, row.names = 1, sep = "\t", check.names = FALSE
  )
  m[is.na(m)] <- 0

  # -------------------------------
  # Gene sets C5:BP
  # -------------------------------
  msig_c5 <- msigdbr(species="Mus musculus", category="C5", subcategory="BP")
  gene_sets <- split(msig_c5$gene_symbol, msig_c5$gs_name)

  # -------------------------------
  # Normalizar con voom
  # -------------------------------
  y <- DGEList(counts = m)
  y <- calcNormFactors(y)
  voom_mat <- voom(y)$E

  # -------------------------------
  # Filtrar gene sets
  # -------------------------------
  gs_filtered <- lapply(gene_sets, function(x) intersect(x, rownames(voom_mat)))
  gs_filtered <- gs_filtered[
    sapply(gs_filtered, length) >= min_genes &
    sapply(gs_filtered, length) <= max_genes
  ]

  # -------------------------------
  # GSVA
  # -------------------------------
  params <- gsvaParam(exprData = voom_mat, geneSets = gs_filtered, kcdf = "Gaussian")
  gs <- gsva(params)

  # Filtrar pathways irrelevantes
  gs <- gs[!grepl(paste(exclude_keywords, collapse="|"), rownames(gs), ignore.case=TRUE), ]

  # -------------------------------
  # Seleccionar top pathways según promedio KO
  # -------------------------------
  KO_cols <- grep("KO", colnames(gs))
  KO_means <- rowMeans(gs[, KO_cols, drop=FALSE])
  KO_sorted <- sort(KO_means, decreasing=TRUE)
  top_pathways <- names(KO_sorted)[1:min(top_n, length(KO_sorted))]
  gs_top <- gs[top_pathways, , drop=FALSE]

  # -------------------------------
  # Escalar filas (z-score)
  # -------------------------------
  gs_scaled <- t(scale(t(gs_top)))

  # -------------------------------
  # Orden correcto de columnas
  # -------------------------------
  col_order <- c(
    "WT1.DsRedN", "WT2_DsRedN",
    "WT1.DsRedP", "WT2_DsRedP",
    "DsRedN.KO2", "DsRedN.KO3",
    "DsRedP.KO2", "DsRedP.KO3"
  )
  col_order <- col_order[col_order %in% colnames(gs_scaled)]

  annotation_col <- data.frame(
    Condition = c(rep("WT", 4), rep("KO", 4)),
    SampleType = rep(c("DsRedN", "DsRedP"), each=2, times=2)
  )
  rownames(annotation_col) <- col_order

  ann_colors <- list(
    Condition = c(WT="orange3", KO="aquamarine4"),
    SampleType = c(DsRedN="#4876FF", DsRedP="#CD4F39")
  )

  # -------------------------------
  # Generar heatmap
  # -------------------------------
  p <- pheatmap(
    gs_scaled[, col_order, drop=FALSE],
    cluster_rows=TRUE,
    cluster_cols=FALSE,
    annotation_col=annotation_col,
    annotation_colors=ann_colors,
    color=rdbu_pal,
    main=paste0("GSVA - ", cl, " Cluster"),
    fontsize_row=6,
    border_color=NA
  )

  # -------------------------------
  # Guardar PDF con etiquetas
  # -------------------------------
  pdf(file.path(outdir, paste0("/GSVA/GSVA_heatmap_", cl, "_8cols_WTvsKO2.pdf")),
      width=14, height=40)

  grid.newpage()
  grid.draw(p$gtable)

  # Etiquetas extra
  grid.text("GSVA Score", x=0.85, y=0.93, gp=gpar(fontsize=10))
  grid.text("GO Terms", x=0.88, y=0.5, rot=90, gp=gpar(fontsize=10))

  dev.off()

  message("✔ Heatmap generado para ", cl)
}






plot_gsva_cluster <- function(cl, cluster_file) {

  message("Procesando cluster: ", cl)

  # -------------------------------
  # Parámetros fijos
  # -------------------------------
  min_genes <- 30
  max_genes <- 500
  top_n <- 40
  exclude_keywords <- c(
    "neuro", "axon", "brain", "neuronal", "synapse",
    "muscle", "myocyte", "cardiac", "heart",
    "eye", "retina", "lens",
    "spinal", "fear", "light"
  )
  rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))

  # -------------------------------
  # Preparar matriz
  # -------------------------------

  # Construir path completo
  full_path <- file.path(rsdir, cluster_file)
  
 voom_mat <- read.table(full_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)

# Convertir todo a numérico
voom_mat <- as.matrix(voom_mat)
mode(voom_mat) <- "numeric"   # fuerza que sea numérica

# Opcional: reemplazar NAs por 0
voom_mat[is.na(voom_mat)] <- 0

  

  # -------------------------------
  # Gene sets C5:BP
  # -------------------------------
  msig_c5 <- msigdbr(species="Mus musculus", category="C5", subcategory="BP")
  gene_sets <- split(msig_c5$gene_symbol, msig_c5$gs_name)

  # Filtrar gene sets
  gs_filtered <- lapply(gene_sets, function(x) intersect(x, rownames(voom_mat)))
  gs_filtered <- gs_filtered[
    sapply(gs_filtered, length) >= min_genes &
    sapply(gs_filtered, length) <= max_genes
  ]

  # -------------------------------
  # GSVA
  # -------------------------------
  params <- gsvaParam(exprData = voom_mat, geneSets = gs_filtered, kcdf = "Gaussian")
  gs <- gsva(params)

  # Filtrar pathways irrelevantes
  gs <- gs[!grepl(paste(exclude_keywords, collapse="|"), rownames(gs), ignore.case=TRUE), ]

  # -------------------------------
  # Seleccionar top pathways según promedio KO
  # -------------------------------
  KO_cols <- grep("KO", colnames(gs))
  KO_means <- rowMeans(gs[, KO_cols, drop=FALSE])
  KO_sorted <- sort(KO_means, decreasing=TRUE)
  top_pathways <- names(KO_sorted)[1:min(top_n, length(KO_sorted))]
  gs_top <- gs[top_pathways, , drop=FALSE]

  # Escalar filas (z-score)
  gs_scaled <- gs_top 

  # -------------------------------
  # Anotación automática de columnas
  # -------------------------------
  samples <- colnames(gs_scaled)
  Condition <- ifelse(grepl("KO", samples), "KO", "WT")
  SampleType <- ifelse(grepl("DsRedP", samples), "DsRedP", "DsRedN")
  annotation_col <- data.frame(Condition=Condition, SampleType=SampleType, row.names=samples)

  ann_colors <- list(
    Condition = c(WT="orange3", KO="aquamarine4"),
    SampleType = c(DsRedN="#4876FF", DsRedP="#CD4F39")
  )

  # -------------------------------
  # Generar heatmap
  # -------------------------------
  p <- pheatmap(
    gs_scaled,
    cluster_rows=TRUE,
    cluster_cols=FALSE,
    annotation_col=annotation_col,
    annotation_colors=ann_colors,
    color=rdbu_pal,
    main=paste0("GSVA - ", cl, " Cluster"),
    fontsize_row=6,
    border_color=NA
  )

  # -------------------------------
  # Guardar PDF
  # -------------------------------
  pdf(file.path(outdir, paste0("/GSVA/GSVA_heatmap_", cl, "_WTvsKO_normalized.pdf")),
      width=14, height=20)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grid::grid.text("GSVA Score", x=0.85, y=0.93, gp=grid::gpar(fontsize=10))
  grid::grid.text("GO Terms", x=0.88, y=0.5, rot=90, gp=grid::gpar(fontsize=10))
  dev.off()

  message("✔ Heatmap generado para ", cl)
}





# Lista de clusters y sus archivos
cluster_files <- c(
  "Ly6cHi_Monocytes"               = "matrix.macros.Ly6cHi Monocytes.tsv",
  "Ly6cLo_Monocytes"               = "matrix.macros.Ly6cLo Monocytes.tsv",
  "Early_IFN_MHCII_TAMs"           = "matrix.macros.Early IFN MHCII TAMs.tsv",
  "Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac" = "matrix.macros.Arg1 Spp1 Mmp12 Mmp19 Il1a Mac.tsv",
  "Trem1_Ptgs2_Plaur_Celc4e_Mac"   = "matrix.macros.Trem1 Ptgs2 Plaur Celc4e Mac.tsv",
  "MHCII_Ccl12_Mac"                = "matrix.macros.MHCII Ccl12 Mac.tsv",
  "MHCII_Siglec_Mac"               = "matrix.macros.MHCII Siglec Mac.tsv",
  "IFN_Mac"                        = "matrix.macros.IFN Mac.tsv",
  "Mmp9_Ctsk_Mac"                  = "matrix.macros.Mmp9 Ctsk Mac.tsv",
  "Mrc1_C1qc_Cbr2_Gas6_Mac"        = "matrix.macros.Mrc1 C1qc Cbr2 Gas6 Mac.tsv",
  "Npr2_Actn1_Mac"                 = "matrix.macros.Npr2 Actn1 Mac.tsv",
  "Fn1_Vegfa_Mac"                  = "matrix.macros.Fn1 Vegfa Mac.tsv",
  "Neutrophils"                    = "matrix.macros.Neutrophils.tsv"
)




lapply(names(cluster_files), function(cl) {
  plot_gsva_cluster(
    cl = cl,
    cluster_file = cluster_files[[cl]]  # aquí sí es solo el nombre del archivo
  )
})



plot_gsva_cluster <- function(cl, cluster_file, rsdir, outdir) {

  message("Procesando cluster: ", cl)

  # -------------------------------
  # Parámetros fijos
  # -------------------------------
  min_genes <- 30
  max_genes <- 500
  top_n <- 40
  exclude_keywords <- c(
    "neuro", "axon", "brain", "neuronal", "synapse",
    "muscle", "myocyte", "cardiac", "heart",
    "eye", "retina", "lens",
    "spinal", "fear", "light"
  )
  rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))

  # -------------------------------
  # Preparar matriz
  # -------------------------------
  full_path <- file.path(rsdir, cluster_file)
  voom_mat <- read.table(full_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  voom_mat <- as.matrix(voom_mat)
  mode(voom_mat) <- "numeric"
  voom_mat[is.na(voom_mat)] <- 0

  # -------------------------------
  # Gene sets C5:BP
  # -------------------------------
  msig_c5 <- msigdbr(species="Mus musculus", category="C5", subcategory="BP")
  gene_sets <- split(msig_c5$gene_symbol, msig_c5$gs_name)
  gs_filtered <- lapply(gene_sets, function(x) intersect(x, rownames(voom_mat)))
  gs_filtered <- gs_filtered[
    sapply(gs_filtered, length) >= min_genes &
    sapply(gs_filtered, length) <= max_genes
  ]

  # -------------------------------
  # GSVA
  # -------------------------------
  params <- gsvaParam(exprData = voom_mat, geneSets = gs_filtered, kcdf = "Gaussian")
  gs <- gsva(params)
  gs <- gs[!grepl(paste(exclude_keywords, collapse="|"), rownames(gs), ignore.case=TRUE), ]

  # -------------------------------
  # Seleccionar top pathways según promedio KO
  # -------------------------------
  KO_cols <- grep("KO", colnames(gs))
  KO_means <- rowMeans(gs[, KO_cols, drop=FALSE])
  KO_sorted <- sort(KO_means, decreasing=TRUE)
  top_pathways <- names(KO_sorted)[1:min(top_n, length(KO_sorted))]
  gs_top <- gs[top_pathways, , drop=FALSE]

  # -------------------------------
  # Agrupar columnas por tipo de muestra
  # -------------------------------
  groups <- list(
    KO_DsRedN = c("DsRedN-KO2", "DsRedN-KO3"),
    KO_DsRedP = c("DsRedP-KO2", "DsRedP-KO3"),
    WT_DsRedN = c("WT1-DsRedN", "WT2_DsRedN"),
    WT_DsRedP = c("WT1-DsRedP", "WT2_DsRed")
  )

  gs_grouped <- sapply(groups, function(cols) {
    cols_exist <- cols[cols %in% colnames(gs_top)]
    if (length(cols_exist) == 0) {
      warning("Ninguna columna encontrada para: ", paste(cols, collapse=", "))
      return(rep(NA, nrow(gs_top)))
    }
    rowMeans(gs_top[, cols_exist, drop=FALSE], na.rm=TRUE)
  })

  gs_grouped <- as.data.frame(gs_grouped)

  # Escalar filas (z-score)
  gs_scaled <- t(scale(t(gs_grouped)))

  # -------------------------------
  # Anotación automática de columnas
  # -------------------------------
  annotation_col <- data.frame(
    Condition = c("KO", "KO", "WT", "WT"),
    SampleType = c("DsRedN", "DsRedP", "DsRedN", "DsRedP"),
    row.names = colnames(gs_scaled)
  )

  ann_colors <- list(
    Condition = c(WT="orange3", KO="aquamarine4"),
    SampleType = c(DsRedN="#4876FF", DsRedP="#CD4F39")
  )

  # -------------------------------
  # Generar heatmap
  # -------------------------------
  p <- pheatmap(
    gs_scaled,
    cluster_rows=TRUE,
    cluster_cols=FALSE,
    annotation_col=annotation_col,
    annotation_colors=ann_colors,
    color=rdbu_pal,
    main=paste0("GSVA - ", cl, " Cluster"),
    fontsize_row=6,
    border_color=NA
  )

  # -------------------------------
  # Guardar PDF
  # -------------------------------
  pdf(file.path(outdir, paste0("GSVA/GSVA_heatmap_", cl, "_WTvsKO_normalized.pdf")),
      width=12, height=20)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grid::grid.text("GSVA Score", x=0.85, y=0.93, gp=grid::gpar(fontsize=10))
  grid::grid.text("GO Terms", x=0.88, y=0.5, rot=90, gp=grid::gpar(fontsize=10))
  dev.off()

  message("✔ Heatmap generado para ", cl)
}


# Lista de clusters y sus archivos
cluster_files <- c(
  "Ly6cHi_Monocytes"               = "matrix.macros.Ly6cHi Monocytes.tsv",
  "Ly6cLo_Monocytes"               = "matrix.macros.Ly6cLo Monocytes.tsv",
  "Early_IFN_MHCII_TAMs"           = "matrix.macros.Early IFN MHCII TAMs.tsv",
  "Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac" = "matrix.macros.Arg1 Spp1 Mmp12 Mmp19 Il1a Mac.tsv",
  "Trem1_Ptgs2_Plaur_Celc4e_Mac"   = "matrix.macros.Trem1 Ptgs2 Plaur Celc4e Mac.tsv",
  "MHCII_Ccl12_Mac"                = "matrix.macros.MHCII Ccl12 Mac.tsv",
  "MHCII_Siglec_Mac"               = "matrix.macros.MHCII Siglec Mac.tsv",
  "IFN_Mac"                        = "matrix.macros.IFN Mac.tsv",
  "Mmp9_Ctsk_Mac"                  = "matrix.macros.Mmp9 Ctsk Mac.tsv",
  "Mrc1_C1qc_Cbr2_Gas6_Mac"        = "matrix.macros.Mrc1 C1qc Cbr2 Gas6 Mac.tsv",
  "Npr2_Actn1_Mac"                 = "matrix.macros.Npr2 Actn1 Mac.tsv",
  "Fn1_Vegfa_Mac"                  = "matrix.macros.Fn1 Vegfa Mac.tsv",
  "Neutrophils"                    = "matrix.macros.Neutrophils.tsv"
)




# Ejecutar para todos los clusters
lapply(names(cluster_files), function(cl) {
  plot_gsva_cluster(
    cl = cl,
    cluster_file = cluster_files[[cl]],
    rsdir = rsdir,
    outdir = outdir
  )
})





library(GSVA)
library(msigdbr)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(dplyr)

plot_gsva_cluster_cytokine <- function(cl, cluster_file, rsdir, outdir) {
  message("\n==============================")
  message("Procesando cluster: ", cl)
  message("==============================\n")
  
  # -------------------------------
  # Parámetros fijos
  # -------------------------------
  min_genes <- 30
  max_genes <- 500
  top_n <- 40
  rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  
  # -------------------------------
  # Leer matriz desde archivo
  # -------------------------------
  full_path <- file.path(rsdir, cluster_file)
  voom_mat <- read.table(full_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  voom_mat <- as.matrix(voom_mat)
  mode(voom_mat) <- "numeric"
  voom_mat[is.na(voom_mat)] <- 0
  
  # -------------------------------
  # Gene sets GO BP relacionados con citoquinas
  # -------------------------------
  msig_go <- msigdbr(species="Mus musculus", category="C5", subcategory="BP")
  
  include_keywords <- c(
    "cytokine", "interleukin", "IL-", "IL[0-9]",
    "chemokine", "CXCL", "CCL", "TNF"
  )
  exclude_keywords <- c("cytokinesis", "cytokinetic")
  
  cytokine_go <- msig_go %>%
    filter(grepl(paste(include_keywords, collapse="|"), gs_name, ignore.case=TRUE)) %>%
    filter(!grepl(paste(exclude_keywords, collapse="|"), gs_name, ignore.case=TRUE))
  
  cytokine_list <- split(cytokine_go$gene_symbol, cytokine_go$gs_name) %>%
    lapply(unique) %>%
    Filter(function(x) length(x) >= min_genes & length(x) <= max_genes, .)
  
  message("Pathways seleccionados: ", length(cytokine_list))
  
  # -------------------------------
  # GSVA
  # -------------------------------
  params <- gsvaParam(exprData = voom_mat, geneSets = cytokine_list, kcdf = "Gaussian")
  gs <- gsva(params)
  
  # -------------------------------
  # Seleccionar top pathways según promedio KO
  # -------------------------------
  KO_cols <- grep("KO", colnames(gs))
  KO_means <- rowMeans(gs[, KO_cols, drop=FALSE])
  KO_sorted <- sort(KO_means, decreasing=TRUE)
  top_pathways <- names(KO_sorted)[1:min(top_n, length(KO_sorted))]
  gs_top <- gs[top_pathways, , drop=FALSE]
  
  # -------------------------------
  # Agrupar columnas por tipo de muestra
  # -------------------------------
  groups <- list(
    KO_DsRedN = c("DsRedN-KO2", "DsRedN-KO3"),
    KO_DsRedP = c("DsRedP-KO2", "DsRedP-KO3"),
    WT_DsRedN = c("WT1-DsRedN", "WT2_DsRedN"),
    WT_DsRedP = c("WT1-DsRedP", "WT2_DsRed")
  )
  
  gs_grouped <- sapply(groups, function(cols) {
    cols_exist <- cols[cols %in% colnames(gs_top)]
    if (length(cols_exist) == 0) {
      warning("Ninguna columna encontrada para: ", paste(cols, collapse=", "))
      return(rep(NA, nrow(gs_top)))
    }
    rowMeans(gs_top[, cols_exist, drop=FALSE], na.rm=TRUE)
  })
  
  gs_grouped <- as.data.frame(gs_grouped)
  
  # Escalar filas (z-score)
  gs_scaled <- t(scale(t(gs_grouped)))
  
  # -------------------------------
  # Anotación automática de columnas
  # -------------------------------
  annotation_col <- data.frame(
    Condition = c("KO", "KO", "WT", "WT"),
    SampleType = c("DsRedN", "DsRedP", "DsRedN", "DsRedP"),
    row.names = colnames(gs_scaled)
  )
  
  ann_colors <- list(
    Condition = c(WT="orange3", KO="aquamarine4"),
    SampleType = c(DsRedN="#4876FF", DsRedP="#CD4F39")
  )
  
  # -------------------------------
  # Crear carpeta GSVA si no existe
  # -------------------------------
  gsva_dir <- file.path(outdir, "GSVA")
  if(!dir.exists(gsva_dir)) dir.create(gsva_dir, recursive = TRUE)
  
  # -------------------------------
  # Generar heatmap
  # -------------------------------
  p <- pheatmap(
    gs_scaled,
    cluster_rows=TRUE,
    cluster_cols=FALSE,
    annotation_col=annotation_col,
    annotation_colors=ann_colors,
    color=rdbu_pal,
    main=paste0("GSVA Cytokine - ", cl),
    fontsize_row=6,
    border_color=NA
  )
  
  # -------------------------------
  # Guardar PDF
  # -------------------------------
  pdf(file.path(gsva_dir, paste0("GSVA_Cytokine_", cl, "_WTvsKO_normalized.pdf")), width=12, height=20)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grid::grid.text("GSVA Score", x=0.85, y=0.93, gp=grid::gpar(fontsize=10))
  grid::grid.text("GO Terms", x=0.88, y=0.5, rot=90, gp=grid::gpar(fontsize=10))
  dev.off()
  
  message("✔ PDF generado para ", cl)
}




lapply(names(cluster_files), function(cl) {
  plot_gsva_cluster_cytokine(
    cl = cl,
    cluster_file = cluster_files[[cl]],
    rsdir = rsdir,
    outdir = outdir
  )
})






plot_gsva_cluster_cytokine <- function(cl, cluster_file) {
  message("\n==============================")
  message("Procesando cluster: ", cl)
  message("==============================\n")
  
  # -------------------------------
  # Parámetros fijos
  # -------------------------------
  min_genes <- 30
  max_genes <- 500
  top_n <- 40
  rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  
  # -------------------------------
  # Leer matriz desde archivo
  # -------------------------------
  full_path <- file.path(rsdir, cluster_file)
  voom_mat <- read.table(full_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  
  # Convertir todo a numérico y reemplazar NAs
  voom_mat <- as.matrix(voom_mat)
  mode(voom_mat) <- "numeric"
  voom_mat[is.na(voom_mat)] <- 0
  
  # -------------------------------
  # Gene sets GO BP relacionados con citoquinas
  # -------------------------------
  msig_go <- msigdbr(species="Mus musculus", category="C5", subcategory="BP")
  
  include_keywords <- c(
    "cytokine", "interleukin", "IL-", "IL[0-9]",
    "chemokine", "CXCL", "CCL", "TNF"
  )
  exclude_keywords <- c("cytokinesis", "cytokinetic")
  
  cytokine_go <- msig_go %>%
    dplyr::filter(grepl(paste(include_keywords, collapse="|"), gs_name, ignore.case=TRUE)) %>%
    dplyr::filter(!grepl(paste(exclude_keywords, collapse="|"), gs_name, ignore.case=TRUE))
  
  cytokine_list <- split(cytokine_go$gene_symbol, cytokine_go$gs_name) %>%
    lapply(unique) %>%
    Filter(function(x) length(x) >= min_genes & length(x) <= max_genes, .)
  
  message("Pathways seleccionados: ", length(cytokine_list))
  
  # -------------------------------
  # GSVA
  # -------------------------------
  params <- gsvaParam(exprData = voom_mat, geneSets = cytokine_list, kcdf = "Gaussian")
  gs <- gsva(params)
  
  # -------------------------------
  # Seleccionar top pathways según promedio KO
  # -------------------------------
  KO_cols <- grep("KO", colnames(gs))
  KO_means <- rowMeans(gs[, KO_cols, drop=FALSE])
  KO_sorted <- sort(KO_means, decreasing=TRUE)
  top_pathways <- names(KO_sorted)[1:min(top_n, length(KO_sorted))]
  gs_top <- gs[top_pathways, , drop=FALSE]
  
  # Escalado (z-score)
  gs_scaled <- gs_top 
  
  # -------------------------------
  # Anotación automática de columnas
  # -------------------------------
  samples <- colnames(gs_scaled)
  Condition <- ifelse(grepl("KO", samples), "KO", "WT")
  SampleType <- ifelse(grepl("DsRedP", samples), "DsRedP", "DsRedN")
  annotation_col <- data.frame(Condition=Condition, SampleType=SampleType, row.names=samples)
  
  ann_colors <- list(
    Condition = c(WT="orange3", KO="aquamarine4"),
    SampleType = c(DsRedN="#4876FF", DsRedP="#CD4F39")
  )
  
  # -------------------------------
  # Crear carpeta GSVA si no existe
  # -------------------------------
  gsva_dir <- file.path(outdir, "GSVA")
  if(!dir.exists(gsva_dir)) dir.create(gsva_dir, recursive = TRUE)
  
  # -------------------------------
  # Generar heatmap
  # -------------------------------
  p <- pheatmap(
    gs_scaled,
    cluster_rows=TRUE,
    cluster_cols=FALSE,
    annotation_col=annotation_col,
    annotation_colors=ann_colors,
    color=rdbu_pal,
    main=paste0("GSVA Cytokine - ", cl),
    fontsize_row=6,
    border_color=NA
  )
  
  # -------------------------------
  # Guardar PDF
  # -------------------------------
  pdf(file.path(gsva_dir, paste0("GSVA_Cytokine_", cl, "_normalized.pdf")), width=14, height=14)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grid::grid.text("GSVA Score", x=0.85, y=0.93, gp=grid::gpar(fontsize=10))
  grid::grid.text("GO Terms", x=0.88, y=0.5, rot=90, gp=grid::gpar(fontsize=10))
  dev.off()
  
  message("✔ PDF generado para ", cl)
}


lapply(names(cluster_files), function(cl) {
  plot_gsva_cluster_cytokine(
    cl = cl,
    cluster_file = cluster_files[[cl]]  # archivo TSV relativo a rsdir
  )
})




plot_gsva_cluster_cytokine <- function(cl, cluster_file) {
  message("\n==============================")
  message("Procesando cluster: ", cl)
  message("==============================\n")
  
  # -------------------------------
  # Parámetros fijos
  # -------------------------------
  min_genes <- 30
  max_genes <- 500
  top_n <- 40
  rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))
  
  # -------------------------------
  # Leer matriz desde archivo
  # -------------------------------
  full_path <- file.path(rsdir, cluster_file)
  voom_mat <- read.table(full_path, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
  
  # Convertir todo a numérico y reemplazar NAs
  voom_mat <- as.matrix(voom_mat)
  mode(voom_mat) <- "numeric"
  voom_mat[is.na(voom_mat)] <- 0
  
  # -------------------------------
  # Gene sets GO BP relacionados con citoquinas
  # -------------------------------
  msig_go <- msigdbr(species="Mus musculus", category="C5", subcategory="BP")
  
  include_keywords <- c(
    "cytokine", "interleukin", "IL-", "IL[0-9]",
    "chemokine", "CXCL", "CCL", "TNF"
  )
  exclude_keywords <- c("cytokinesis", "cytokinetic")
  
  cytokine_go <- msig_go %>%
    dplyr::filter(grepl(paste(include_keywords, collapse="|"), gs_name, ignore.case=TRUE)) %>%
    dplyr::filter(!grepl(paste(exclude_keywords, collapse="|"), gs_name, ignore.case=TRUE))
  
  cytokine_list <- split(cytokine_go$gene_symbol, cytokine_go$gs_name) %>%
    lapply(unique) %>%
    Filter(function(x) length(x) >= min_genes & length(x) <= max_genes, .)
  
  message("Pathways seleccionados: ", length(cytokine_list))
  
  # -------------------------------
  # GSVA
  # -------------------------------
  params <- gsvaParam(exprData = voom_mat, geneSets = cytokine_list, kcdf = "Gaussian")
  gs <- gsva(params)
  
  # -------------------------------
  # Seleccionar top pathways según promedio KO
  # -------------------------------
  KO_cols <- grep("KO", colnames(gs))
  KO_means <- rowMeans(gs[, KO_cols, drop=FALSE])
  KO_sorted <- sort(KO_means, decreasing=TRUE)
  top_pathways <- names(KO_sorted)[1:min(top_n, length(KO_sorted))]
  gs_top <- gs[top_pathways, , drop=FALSE]
  
  # Escalado (z-score)
  gs_scaled <- gs_top 
  
  # -------------------------------
  # Anotación automática de columnas
  # -------------------------------
  samples <- colnames(gs_scaled)
  Condition <- ifelse(grepl("KO", samples), "KO", "WT")
  SampleType <- ifelse(grepl("DsRedP", samples), "DsRedP", "DsRedN")
  annotation_col <- data.frame(Condition=Condition, SampleType=SampleType, row.names=samples)
  
  ann_colors <- list(
    Condition = c(WT="orange3", KO="aquamarine4"),
    SampleType = c(DsRedN="#4876FF", DsRedP="#CD4F39")
  )
  
  # -------------------------------
  # Crear carpeta GSVA si no existe
  # -------------------------------
  gsva_dir <- file.path(outdir, "GSVA")
  if(!dir.exists(gsva_dir)) dir.create(gsva_dir, recursive = TRUE)
  
  # -------------------------------
  # Generar heatmap
  # -------------------------------
  p <- pheatmap(
    gs_scaled,
    cluster_rows=TRUE,
    cluster_cols=FALSE,
    annotation_col=annotation_col,
    annotation_colors=ann_colors,
    color=rdbu_pal,
    main=paste0("GSVA Cytokine - ", cl),
    fontsize_row=6,
    border_color=NA
  )
  
  # -------------------------------
  # Guardar PDF
  # -------------------------------
  pdf(file.path(gsva_dir, paste0("GSVA_Cytokine_", cl, "_normalized.pdf")), width=14, height=14)
  grid::grid.newpage()
  grid::grid.draw(p$gtable)
  grid::grid.text("GSVA Score", x=0.85, y=0.93, gp=grid::gpar(fontsize=10))
  grid::grid.text("GO Terms", x=0.88, y=0.5, rot=90, gp=grid::gpar(fontsize=10))
  dev.off()
  
  message("✔ PDF generado para ", cl)
}


lapply(names(cluster_files), function(cl) {
  plot_gsva_cluster_cytokine(
    cl = cl,
    cluster_file = cluster_files[[cl]]  # archivo TSV relativo a rsdir
  )
})





plot_gsva_cluster_cytokine <- function(cl, cluster_file) {
  message("\n==============================")
  message("Procesando cluster: ", cl)
  message("==============================\n")

  # Parámetros
  min_genes <- 30
  max_genes <- 500
  top_n <- 40
  rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))

  # 1. Leer matriz cruda
  m <- read.table(paste0(rsdir, cluster_file),
                  header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  m[is.na(m)] <- 0

  # 2. Obtener pathways GO BP
  msig_go <- msigdbr(
    species = "Mus musculus",
    category = "C5",
    subcategory = "BP"
  )

  include_keywords <- c(
    "cytokine",
    "interleukin", "IL-", "IL[0-9]",
    "chemokine", "CXCL", "CCL",
    "TNF"
  )

  exclude_keywords <- c("cytokinesis", "cytokinetic")

  cytokine_go <- msig_go %>%
    dplyr::filter(grepl(paste(include_keywords, collapse = "|"), gs_name, ignore.case = TRUE)) %>%
    dplyr::filter(!grepl(paste(exclude_keywords, collapse = "|"), gs_name, ignore.case = TRUE))

  cytokine_list <- split(cytokine_go$gene_symbol, cytokine_go$gs_name) %>%
    lapply(unique) %>%
    Filter(function(x) length(x) >= min_genes & length(x) <= max_genes, .)

  message("Pathways seleccionados: ", length(cytokine_list))

  # 3. Normalizar
  y <- DGEList(counts = m)
  y <- calcNormFactors(y)
  voom_mat <- voom(y)$E

  # 4. GSVA
  params <- gsvaParam(exprData = voom_mat, geneSets = cytokine_list, kcdf = "Gaussian")
  gs <- gsva(params)

  # 5. Selección top pathways según KO
  KO_cols <- grep("KO", colnames(gs))
  KO_means <- rowMeans(gs[, KO_cols, drop = FALSE])
  KO_sorted <- sort(KO_means, decreasing = TRUE)
  top_pathways <- names(KO_sorted)[1:min(top_n, length(KO_sorted))]
  gs_top <- gs[top_pathways, , drop = FALSE]

  # 6. Escalado
  gs_scaled <- t(scale(t(gs_top)))

  # 7. Orden columnas
  col_order <- c(
    "WT1.DsRedN", "WT2_DsRedN",
    "WT1.DsRedP", "WT2_DsRedP",
    "DsRedN.KO2", "DsRedN.KO3",
    "DsRedP.KO2", "DsRedP.KO3"
  )
  col_order <- col_order[col_order %in% colnames(gs_scaled)]

  annotation_col <- data.frame(
    Condition = c(rep("WT", 4), rep("KO", 4)),
    SampleType = rep(c("DsRedN", "DsRedP"), each = 2, times = 2)
  )
  rownames(annotation_col) <- col_order

  ann_colors <- list(
    Condition = c(WT = "orange3", KO = "aquamarine4"),
    SampleType = c(DsRedN = "#4876FF", DsRedP = "#CD4F39")
  )

  # 8. Heatmap
  p <- pheatmap(
    gs_scaled[, col_order, drop = FALSE],
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    annotation_col = annotation_col,
    annotation_colors = ann_colors,
    color = rdbu_pal,
    main = paste0("GSVA Cytokine - ", cl),
    fontsize_row = 6,
    border_color = NA
  )

  # 9. Guardar PDF
  pdf(file.path(outdir, paste0("/GSVA/GSVA_Cytokine_", cl, ".pdf")), width = 14, height = 14)
  grid.newpage()
  grid.draw(p$gtable)
  grid.text("GSVA Score", x = 0.85, y = 0.93, gp = gpar(fontsize = 10))
  grid.text("GO Terms", x = 0.88, y = 0.5, rot = 90, gp = gpar(fontsize = 10))
  dev.off()

  message("✔ PDF generado para ", cl)
}





lapply(names(cluster_files), function(cl) {
  plot_gsva_cluster_cytokine(
    cl = cl,
    cluster_file = cluster_files[[cl]]
  )
})





# -------------------------------
# Librerías necesarias
# -------------------------------
library(GSVA)
library(pheatmap)
library(msigdbr)
library(dplyr)
library(grid)
library(RColorBrewer)

# -------------------------------
# Parámetros generales
# -------------------------------
min_genes <- 5
max_genes <- 500

# -------------------------------
# 1. Preparar lista de pathways GO BP
# -------------------------------
# -------------------------------
# 1. Preparar lista de pathways GO BP (solo GO, sin Hallmark)
# -------------------------------
msig_go <- msigdbr(
  species = "Mus musculus",
  category = "C5",
  subcategory = "BP"
)

# Palabras clave para seleccionar pathways deseados
include_keywords <- c(
  "cytokine",                          # Citoquinas
  "interleukin", "IL-", "IL[0-9]",     # Interleukinas
  "chemokine", "CXCL", "CCL",          # Quimiocinas
  "TNF"                                 # Vías de TNF
)

# Palabras clave que quieres **excluir**
exclude_keywords <- c(
  "cytokinesis", "cytokinetic"          # Eliminado explícitamente
)

# 1) Seleccionar pathways relacionados con citoquinas/interleukinas/quimiocinas
cytokine_go <- msig_go %>%
  filter(grepl(paste(include_keywords, collapse="|"),
               gs_name, ignore.case = TRUE))

# 2) Eliminar los pathways que contengan "cytokinesis" o "citokinesis"
cytokine_go <- cytokine_go %>%
  filter(!grepl(paste(exclude_keywords, collapse="|"),
                gs_name, ignore.case = TRUE))

# 3) Convertir a lista para GSVA
cytokine_list <- split(cytokine_go$gene_symbol, cytokine_go$gs_name) %>%
  lapply(unique) %>%
  Filter(function(x) length(x) >= min_genes & length(x) <= max_genes, .)

message("Pathways GO seleccionados (excluyendo citokinesis): ", 
        length(cytokine_list))


# -------------------------------
# 2. Crear listas para almacenar GSVA medios por cluster
# -------------------------------
WT_list <- list()
KO_list <- list()

# -------------------------------
# 3. Loop sobre clusters
# -------------------------------
for (cl in names(cluster_files)) {

  message("Cluster: ", cl)

  m <- read.table(
    paste0(rsdir, cluster_files[[cl]]),
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  )

  WT_cols <- grep("WT", colnames(m))
  KO_cols <- grep("KO", colnames(m))

  if (length(WT_cols) == 0 & length(KO_cols) == 0) {
    warning("No se encontraron columnas WT o KO en el cluster ", cl, ". Se omite.")
    next
  }

  go_list2 <- lapply(cytokine_list, function(x) intersect(x, rownames(m)))
  go_list2 <- go_list2[sapply(go_list2, length) >= min_genes]

  message("Número de pathways de citoquinas usados en este cluster: ", length(go_list2))

  params <- gsvaParam(
    exprData = as.matrix(m),
    geneSets = go_list2,
    kcdf     = "Gaussian",
    minSize  = 1,
    maxSize  = 500
  )

  gs <- gsva(params, verbose = FALSE)

  if (length(WT_cols) > 0) {
    WT_list[[cl]] <- rowMeans(gs[, WT_cols, drop = FALSE])
  } else {
    WT_list[[cl]] <- NA
  }

  if (length(KO_cols) > 0) {
    KO_list[[cl]] <- rowMeans(gs[, KO_cols, drop = FALSE])
  } else {
    KO_list[[cl]] <- NA
  }
}

# -------------------------------
# 4. Convertir a matrices
# -------------------------------
GSVA_WT <- do.call(cbind, WT_list)
GSVA_KO <- do.call(cbind, KO_list)

colnames(GSVA_WT) <- paste0(names(WT_list), "_WT")
colnames(GSVA_KO) <- paste0(names(KO_list), "_KO")

# -------------------------------
# 5. Unir WT + KO
# -------------------------------
GSVA_combined <- cbind(GSVA_WT, GSVA_KO)

# -------------------------------
# 6. Anotaciones de columnas
# -------------------------------
annotation_col <- data.frame(
  Condition = c(rep("WT", ncol(GSVA_WT)), rep("KO", ncol(GSVA_KO))),
  Cluster   = rep(names(cluster_files), 2)
)
rownames(annotation_col) <- colnames(GSVA_combined)

nora.colors2 <- c(
  "Ly6cHi_Monocytes"                = "#FF0000",
  "Ly6cLo_Monocytes"                = "#FF6A6A",
  "Early_IFN_MHCII_TAMs"            = "#B22222",
  "Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac" = "#4DAF4A",
  "Trem1_Ptgs2_Plaur_Celc4e_Mac"    = "#EE7942",
  "MHCII_Ccl12_Mac"                 = "#4682B4",
  "MHCII_Siglec_Mac"                = "#1E90FF",
  "IFN_Mac"                          = "#C080FF",
  "Mmp9_Ctsk_Mac"                    = "#00723F",
  "Mrc1_C1qc_Cbr2_Gas6_Mac"         = "#FFD92F",
  "Npr2_Actn1_Mac"                   = "#A6D854",
  "Fn1_Vegfa_Mac"                    = "#FFA500",
  "Neutrophils"                      = "#4876FF"
)

annotation_colors <- list(
  Condition = c(WT = "#1F78B4", KO = "#E31A1C"),
  Cluster   = nora.colors2
)

# -------------------------------
# 7. Top 50 pathways por KO–WT
# -------------------------------
FC <- rowMeans(GSVA_KO, na.rm = TRUE) - rowMeans(GSVA_WT, na.rm = TRUE)
top_cytokine_pathways <- names(sort(FC, decreasing = TRUE))[1:50]

GSVA_combined_top <- GSVA_combined[top_cytokine_pathways, ]

# -------------------------------
# 8. Clustering de columnas CON Pearson + ward.D2
# -------------------------------
KO_cols_idx <- grep("_KO$", colnames(GSVA_combined_top))

# correlación entre columnas KO
cor_cols <- cor(GSVA_combined_top[, KO_cols_idx], method = "pearson")

dist_cols <- as.dist(1 - cor_cols)

hc_cols <- hclust(dist_cols, method = "ward.D2")

KO_order <- KO_cols_idx[hc_cols$order]
KO_colnames_ordered <- colnames(GSVA_combined_top)[KO_order]

WT_colnames_ordered <- sub("_KO$", "_WT", KO_colnames_ordered)

combined_col_order <- c(WT_colnames_ordered, KO_colnames_ordered)

# -------------------------------
# 9. Clustering de FILAS (pathways) con Pearson + ward.D2
# -------------------------------
cor_rows <- cor(t(GSVA_combined_top), method = "pearson")
dist_rows <- as.dist(1 - cor_rows)
hc_rows <- hclust(dist_rows, method = "ward.D2")

row_order <- hc_rows$order

# -------------------------------
# 10. Paleta de color
# -------------------------------
rdbu_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(100))

# -------------------------------
# 11. PDF final
# -------------------------------
pdf(file = file.path(outdir, "GSVA", "GSVA_WT_vs_KO_cytokines_top50_PEARSON_WARD_rows_cols.pdf"),
    width = 20, height = 12)

pheatmap(
  GSVA_combined_top[row_order, combined_col_order],
  main = "GSVA - Top Cytokines 50 pathways",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_col[combined_col_order, ],
  annotation_colors = annotation_colors,
  fontsize_row = 6,
  show_colnames = FALSE,
  border_color = NA,
  gaps_col = length(WT_colnames_ordered),
  color = rdbu_pal,
  scale = "row"
)

grid.text("GSVA score",
          x = 0.80, y = 0.89,
          gp = gpar(fontsize = 12))

dev.off()





library(GSVA)
library(pheatmap)
library(msigdbr)
library(dplyr)
library(grid)

# -------------------------------
# Preparar lista de pathways GO BP
# -------------------------------
msig_go <- msigdbr(species = "Mus musculus", category = "C5", subcategory = "BP")

# Filtrar por tamaño de pathway
go_list_all <- split(msig_go$gene_symbol, msig_go$gs_name) %>%
  lapply(unique) %>%
  Filter(function(x) length(x) >= min_genes & length(x) <= max_genes, .)

# -------------------------------
# 1. Crear listas para almacenar GSVA medios por cluster
# -------------------------------
WT_list <- list()
KO_list <- list()

for (cl in names(cluster_files)) {

  message("Cluster: ", cl)

  # Leer pseudobulk del cluster
  m <- read.table(
    paste0(rsdir, cluster_files[[cl]]),
    header = TRUE,
    row.names = 1,
    sep = "\t",
    check.names = FALSE
  )

  WT_cols <- grep("WT", colnames(m))
  KO_cols <- grep("KO", colnames(m))

  if (length(WT_cols) == 0 & length(KO_cols) == 0) next

  # Filtrar pathways según genes presentes en el cluster
  go_list2 <- lapply(go_list_all, function(x) intersect(x, rownames(m)))
  go_list2 <- go_list2[sapply(go_list2, length) >= min_genes]

  # GSVA por cluster
  params <- gsvaParam(
    exprData = as.matrix(m),
    geneSets = go_list2,
    kcdf     = "Gaussian",
    minSize  = 1,
    maxSize  = 500
  )
  gs <- gsva(params, verbose = FALSE)

  # Promedios por condición
  if (length(WT_cols) > 0) WT_list[[cl]] <- rowMeans(gs[, WT_cols, drop = FALSE])
  if (length(KO_cols) > 0) KO_list[[cl]] <- rowMeans(gs[, KO_cols, drop = FALSE])
}

GSVA_WT <- do.call(cbind, WT_list)
GSVA_KO <- do.call(cbind, KO_list)

colnames(GSVA_WT) <- paste0(names(WT_list), "_WT")
colnames(GSVA_KO) <- paste0(names(KO_list), "_KO")

# -------------------------------
# 2. Matriz de diferencias KO-WT
# -------------------------------
GSVA_delta <- GSVA_KO - GSVA_WT

# -------------------------------
# 3. Seleccionar top pathways según promedio KO
# -------------------------------
KO_means <- rowMeans(GSVA_KO, na.rm = TRUE)
top_pathways <- names(sort(KO_means, decreasing = TRUE))[1:50]

# Filtrar matrices
GSVA_WT_top <- GSVA_WT[top_pathways, ]
GSVA_KO_top <- GSVA_KO[top_pathways, ]
GSVA_delta_top <- GSVA_delta[top_pathways, ]

# -------------------------------
# 4. Crear heatmap combinado (WT + KO) con anotación de Delta
# -------------------------------
GSVA_combined_top <- cbind(GSVA_WT_top, GSVA_KO_top)

annotation_col <- data.frame(
  Condition = c(rep("WT", ncol(GSVA_WT_top)), rep("KO", ncol(GSVA_KO_top))),
  Cluster   = rep(names(cluster_files), 2)
)
rownames(annotation_col) <- colnames(GSVA_combined_top)

# Colores de clusters
nora.colors2 <- c(
  "Ly6cHi_Monocytes"               = "#FF0000",
  "Ly6cLo_Monocytes"               = "#FF6A6A",
  "Early_IFN_MHCII_TAMs"           = "#B22222",
  "Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac" = "#4DAF4A",
  "Trem1_Ptgs2_Plaur_Celc4e_Mac"   = "#EE7942",
  "MHCII_Ccl12_Mac"                = "#4682B4",
  "MHCII_Siglec_Mac"               = "#1E90FF",
  "IFN_Mac"                         = "#C080FF",
  "Mmp9_Ctsk_Mac"                   = "#00723F",
  "Mrc1_C1qc_Cbr2_Gas6_Mac"        = "#FFD92F",
  "Npr2_Actn1_Mac"                  = "#A6D854",
  "Fn1_Vegfa_Mac"                   = "#FFA500",
  "Neutrophils"                     = "#4876FF"
)

annotation_colors <- list(
  Condition = c(WT = "#1F78B4", KO = "#E31A1C"),
  Cluster   = nora.colors2
)

# -------------------------------
# 5. Heatmap WT + KO
# -------------------------------
pdf(file = file.path(outdir, "GSVA", "GSVA_WT_vs_KO_byCluster_top50.pdf"),
    width = 20, height = 12)

pheatmap(
  GSVA_combined_top,
  main = "GSVA - Top 50 pathways by Tet2 KO average",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 6,
  show_colnames = FALSE,
  border_color = NA,
  gaps_col = ncol(GSVA_WT_top)  # línea separadora entre WT y KO
)

# Añadir texto sobre la barra de color
grid.text("GSVA score", x = 0.79, y = 0.93, gp = gpar(fontsize = 12))
dev.off()

# -------------------------------
# 6. Heatmap de diferencias KO-WT (opcional)
# -------------------------------
pdf(file = file.path(outdir, "GSVA", "GSVA_delta_KO_minus_WT_top50.pdf"),
    width = 15, height = 12)

pheatmap(
  GSVA_delta_top,
  main = "GSVA Delta (KO - WT) - Top 50 pathways",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  fontsize_row = 6,
  show_colnames = FALSE,
  border_color = NA,
  gaps_col = ncol(GSVA_WT_top)
)
dev.off()
