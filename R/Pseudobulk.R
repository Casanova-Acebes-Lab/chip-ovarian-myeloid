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


run_pseudobulk_analysis <- function(data, cluster_name, outdir, rsdir) {
  # Ejecutar pseudobulk para un cluster específico
  res <- pseudobulk_cluster(data, cluster_name = cluster_name)
  pb_mat <- res$pb_mat
  n_cells <- res$n_cells
  n_counts <- res$n_counts

  # Crear metadatos
  sample_meta <- make_sample_metadata(colnames(pb_mat))

  # Ejecutar limma
  res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)

  # Clasificación de genes diferencialmente expresados
  res_limma$diffexpressed <- "NO"
  res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
  res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"

  # Top genes
  top_genes <- res_limma %>%
    arrange(adj.P.Val) %>%
    slice_head(n = 50) %>%
    as.data.frame()

  # Guardar tabla completa
  write.table(
    res_limma,
    file = paste0(rsdir, "/table.macros.", cluster_name, ".tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Crear Volcano plot
  p <- Volcano2(data = res_limma, legend = paste0("KO vs WT samples. ", cluster_name))

  # Guardar el plot en PDF
  pdf(paste0(outdir, "/Pseudobulk/Volcano.KO_vs_WT_", cluster_name, ".pdf"),
      width = 16, height = 12)
  print(p)
  dev.off()

  # Devolver resultados
  return(list(
    cluster = cluster_name,
    res_limma = res_limma,
    top_genes = top_genes,
    n_cells = n_cells,
    n_counts = n_counts,
    plot = p
  ))
}




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



# Ly6c|Ms4ac Monocytes

res <- pseudobulk_cluster(data, cluster_name="Ly6cHi Monocytes")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)


p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Ly6cHi Monocytes")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Ly6cHi Monocytes.pdf"), width=16, height=12)
print(p)
dev.off()



write.table(res_limma, paste0(rsdir,"table.macros.Ly6cHi Monocytes.tsv"), sep='\t')


# Trem1|Ptgs2|Plaur|Celc4e


res <- pseudobulk_cluster(data, cluster_name="Trem1|Ptgs2|Plaur|Celc4e Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples.Trem1|Ptgs2|Plaur|Celc4e Mac")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Trem1|Ptgs2|Plaur|Celc4e Mac.pdf"), width=16, height=12)
print(p)
dev.off()

write.table(res_limma, paste0(rsdir,"table.macros.Trem1|Ptgs2|Plaur|Celc4e Mac.tsv"), sep='\t')




# Mrc1|C1Pseudobulk|Cbr2|Gas6 Mac


res <- pseudobulk_cluster(data, cluster_name="Mrc1|C1Pseudobulk|Cbr2|Gas6 Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples.Mrc1|C1Pseudobulk|Cbr2|Gas6 Mac")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Mrc1|C1Pseudobulk|Cbr2|Gas6 Mac.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Mrc1|C1Pseudobulk|Cbr2|Gas6 Mac.tsv"), sep='\t')






# Arg1|Spp1|Mmp12|Mmp19|Il1a Mac


res <- pseudobulk_cluster(data, cluster_name="Arg1|Spp1|Mmp12|Mmp19|Il1a Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)




p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Arg1|Spp1|Mmp12|Mmp19|Il1a Mac")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Arg1|Spp1|Mmp12|Mmp19|Il1a Mac.pdf"), width=16, height=12)
print(p)
dev.off()

write.table(res_limma, paste0(rsdir,"table.macros.Arg1|Spp1|Mmp12|Mmp19|Il1a Mac.tsv"), sep='\t')




# Npr2|Actn1 Mac


res <- pseudobulk_cluster(data, cluster_name="Npr2|Actn1 Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Npr2|Actn1 Mac")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Npr2|Actn1 Mac.pdf"), width=16, height=12)
print(p)
dev.off()



write.table(res_limma, paste0(rsdir,"table.macros.Npr2|Actn1 Mac.tsv"), sep='\t')




# Mmp9|Ctsk Mac


res <- pseudobulk_cluster(data, cluster_name="Mmp9|Ctsk Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)


p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Mmp9|Ctsk Mac")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Mmp9|Ctsk Mac.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Mmp9|Ctsk Mac.tsv"), sep='\t')




# Ciita|Siglec Mac 


res <- pseudobulk_cluster(data, cluster_name="Ciita|Siglec Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Ciita|Siglec Mac")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Ciita|Siglec Mac.pdf"), width=16, height=12)
print(p)
dev.off()



write.table(res_limma, paste0(rsdir,"table.macros.Ciita|Siglec Mac.tsv"), sep='\t')




 # Ciita|Ccl12 Mac

res <- pseudobulk_cluster(data, cluster_name="Ciita|Ccl12 Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Ciita|Ccl12 Mac")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Ciita|Ccl12 Mac.pdf"), width=16, height=12)
print(p)
dev.off()




write.table(res_limma, paste0(rsdir,"table.macros.Ciita|Ccl12 Mac.tsv"), sep='\t')


# IFN Mac 


res <- pseudobulk_cluster(data, cluster_name="IFN Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. IFN Mac")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.IFN Mac.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.IFN Mac.tsv"), sep='\t')




# Slamf Mac 


res <- pseudobulk_cluster(data, cluster_name="Slamf Monocytes")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Slamf Monocytes")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Slamf Monocytes.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Slamf Monocytes.tsv"), sep='\t')






# Fn1 Mac 


res <- pseudobulk_cluster(data, cluster_name="Fn1|Vegfa Mac")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Fn1|Vegfa Mac")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Fn1|Vegfa Mac.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Fn1|Vegfa Mac.tsv"), sep='\t')




# Neutrophils 


res <- pseudobulk_cluster(data, cluster_name="Neutrophils")
pb_mat <- res$pb_mat
n_cells <- res$n_cells
n_counts <- res$n_counts

sample_meta <- make_sample_metadata(colnames(pb_mat))

res_limma <- run_limma_cluster(pb_mat, sample_meta, n_cells)


res_limma$diffexpressed <- "NO"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC > 0.5] <- "Up"
res_limma$diffexpressed[res_limma$adj.P.Val < 0.05 & res_limma$logFC < -0.5] <- "Down"


# Top 20 genes
top_genes <- res_limma %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 50) %>%
  as.data.frame()

top_genes
table(res_limma$diffexpressed)



p <- Volcano2(data = res_limma, legend = "KO vs WT samples. Neutrophils")


pdf(paste0(outdir,"/Pseudobulk/Volcano.KO vs WT Macrophages.Neutrophils.pdf"), width=16, height=12)
print(p)
dev.off()


write.table(res_limma, paste0(rsdir,"table.macros.Neutrophils.tsv"), sep='\t')







# All macros clusters 

library(Matrix)
library(limma)
library(edgeR)
library(Matrix.utils)
library(dplyr)

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
# 🚀 Ejemplo: todos los clusters de macrófagos
# ---------------------------


# Lista de clusters
macro_clusters <- c(
  "Ly6c|Ms4ac Monocytes",
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac",
  "Ciita|Ccl12 Mac",
  "Ciita|Siglec Mac",
  "Fn1|Vegfa Mac",
  "IFN Mac",
  "Mmp9|Ctsk Mac",
  "Mrc1|C1Pseudobulk|Cbr2|Gas6 Mac",
  "Neutrophils",
  "Npr2|Actn1 Mac",
  "Slamf Monocytes",
  "Trem1|Ptgs2|Plaur|Celc4e Mac"
)

pb_out <- pseudobulk_group(data, clusters = macro_clusters, cluster_col="Cluster", sample_col="tag")
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

