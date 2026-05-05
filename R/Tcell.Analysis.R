
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

data <- readRDS(paste0(rsdir,"tcells.clusterized.rds"))






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
  subset = Subcluster%in% c(
    "NK",
    "ILC"
  ), invert=TRUE
)


## Pseudobulk


pb <- AggregateExpression(
  data,
  assays = "RNA",   # usa RNA, no SCT
  slot = "counts",
  group.by = "tag"
)


# Counts
pb_counts <- pb$RNA
# Redondear a enteros
pb_counts <- round(pb_counts)

# Opcional: asegurarte de que sean numeric
pb_counts <- as.matrix(pb_counts)



# Metadata
meta <- data@meta.data

sample_meta <- meta[, c("tag", "group", "DsRed")] |>
  unique()

rownames(sample_meta) <- sample_meta$tag

sample_meta$tag <- gsub("_", "-", sample_meta$tag)
rownames(sample_meta) <- sample_meta$tag

sample_meta <- sample_meta[colnames(pb_counts), ]




# DESed2 Object 



dds <- DESeqDataSetFromMatrix(
  countData = pb_counts,
  colData = sample_meta,
  design = ~ group
)

dds <- dds[rowSums(counts(dds)) > 20, ]
dds <- DESeq(dds)



# PCA 

# Sumar counts y convertir a CPM (counts por millón) para normalizar tamaño de librería
library(edgeR)

# Crear objeto DGEList
dge <- DGEList(counts = pb_counts)
dge <- calcNormFactors(dge)  # TMM normalization

# CPM log-transformado
log_cpm <- cpm(dge, log=TRUE, prior.count=1)


var_genes <- apply(log_cpm, 1, var)
top_genes <- names(sort(var_genes, decreasing = TRUE))[1:2000]
pca_res <- prcomp(t(log_cpm[top_genes, ]))
group = sample_meta$group

library(ggplot2)

# Preparar dataframe para plot
pca_df <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  group = sample_meta$group,
  tag = sample_meta$tag
)
pdf(paste0(outdir,"/Tcell.Analysis/PCA.pdf"), width=16, height=16)
ggplot(pca_df, aes(x=PC1, y=PC2, color=group, label=tag)) +
  geom_point(size=4) +
  geom_text(vjust=1.5, hjust=1.1) +
  theme_classic() +
  labs(title="PCA de muestras pseudobulk", x="PC1", y="PC2")
dev.off()



# -------------------------------
# Library sizes (suma de counts por muestra)
# -------------------------------
library_sizes <- colSums(pb_counts)

# -------------------------------
# Data.frame para graficar
# -------------------------------
lib_df <- data.frame(
  Sample = names(library_sizes),
  LibrarySize = library_sizes,
  Group = sample_meta$group,
  DsRed = sample_meta$DsRed
)

# -------------------------------
# Plot
# -------------------------------
pdf(paste0(outdir,"/Tcell.Analysis/librarySize.pdf"), width=16, height=16)
ggplot(lib_df, aes(x=Sample, y=LibrarySize, fill=DsRed)) +
  geom_bar(stat="identity") +
  geom_text(aes(label=round(LibrarySize,0)), vjust=-0.5, size=3) +
  theme_classic() +
  labs(title="Library size por muestra (pseudobulk counts)",
       y="Library size (total counts)",
       x="Muestra") +
  scale_fill_manual(values=c("steelblue","tomato")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
  dev.off()




res <- results(dds, contrast = c("group", "KO", "WT"))
head(as.data.frame(res[order(res$padj), ]), n=100)

result <- as.data.frame(res)


library(dplyr)

result_ordered <- result %>%
  filter(!is.na(padj)) %>%
  arrange(padj)

head(result_ordered, 20)




  # Clasificar DEGs
result_ordered$diffexpressed <- "NO"
result_ordered$diffexpressed[result_ordered$padj < 0.05 & result_ordered$log2FoldChange > 1] <- "Up"
result_ordered$diffexpressed[result_ordered$padj < 0.05 & result_ordered$log2FoldChange < -1] <- "Down"

write.table(result_ordered, paste0(rsdir,"table.DEG.Tcells.KOvsWT.tsv"), sep='\t')



Volcano3 <- function(data, legend, logFC_threshold = 1, top_n = 40) {


  data$gene <- rownames(data)
  data <- data[!is.na(data$padj), ]

  # Clasificar DEGs
  data$diffexpressed <- "NO"
  data$diffexpressed[data$padj < 0.05 & data$log2FoldChange > logFC_threshold] <- "Up"
  data$diffexpressed[data$padj < 0.05 & data$log2FoldChange < -logFC_threshold] <- "Down"

  # Seleccionar los top_n genes por bando para etiquetar
  top_up <- data %>%
    filter(diffexpressed == "Up") %>%
    arrange(padj) %>%
    head(top_n)

  top_down <- data %>%
    filter(diffexpressed == "Down") %>%
    arrange(padj) %>%
    head(top_n)

  top_genes <- rbind(top_up, top_down)

  # Añadir columna de label
  data$label <- ifelse(data$gene %in% top_genes$gene, data$gene, NA)

  # Volcano plot
  plot <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj), color = diffexpressed)) +
    geom_point() +
    geom_text_repel(
  data = subset(data, !is.na(label)),
  aes(label = label),
  size = 3.5,               # tamaño de letra
  segment.color = "grey50", # color de líneas que conectan puntos con labels
  box.padding = 0.8,        # espacio alrededor de cada label
  point.padding = 0.5,      # distancia del punto
  max.overlaps = Inf,        # fuerza máxima para mostrar labels
  nudge_y = 0.5             # desplaza ligeramente los labels verticalmente
) +
    geom_vline(xintercept = c(-logFC_threshold, logFC_threshold),
               linetype = "dashed", color = "black") +
    theme_minimal() +
    scale_color_manual(values = c("Down" = "dodgerblue4",
                                  "NO"   = "dimgrey",
                                  "Up"   = "brown1")) +
    ggtitle(paste("DEG Volcano:", legend)) +
    labs(x = "log2FC", y = "-log10(adj.pval)")

  return(plot)
}

volcano <-Volcano3(result, "DEG. KOvsWT Tet2 Chimeras")
pdf(paste0(outdir,"/Tcell.Analysis/Vocano.KO-WT.Tcells.pdf"), width=22, height=12)
print(volcano)
dev.off()






GSEA.pseudo <- function(data) {
 
data <- result
  
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
  

  data <- data[!is.na(data$log2FoldChange) & !is.na(data$padj), ]
  # Reemplazar p_val == 0 por un valor muy pequeño para evitar -Inf
  epsilon <- 1e-300
  data$padj[data$padj == 0] <- epsilon
  
  # Calcular la métrica de ranking
  data$metric <- data$log2FoldChange * -log10(data$padj + 1e-8)

  
  # Ordenar por la métrica de ranking
  data <- data[order(data$metric, decreasing = TRUE), ]
  
  # Crear el vector nombrado para GSEA
  gene_metric <- data$metric
  names(gene_metric) <- rownames(data)
  
  # Ejecutar GSEA (ontología biológica por defecto)
  gseGO_result <- gseGO(geneList = gene_metric,
                        ont = "BP",
                        OrgDb = org.Mm.eg.db,
                        minGSSize = 25,
                        maxGSSize = 500,
                        eps = 1e-20,
                        nPermSimple = 10000,
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.01,
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



GSEA <- GSEA.pseudo(result)
GSEA@result




# =========================
# 0. Librerías
# =========================
library(clusterProfiler)
library(msigdbr)
library(dplyr)
library(enrichplot)

# =========================
# 1. Limpiar datos
# =========================
# Quitar NA en pvalue o logFC
res_clean <- res[!is.na(res$pvalue) & !is.na(res$log2FoldChange), ]

# Quitar duplicados si los hay
res_clean <- res_clean[!duplicated(rownames(res_clean)), ]

# =========================
# 2. Descargar MSigDB Hallmark
# =========================
msig <- msigdbr(species = "Mus musculus", category = "H")

# Si es humano:
# msig <- msigdbr(species = "Homo sapiens", category = "H")

term2gene <- msig[, c("gs_name", "gene_symbol")]

# =========================
# 3. Ranking 1: logFC
# =========================
gene_list_fc <- res_clean$log2FoldChange
names(gene_list_fc) <- rownames(res_clean)

gene_list_fc <- sort(gene_list_fc, decreasing = TRUE)

# =========================
# 4. Ranking 2: sign * -log10(pvalue)
# =========================
gene_list_pval <- sign(res_clean$log2FoldChange) * -log10(res_clean$pvalue)
names(gene_list_pval) <- rownames(res_clean)

gene_list_pval <- sort(gene_list_pval, decreasing = TRUE)

# =========================
# 5. GSEA con logFC
# =========================
gsea_fc <- clusterProfiler::GSEA(
  geneList = gene_list_fc,
  TERM2GENE = term2gene,
  pvalueCutoff = 0.05,
  verbose = FALSE
)
# =========================
# 6. GSEA con ranking combinado
# =========================
gsea_pval <- clusterProfiler::GSEA(
  geneList = gene_list_pval,
  TERM2GENE = term2gene,
  pvalueCutoff = 0.05,
  verbose = FALSE,
)

# =========================
# 7. Visualización
# =========================

# Dotplot
pdf(paste0(outdir,"/Tcell.Analysis/GSEA.KO-WT.Tcells.Hallmark.LogFC2.pdf"), width=12, height=12)
dotplot(gsea_fc, showCategory = 20, title = "GSEA Hallmark (logFC)")
dev.off()
dotplot(gsea_pval, showCategory = 20, title = "GSEA Hallmark (sign * -log10 pvalue)")

# Ridgeplot
ridgeplot(gsea_fc)
ridgeplot(gsea_pval)

# =========================
# 8. Ver resultados
# =========================
head(as.data.frame(gsea_fc))
head(as.data.frame(gsea_pval))

# =========================
# 9. Guardar resultados
# =========================
write.csv(as.data.frame(gsea_fc), "GSEA_Hallmark_logFC.csv")
write.csv(as.data.frame(gsea_pval), "GSEA_Hallmark_pval.csv")


write.table(gsea_fc, paste0(rsdir,"table.GSEA.Hallmark.Tcells.KOvsWT.tsv"), sep='\t')


library(ggplot2)
library(dplyr)

# =========================
# 1. Preparar datos
# =========================
df <- as.data.frame(gsea_fc@result)

# Limpiar nombres
df$Description <- gsub("HALLMARK_", "", df$Description)
df$Description <- gsub("_", " ", df$Description)

# Filtrar significativos (opcional pero recomendado)
df <- df %>% 
  filter(p.adjust < 0.05)

# Ordenar por NES
df <- df %>% 
  arrange(NES)

# Definir grupo
df$group <- ifelse(df$NES > 0, "KO", "WT")

# Significancia
df$signif <- ifelse(df$p.adjust < 0.001, "***",
              ifelse(df$p.adjust < 0.01, "**",
              ifelse(df$p.adjust < 0.05, "*", "")))

# Factor ordenado
df$Description <- factor(df$Description, levels = df$Description)

# =========================
# 2. Plot a PDF
# =========================
pdf(paste0(outdir,"/Tcell.Analysis/GSEA.KO-WT.Tcells.Hallmark2.pdf"),
    width=12, height=8)

ggplot(df, aes(x = NES, y = Description, fill = group)) +
  
  # Barras
  geom_col(width = 0.7) +
  
  # Línea central
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.8) +
  
  # Significancia
  geom_text(aes(label = signif),
            hjust = ifelse(df$NES > 0, -0.3, 1.2),
            size = 5) +
  
  # Colores PRO
  scale_fill_manual(values = c(
    "WT" = "#2563EB",   # azul
    "KO" = "#DC2626"    # rojo
  )) +
  
  # Expandir límites para que quepan los asteriscos
  expand_limits(x = c(min(df$NES) - 0.5, max(df$NES) + 0.5)) +
  
  # Tema limpio tipo paper
  theme_minimal(base_size = 16) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.text.x = element_text(size = 12),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5)
  ) +
  
  labs(
    x = "Normalized Enrichment Score (NES)",
    y = "",
    title = "GSEA Hallmark: KO vs WT T cells"
  )

dev.off()




### Hipoxia genes 


vsd <- vst(dds, blind = FALSE)

expr_matrix <- assay(vsd)

hypoxia_genes <- c(
  "Mxi1","Hs3st1","Xpnpep1","Chst3","Slc2a3","Nfil3","Noct","Ldhc",
  "Sult2b1","Ets1","S100a4","Ero1a","Edn2","Adm","Kdm3a","Ndrg1",
  "Cdkn1c","Irs2","Klf7","Pnrc1","Prkca","Rbpj","Fosl2","Akap12",
  "Siah2","Ids","Sdc4","Wsb1","Chst2","Ppp1r15a","Map3k1","Cdkn1a",
  "Bnip3l","Nr3c1","Hk2","Pim1","Cited2","Pygm","Tnfaip3","Bhlhe40",
  "Errfi1","Ext1","Rora"
)


genes_use <- intersect(rownames(expr_matrix), hypoxia_genes)

expr_sub <- expr_matrix[genes_use, ]


hypoxia_score <- colMeans(expr_sub)

df <- data.frame(
  sample = names(hypoxia_score),
  score = hypoxia_score
) %>%
  left_join(coldata, by = "sample")


expr_scaled <- t(scale(t(expr_sub)))

hypoxia_score <- colMeans(expr_scaled)




library(ggplot2)
library(ggpubr)


pdf(paste0(outdir,"/Tcell.Analysis/Hipoxia.Score.pdf"),
    width=12, height=8)
ggplot(df, aes(x = group, y = score, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  
  stat_compare_means(
    method = "wilcox.test", 
    label = "p.format"   # o "p.signif" para ****
  ) +
  
  scale_fill_manual(values = c("WT" = "#2563EB", "KO" = "#DC2626")) +
  
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    title = "Hypoxia signature (pseudobulk)",
    x = "",
    y = "VST expression"
  )
  dev.off()



  ###

library(tidyverse)
library(SummarizedExperiment)  # para trabajar con colData

# Convertir colData a data.frame
coldata <- as.data.frame(dds@colData)
coldata$sample <- coldata$tag  # tu variable de sample

# Subconjunto de genes
genes_use <- intersect(rownames(expr_matrix), hypoxia_genes)
expr_sub <- expr_matrix[genes_use, ]

# Transformar a long format
df_long <- as.data.frame(expr_sub) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  left_join(coldata, by = "sample")  # ahora funciona porque ambos son data.frame

# Asegurar orden de condiciones
df_long$group <- factor(df_long$group, levels = c("WT", "KO"))


library(ggpubr)
library(dplyr)

# Comparación entre grupos
my_comparisons <- list(c("WT", "KO"))

# Calcular el máximo de expresión por gen para colocar los asteriscos
y_positions <- df_long %>%
  group_by(gene) %>%
  summarise(y.position = max(expression, na.rm = TRUE) * 1.1)

# Añadir la posición de y a df_long por gen
df_long <- df_long %>%
  left_join(y_positions, by = "gene")

# Boxplot faceteado por gen
p <- ggboxplot(
  df_long,
  x = "group",
  y = "expression",
  color = "group",
  palette = c("blue", "brown1"),
  add = "jitter"
) +
  stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    method = "wilcox.test",
    exact = FALSE,
    aes(y = y.position)  # ahora y.position está en df_long
  ) +
  facet_wrap(~ gene, scales = "free_y", ncol = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

# Guardar en PDF
pdf(paste0(outdir, "/Tcell.Analysis/Hipoxia.Score.genes.pdf"),
    width = 18, height = 22)
print(p)
dev.off()




### E2F Targets


library(tidyverse)
library(SummarizedExperiment)
library(ggpubr)

# Lista de genes E2F targets
e2f_genes <- c(
  "Fas","Fasl"
)

# Asegurar colData como data.frame
coldata <- as.data.frame(dds@colData)
coldata$sample <- coldata$tag

# Subconjunto de genes
genes_use <- intersect(rownames(expr_matrix), e2f_genes)
expr_sub <- expr_matrix[genes_use, ]

# Calcular "E2F score" promedio por sample
expr_scaled <- t(scale(t(expr_sub)))  # opcional: escalar por gene
e2f_score <- colMeans(expr_scaled)

df <- data.frame(
  sample = names(e2f_score),
  score = e2f_score
) %>%
  left_join(coldata, by = "sample")

# Boxplot general del score
pdf(paste0(outdir, "/Tcell.Analysis/E2F.Score.pdf"),
    width = 12, height = 8)
ggplot(df, aes(x = group, y = score, fill = group)) +
  geom_boxplot(width = 0.6, outlier.shape = NA) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +
  stat_compare_means(method = "wilcox.test", label = "p.format") +
  scale_fill_manual(values = c("WT" = "#2563EB", "KO" = "#DC2626")) +
  theme_minimal(base_size = 16) +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "E2F targets signature (pseudobulk)",
    x = "",
    y = "Scaled VST expression"
  )
dev.off()

# Transformar a long format para faceteado por gen
df_long <- as.data.frame(expr_sub) %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expression") %>%
  left_join(coldata, by = "sample")

df_long$group <- factor(df_long$group, levels = c("WT", "KO"))

# Calcular posición de los asteriscos
my_comparisons <- list(c("WT","KO"))
y_positions <- df_long %>%
  group_by(gene) %>%
  summarise(y.position = max(expression, na.rm = TRUE) * 1.1)

df_long <- df_long %>% left_join(y_positions, by = "gene")

# Boxplot faceteado por gen
p <- ggboxplot(
  df_long,
  x = "group",
  y = "expression",
  color = "group",
  palette = c("blue", "brown1"),
  add = "jitter"
) +
  stat_compare_means(
    comparisons = my_comparisons,
    label = "p.signif",
    method = "wilcox.test",
    exact = FALSE,
    aes(y = y.position)
  ) +
  facet_wrap(~ gene, scales = "free_y", ncol = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(size = 10)
  )

# Guardar PDF
pdf(paste0(outdir, "/Tcell.Analysis/Fas.Score.genes.pdf"),
    width = 12, height = 12)
print(p)
dev.off()




#### Pseudobulk cluster basis


# counts (IMPORTANTE: no usar nombre "counts")
counts_mat <- GetAssayData(data, layer = "counts")

meta <- data@meta.data

# checks rápidos
table(meta$group)   # KO vs WT
table(meta$tag)     # muestras (deberían ser 8)
table(meta$Subcluster)


clusters <- unique(meta$Subcluster)

results_list <- list()
dds_list <- list()

for (cl in clusters) {
  
  message("Processing: ", cl)
  
  # -------------------------
  # Subset cluster
  # -------------------------
  cells <- meta$Subcluster == cl
  
  sub_counts <- counts_mat[, cells]
  sub_meta <- meta[cells, ]
  
  # sample = tag (CLAVE)
  sub_meta$sample_id <- sub_meta$tag
  
  # -------------------------
  # Filtrar muestras con pocas células
  # -------------------------
  cell_counts <- table(sub_meta$sample_id)
  keep_samples <- names(cell_counts[cell_counts >= 20])  # ajusta si quieres
  
  sub_meta <- sub_meta[sub_meta$sample_id %in% keep_samples, ]
  sub_counts <- sub_counts[, colnames(sub_counts) %in% rownames(sub_meta)]
  
  # si después del filtro no queda nada → skip
  if (ncol(sub_counts) == 0) {
    message("Skipping ", cl, " (no cells after filtering)")
    next
  }
  
  # -------------------------
  # Pseudobulk (sumar counts por muestra)
  # -------------------------
  pb_counts <- rowsum(
    t(as.matrix(sub_counts)),
    group = sub_meta$sample_id
  )
  
  pb_counts <- t(pb_counts)
  
  # -------------------------
  # Metadata pseudobulk
  # -------------------------
  pb_meta <- sub_meta %>%
    group_by(sample_id) %>%
    summarise(group = unique(group))
  
  pb_meta <- as.data.frame(pb_meta)
  rownames(pb_meta) <- pb_meta$sample_id
  
  # asegurar mismo orden
  pb_counts <- pb_counts[, rownames(pb_meta)]
  
  # -------------------------
  # Checks de calidad
  # -------------------------
  if (length(unique(pb_meta$group)) < 2) {
    message("Skipping ", cl, " (no KO/WT balance)")
    next
  }
  
  if (any(table(pb_meta$group) < 2)) {
    message("Skipping ", cl, " (not enough replicates)")
    next
  }
  
  # -------------------------
  # DESeq2
  # -------------------------
  dds <- DESeqDataSetFromMatrix(
    countData = round(pb_counts),
    colData = pb_meta,
    design = ~ group
  )
  
  # filtrar genes poco expresados
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = c("group", "KO", "WT"))
  
  # guardar
  results_list[[cl]] <- res
  dds_list[[cl]] <- dds
}





top20_list <- list()

for (cl in names(results_list)) {
  
  res <- as.data.frame(results_list[[cl]])
  
  # quitar NA (muy importante)
  res <- res[!is.na(res$log2FoldChange), ]
  
  # ordenar por valor absoluto de logFC
  res <- res[order(abs(res$log2FoldChange), decreasing = TRUE), ]
  
  # coger top 20
  top20 <- head(res, 20)
  
  # añadir nombre del gen
  top20$gene <- rownames(top20)
  
  # guardar
  top20_list[[cl]] <- top20
}

library(dplyr)

top20_df <- bind_rows(
  lapply(names(top20_list), function(cl) {
    
    df <- top20_list[[cl]]
    df$cluster <- cl
    return(df)
    
  })
)


## Quantification




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
pdf(paste0(outdir,"/Tcell.Analysis/quantification.tags.pdf"), width=24, height=12)
ggplot(prop_df, aes(x = tag, y = prop, fill = Subcluster)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(x = "Tag", y = "Cell proportion (Downsampled to 8k)", fill = "Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

