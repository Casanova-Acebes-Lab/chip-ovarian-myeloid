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
library(variancePartition)
library(SummarizedExperiment)
library(BiocParallel)
library(Matrix.utils)



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

data <- readRDS(paste0(rsdir,"objects/data.macrophages.Clusterized2.rds"))




### Pseudobulk


data@meta.data$sample <- data@meta.data$tag
table(data@meta.data$sample)

#####


# Dimensionality reduction plots
library(sva)
mat_corrected <- ComBat(dat = as.matrix(Ly6cHi_Monocytes[[3]]), 
                        batch = Ly6cHi_Monocytes[[2]]$batch, 
                        mod = model.matrix(~ genotype, data = Ly6cHi_Monocytes[[2]]))






# expr_corrected: genes x samples (como en tu función residualize_for_covariates)
mat <- Ly6cHi_Monocytes[[7]]




# 1️⃣ Calcular distancia euclidiana entre columnas (muestras)
dist_mat <- dist(t(mat))  # transponemos para que dist calcule entre muestras

# 2️⃣ MDS
mds <- cmdscale(dist_mat, k = 2)

# 3️⃣ Convertir a data.frame para ggplot
mds_df <- data.frame(
  MDS1 = mds[,1],
  MDS2 = mds[,2],
  sample = colnames(mat)
)

# 4️⃣ Revisar
head(mds_df)


mds_df$genotype <- Ly6cHi_Monocytes[[2]]$genotype[match(mds_df$sample, Ly6cHi_Monocytes[[2]]$sample)]


pdf(paste0(outdir,"/Pseudobulk2/MDS.Ly6cHi_Monocytes2.pdf"), width=16, height=12)
ggplot(mds_df, aes(x = MDS1, y = MDS2, color = sample)) +
  geom_point(size = 4) +
  geom_text(aes(label = sample), vjust = -1.2) +
  theme_minimal() +
  ggtitle("MDS de pseudobulk residualizado")
dev.off()





library_entropy <- apply(Ly6cHi_Monocytes[[8]], 2, function(x) {
  p <- x / sum(x)
  -sum(p * log(p + 1e-12))
})



# --------------------------------------------
# Lets go
# --------------------------------------------


 Ly6cHi_Monocytes  <- downsample_pseudobulk_iterations(data, "Ly6cHi Monocytes")
# Ordenar por p-valor ajustado
tt_ordered <- Ly6cHi_Monocytes[[6]][order(Ly6cHi_Monocytes[[6]]$adj.P.Val), ]

# Seleccionar las 50 primeras filas
top50 <- head(tt_ordered, 100)

# Mostrar
top50


pdf(paste0(outdir,"/Pseudobulk2/Volcano.Ly6cHi_Monocytes.pdf"), width=16, height=12)
p <- Volcano2(Ly6cHi_Monocytes[[6]], "Monocytes")
p
dev.off()



# Ly6cHi_Monocytes
Ly6cHi_Monocytes <- run_pseudobulk_analysis(data, "Ly6cHi Monocytes", outdir, rsdir)

Ly6cHi_Monocytes[[3]]
head(Ly6cHi_Monocytes[[2]], n=200)
Ly6cHi_Monocytes[[4]]
head(Ly6cHi_Monocytes[[8]], n=50)
head(Ly6cHi_Monocytes[[9]], n=50)





# Ly6cLo_Monocytes
Ly6cLo_Monocytes <- run_pseudobulk_analysis(data, "Ly6cLo Monocytes", outdir, rsdir)

Ly6cLo_Monocytes[[1]]
head(Ly6cLo_Monocytes[[2]], n=50)
Ly6cLo_Monocytes[[3]]
Ly6cLo_Monocytes[[4]]
Ly6cLo_Monocytes[[9]]




# Early IFN|MHCII-TAMs

Early_IFN_MHCII_TAMs <- pseudobulk_limma_SVA(data, "Early IFN|MHCII-TAMs")  
Early_IFN_MHCII_TAMs[[1]]       
head(Early_IFN_MHCII_TAMs[[2]], n=50)
Early_IFN_MHCII_TAMs[[3]]
Early_IFN_MHCII_TAMs[[4]]
Early_IFN_MHCII_TAMs[[9]]



tt_ordered <- Early_IFN_MHCII_TAMs[[9]][order(Early_IFN_MHCII_TAMs[[9]]$adj.P.Val), ]

# Seleccionar las 50 primeras filas
top50 <- head(tt_ordered, 100)

# Mostrar
top50


# Arg1|Spp1|Mmp12|Mmp19|Il1a Mac
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac <- run_pseudobulk_analysis(data, "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac", outdir, rsdir)
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[1]]
head(Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[2]], n=50)
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[3]]
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[4]]
Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[9]]

# Trem1|Ptgs2|Plaur|Celc4e Mac
Trem1_Ptgs2_Plaur_Celc4e_Mac <- run_pseudobulk_analysis(data, "Trem1|Ptgs2|Plaur|Celc4e Mac", outdir, rsdir)    
Trem1_Ptgs2_Plaur_Celc4e_Mac[[1]]
head(Trem1_Ptgs2_Plaur_Celc4e_Mac[[2]], n=50)
Trem1_Ptgs2_Plaur_Celc4e_Mac[[3]]
Trem1_Ptgs2_Plaur_Celc4e_Mac[[4]]
Trem1_Ptgs2_Plaur_Celc4e_Mac[[9]]

# MHCII|Ccl12 Mac
MHCII_Ccl12_Mac <- (data, "MHCII|Ccl12 Mac", outdir, rsdir)
MHCII_Ccl12_Mac[[1]]
head(MHCII_Ccl12_Mac[[2]], n=50)
MHCII_Ccl12_Mac[[3]]
MHCII_Ccl12_Mac[[4]]
MHCII_Ccl12_Mac[[9]]


# MHCII|Siglec Mac
MHCII_Siglec_Mac <- pseudobulk_limma(data, "MHCII|Siglec Mac")
MHCII_Siglec_Mac[[1]]
head(MHCII_Siglec_Mac[[9]], n=50)
MHCII_Siglec_Mac[[3]]
MHCII_Siglec_Mac[[4]]
MHCII_Siglec_Mac[[9]]


tt_ordered <- MHCII_Siglec_Mac[[6]][order(MHCII_Siglec_Mac[[6]]$adj.P.Val), ]

# Seleccionar las 50 primeras filas
top50 <- head(tt_ordered, 100)

# Mostrar
top50




pdf(paste0(outdir,"/Pseudobulk2/MHCII|Siglec Mac.pdf"), width=16, height=12)
p <- Volcano2(MHCII_Siglec_Mac[[6]], "MHCII|Siglec Mac")
p
dev.off()



# IFN Mac
IFN_Mac <- pseudobulk_limma_avg_per_cell_with_umis(data, "IFN Mac")
IFN_Mac[[1]]
head(IFN_Mac[[9]], n=50)
IFN_Mac[[3]]
IFN_Mac[[4]]
IFN_Mac[[9]]



tt_ordered <- IFN_Mac[[6]][order(IFN_Mac[[6]]$adj.P.Val), ]

# Seleccionar las 50 primeras filas
top50 <- head(tt_ordered, 100)

# Mostrar
top50


# Mmp9|Ctsk Mac
Mmp9_Ctsk_Mac <- run_pseudobulk_analysis(data, "Mmp9|Ctsk Mac", outdir, rsdir)
Mmp9_Ctsk_Mac[[1]]
head(Mmp9_Ctsk_Mac[[2]], n=50)
Mmp9_Ctsk_Mac[[3]]
Mmp9_Ctsk_Mac[[4]]
Mmp9_Ctsk_Mac[[9]]


# Mrc1|C1qc|Cbr2|Gas6 Mac
Mrc1_C1qc_Cbr2_Gas6_Mac <- pseudobulk_limma_avg_per_cell_with_umis(data, "Mrc1|C1qc|Cbr2|Gas6 Mac")
Mrc1_C1qc_Cbr2_Gas6_Mac[[1]]
head(Mrc1_C1qc_Cbr2_Gas6_Mac[[2]], n=50)
Mrc1_C1qc_Cbr2_Gas6_Mac[[3]]
Mrc1_C1qc_Cbr2_Gas6_Mac[[4]]
Mrc1_C1qc_Cbr2_Gas6_Mac[[9]]

tt_ordered <- Mrc1_C1qc_Cbr2_Gas6_Mac[[6]][order(Mrc1_C1qc_Cbr2_Gas6_Mac[[6]]$adj.P.Val), ]

# Seleccionar las 50 primeras filas
top50 <- head(tt_ordered, 100)

# Mostrar
top50



# Npr2|Actn1 Mac
Npr2_Actn1_Mac <- run_pseudobulk_analysis(data, "Npr2|Actn1 Mac", outdir, rsdir)
Npr2_Actn1_Mac[[1]]
head(Npr2_Actn1_Mac[[2]], n=50)
Npr2_Actn1_Mac[[3]]
Npr2_Actn1_Mac[[4]]
Npr2_Actn1_Mac[[9]]

# Fn1|Vegfa Mac
Fn1_Vegfa_Mac <- run_pseudobulk_analysis(data, "Fn1|Vegfa Mac", outdir, rsdir)
Fn1_Vegfa_Mac[[1]]
head(Fn1_Vegfa_Mac[[2]], n=50)
Fn1_Vegfa_Mac[[3]]
Fn1_Vegfa_Mac[[4]]
Fn1_Vegfa_Mac[[9]]

# Neutrophils
Neutrophils <- run_pseudobulk_analysis(data, "Neutrophils", outdir, rsdir)
Neutrophils[[1]]
head(Neutrophils[[2]], n=50)
Neutrophils[[3]]
Neutrophils[[4]]
Neutrophils[[9]]



# Save tables
write.table(Ly6cHi_Monocytes[[2]], paste0(rsdir,"table.macros.Ly6cHi Monocytes.tsv"), sep='\t')
write.table(Ly6cLo_Monocytes[[2]], paste0(rsdir,"table.macros.Ly6cLo Monocytes.tsv"), sep='\t')
write.table(Early_IFN_MHCII_TAMs[[2]], paste0(rsdir,"table.macros.Early IFN MHCII TAMs.tsv"), sep='\t')
write.table(Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[2]], paste0(rsdir,"table.macros.Arg1 Spp1 Mmp12 Mmp19 Il1a Mac.tsv"), sep='\t')
write.table(Trem1_Ptgs2_Plaur_Celc4e_Mac[[2]], paste0(rsdir,"table.macros.Trem1 Ptgs2 Plaur Celc4e Mac.tsv"), sep='\t')
write.table(MHCII_Ccl12_Mac[[2]], paste0(rsdir,"table.macros.MHCII Ccl12 Mac.tsv"), sep='\t')
write.table(MHCII_Siglec_Mac[[2]], paste0(rsdir,"table.macros.MHCII Siglec Mac.tsv"), sep='\t')
write.table(IFN_Mac[[2]], paste0(rsdir,"table.macros.IFN Mac.tsv"), sep='\t')
write.table(Mmp9_Ctsk_Mac[[2]], paste0(rsdir,"table.macros.Mmp9 Ctsk Mac.tsv"), sep='\t')
write.table(Mrc1_C1qc_Cbr2_Gas6_Mac[[2]], paste0(rsdir,"table.macros.Mrc1 C1qc Cbr2 Gas6 Mac.tsv"), sep='\t')
write.table(Npr2_Actn1_Mac[[2]], paste0(rsdir,"table.macros.Npr2 Actn1 Mac.tsv"), sep='\t')
write.table(Fn1_Vegfa_Mac[[2]], paste0(rsdir,"table.macros.Fn1 Vegfa Mac.tsv"), sep='\t')
write.table(Neutrophils[[2]], paste0(rsdir,"table.macros.Neutrophils.tsv"), sep='\t')      







write.table(Ly6cHi_Monocytes[[9]], paste0(rsdir,"matrix.macros.Ly6cHi Monocytes.tsv"), sep='\t')
write.table(Ly6cLo_Monocytes[[9]], paste0(rsdir,"matrix.macros.Ly6cLo Monocytes.tsv"), sep='\t')
write.table(Early_IFN_MHCII_TAMs[[9]], paste0(rsdir,"matrix.macros.Early IFN MHCII TAMs.tsv"), sep='\t')
write.table(Arg1_Spp1_Mmp12_Mmp19_Il1a_Mac[[9]], paste0(rsdir,"matrix.macros.Arg1 Spp1 Mmp12 Mmp19 Il1a Mac.tsv"), sep='\t')
write.table(Trem1_Ptgs2_Plaur_Celc4e_Mac[[9]], paste0(rsdir,"matrix.macros.Trem1 Ptgs2 Plaur Celc4e Mac.tsv"), sep='\t')
write.table(MHCII_Ccl12_Mac[[9]], paste0(rsdir,"matrix.macros.MHCII Ccl12 Mac.tsv"), sep='\t')
write.table(MHCII_Siglec_Mac[[9]], paste0(rsdir,"matrix.macros.MHCII Siglec Mac.tsv"), sep='\t')
write.table(IFN_Mac[[9]], paste0(rsdir,"matrix.macros.IFN Mac.tsv"), sep='\t')
write.table(Mmp9_Ctsk_Mac[[9]], paste0(rsdir,"matrix.macros.Mmp9 Ctsk Mac.tsv"), sep='\t')
write.table(Mrc1_C1qc_Cbr2_Gas6_Mac[[9]], paste0(rsdir,"matrix.macros.Mrc1 C1qc Cbr2 Gas6 Mac.tsv"), sep='\t')
write.table(Npr2_Actn1_Mac[[9]], paste0(rsdir,"matrix.macros.Npr2 Actn1 Mac.tsv"), sep='\t')
write.table(Fn1_Vegfa_Mac[[9]], paste0(rsdir,"matrix.macros.Fn1 Vegfa Mac.tsv"), sep='\t')
write.table(Neutrophils[[9]], paste0(rsdir,"matrix.macros.Neutrophils.tsv"), sep='\t')      



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




# --- Definir epsilon por si acaso hay NA ---
epsilon <- 1e-8


all_clusters <- all_clusters %>%
  mutate(cluster = case_when(
    cluster == "Early_IFN_MHCII_TAMs" ~ "Early IFN|MHCII-TAMs",
    TRUE ~ cluster
  ))


# --- Preparar los datos ---
markers_macro <- all_clusters %>%
  filter(cluster %in% macro_clusters) %>%        # filtrar clusters de interés
  mutate(
    # NUEVA MÉTRICA: logFC ponderado por especificidad (% diferencia)
    metric = logFC * -log10(adj.P.Val + epsilon),
    cluster = factor(cluster, levels = macro_clusters)
    )  # mantener orden
  

# ============================
# 2. SELECCIONAR TOP 10 ↑ Y 10 ↓ POR CLUSTER
# ============================



# --- Filtrar solo genes diferencialmente expresados ---
markers_macro_filtered <- markers_macro %>%
  filter(diffexpressed %in% c("Up","Down")) %>%
  group_by(cluster, gene) %>%
  slice_max(order_by = abs(metric), n = 1) %>%   # nos quedamos con la fila "más significativa"
  ungroup()


# --- Selección de top genes por cluster ---
# --- Selección de top genes por cluster (10 Up + 10 Down sin duplicados) ---
top_genes2 <- markers_macro_filtered %>%
  group_by(cluster) %>%
  # top 10 Up
  filter(diffexpressed == "Up") %>%
  slice_max(metric, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  bind_rows(
    markers_macro_filtered %>%
      group_by(cluster) %>%
      # top 10 Down
      filter(diffexpressed == "Down") %>%
      slice_min(metric, n = 10, with_ties = FALSE) %>%
      ungroup()
  ) %>%
  mutate(
    color = ifelse(diffexpressed == "Up", "red", "blue")
  )
  
  top_genes2 <- as.data.frame(top_genes2)

# ============================
# 3. CAJAS DE COLOR BAJO CADA CLUSTER
# ============================

clusters2 <- macro_clusters


cluster_boxes2 <- data.frame(
  cluster = clusters2,
  xmin = seq_along(clusters2) - 0.4,
  xmax = seq_along(clusters2) + 0.4,
  ymin = -0.15,
  ymax = 0.15
)

cluster_boxes2$cluster <- factor(cluster_boxes2$cluster, levels = clusters2)
cluster_boxes2$fill <- nora.colors2[as.character(cluster_boxes2$cluster)]

# ============================
# 4. PLOT FINAL
# ============================



pdf(paste0(outdir, "/Pseudobulk/KOvsWT.DEGs.All_clusters.pdf"), width = 26, height = 14)

ggplot(markers_macro, aes(x = cluster, y = metric)) +

  # cajas de color
  geom_rect(data = cluster_boxes2,
            aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = cluster),
            inherit.aes = FALSE, alpha = 0.8) +

  # todos los puntos
  geom_jitter(width = 0.25, alpha = 0.4, color = "grey") +

  # top genes
  geom_point(
    data = top_genes2,
    aes(x = cluster, y = metric, color = diffexpressed),
    size = 2,
    position = position_jitter(width = 0.25)
  ) +

  # labels de top genes
  geom_text_repel(
    data = top_genes2,
    aes(x = cluster, y = metric, label = gene),
    size = 4,
    max.overlaps = 60,
    segment.size = 0.3,
    position = position_jitter(width = 0.25)
  ) +

  # definir colores para Up y Down
  scale_color_manual(name = "Regulation", values = c("Up" = "red", "Down" = "blue")) +
  scale_fill_manual(values = nora.colors2) +

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
    title = "KO vs WT Chimeras. Top differentially expressed genes per cluster (limma DE)"
  )

dev.off()






#### Seurat MAST


library(Seurat)
library(dplyr)

# ---- PARÁMETROS ----
cluster_col <- "Clustering.Round2"
cluster_of_interest <- "Ly6cHi Monocytes"

genotype_col <- "group"   # WT / KO
latent_var <- "nCount_RNA"

# ---- 1. Seleccionar células del cluster ----
cells_use <- rownames(data@meta.data)[ data@meta.data[[cluster_col]] == cluster_of_interest ]

if (length(cells_use) == 0) stop("No hay células del cluster seleccionado")

message("Células seleccionadas: ", length(cells_use))

# ---- 2. Crear un objeto SOLO con el cluster ----
sub <- subset(data, cells = cells_use)

# ---- 3. Asignar las identidades al GENOTIPO ----
Idents(sub) <- sub[[genotype_col]][,1]

levels(Idents(sub))
# Debe mostrar: "WT" "KO"

# ---- 4. Ejecutar FindMarkers dentro del cluster ----
markers <- FindMarkers(
    sub,
    ident.1 = "KO",
    ident.2 = "WT",
    test.use = "MAST",
    latent.vars = latent_var
)

# ---- 5. Ordenar y añadir nombres ----
markers$gene <- rownames(markers)
markers <- markers %>% arrange(p_val_adj)

# ---- 6. Ver resultados ----
head(markers, 50)
library(dplyr)

# Proporción de células expresando cada gen
expr_binary <- expr > 0
pct_expr <- data.frame(
  gene = rownames(expr),
  pct_WT = rowSums(expr_binary[, cond == "WT"]) / sum(cond == "WT"),
  pct_KO = rowSums(expr_binary[, cond == "KO"]) / sum(cond == "KO")
)

# Filtrar genes “problemáticos” suavizado
df_problematic <- pct_expr %>%
  filter(
    (pct_WT < 0.1 & pct_KO > 0.3) |   # <10% vs >30%
    (pct_KO < 0.1 & pct_WT > 0.3)
  ) %>%
  arrange(desc(abs(pct_WT - pct_KO)))

# Ver top 20
head(df_problematic, 20)




library(Seurat)
library(dplyr)

# ---- PARÁMETROS ----
cluster_col <- "Clustering.Round2"
cluster_of_interest <- "Trem1|Ptgs2|Plaur|Celc4e Mac"
genotype_col <- "group"   # WT / KO
latent_var <- c("batch")

# ---- 1. Seleccionar células del cluster ----
cells_use <- rownames(data@meta.data)[ data@meta.data[[cluster_col]] == cluster_of_interest ]
if (length(cells_use) == 0) stop("No hay células del cluster seleccionado")
message("Número de células seleccionadas: ", length(cells_use))

# ---- 2. Crear un objeto SOLO con el cluster ----
sub <- subset(data, cells = cells_use)

# ---- 3. Asignar las identidades al GENOTIPO ----
Idents(sub) <- factor(sub[[genotype_col]][,1])
levels(Idents(sub))  # Debe mostrar: "WT" "KO"

# ---- 4. Downsamplear al grupo más pequeño ----
set.seed(123)  # reproducibilidad
group_counts <- table(Idents(sub))
min_cells <- min(group_counts)

cells_downsampled <- unlist(lapply(names(group_counts), function(g) {
  sample(which(Idents(sub) == g), min_cells)
}))

sub_down <- subset(sub, cells = cells_downsampled)
message("Células después del downsampling: ")
table(Idents(sub_down))

# ---- 5. Ejecutar FindMarkers dentro del cluster con MAST ----

DefaultAssay(data) <- "SCT"
markers <- FindMarkers(
  data,
  ident.1 = "KO",
  ident.2 = "WT",
  test.use = "MAST",
  min.pct = 0.2
)

# ---- 6. Ver resultados ----
head(markers, 100)



DefaultAssay(data) <- "RNA"
colnames(data@meta.data)




pdf(paste0(outdir, "/Pseudobulk2/violin.Il1b.pdf"), width = 26, height = 14)
VlnPlot(
  data,
  features = "Il1b",
  group.by = "Clustering.Round2",
  pt.size = 0,
  split.by="chimera"
)
dev.off()