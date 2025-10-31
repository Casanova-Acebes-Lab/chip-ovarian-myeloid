

# Enviroment R.4.3.3
library(Seurat)
library(dplyr)
library(cowplot)
library(ggplot2)
library(enrichplot)
library(org.Mm.eg.db)
library(clusterProfiler)
library(biomaRt)





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


########### Read data


cluster.Arg1 <- read.table(paste0(rsdir,"table.macros.Arg1+ Mac.tsv"), sep='\t', header=T)
cluster.Mrc1 <- read.table(paste0(rsdir,"table.macros.Mrc1+ Mac.tsv"), sep='\t', header=T)
cluster.Mki67 <- read.table(paste0(rsdir,"table.macros.Mki67+ Mac.tsv"), sep='\t', header=T)
cluster.Mmp9 <- read.table(paste0(rsdir,"table.macros.Mmp9+ Mac.tsv"), sep='\t', header=T)
cluster.Gpnmb <- read.table(paste0(rsdir,"table.macros.Gpnmb+ Mac.tsv"), sep='\t', header=T)
cluster.Arg1.Mrc1 <- read.table(paste0(rsdir,"table.macros.Arg1|Mrc1 Mac.tsv"), sep='\t', header=T)
cluster.Ciita.Mrc1.MHCII <- read.table(paste0(rsdir,"table.macros.Ciita|Mrc1|MHC-II Mac.tsv"), sep='\t', header=T)
cluster.Monocytes.Ly6c2Hi <- read.table(paste0(rsdir,"table.macros.Monocytes Ly6c2Hi.tsv"), sep='\t', header=T)
cluster.IFN <- read.table(paste0(rsdir,"table.macros.IFN Mac.tsv"), sep='\t', header=T)
cluster.Neutrophils <- read.table(paste0(rsdir,"table.macros.Neutrophils.tsv"), sep='\t', header=T)
cluster.ALL <- read.table(paste0(rsdir,"table.macros.ALL Mac.tsv"), sep='\t', header=T)





# GSEA



GSEA.Arg1 <- GSEA(cluster.Arg1)
head(GSEA.Arg1,n=20)


GSEA.Mrc1 <- GSEA(cluster.Mrc1)
head(GSEA.Mrc1,n=20)


GSEA.Mki67 <- GSEA(cluster.Mki67)
head(GSEA.Mki67,n=20)


GSEA.Mmp9 <- GSEA(cluster.Mmp9)
head(GSEA.Mmp9,n=20)


GSEA.Gpnmb <- GSEA(cluster.Gpnmb)
head(GSEA.Gpnmb,n=20)



GSEA.Arg1.Mrc1 <- GSEA(cluster.Arg1.Mrc1)
head(GSEA.Arg1.Mrc1,n=20)


GSEA.Ciita.Mrc1.MHCII <- GSEA(cluster.Ciita.Mrc1.MHCII)
head(GSEA.Ciita.Mrc1.MHCII,n=20)


GSEA.Monocytes.Ly6c2Hi <- GSEA(cluster.Monocytes.Ly6c2Hi)
head(GSEA.Monocytes.Ly6c2Hi,n=20)



GSEA.IFN <- GSEA(cluster.IFN)
head(GSEA.IFN,n=20)



GSEA.Neutrophils <- GSEA(cluster.Neutrophils)
head(GSEA.Neutrophils,n=20)




GSEA.Macros <- GSEA(cluster.ALL)
head(GSEA.Macros,n=20)



#Ploting Things


p <- plot_gsea(GSEA.Arg1@result, title = "Macrophages Arg1+. KO vs WT")


pdf(paste0(outdir,"/GSEA/Arg1.pdf"), width = 10, height = 12)
print(p)
dev.off()


p <- plot_gsea(GSEA.Mrc1@result, title = "Macrophages Mrc1+. KO vs WT")

pdf(paste0(outdir,"/GSEA/Mrc1.pdf"), width = 10, height = 12)
print(p)
dev.off()


p <- plot_gsea(GSEA.Mki67@result, title = "Macrophages Mki67+. KO vs WT")

pdf(paste0(outdir,"/GSEA/Mki67.pdf"), width = 17, height = 18)
print(p)
dev.off()



p <- plot_gsea(GSEA.Mmp9@result, title = "Macrophages Mmp9+. KO vs WT")

pdf(paste0(outdir,"/GSEA/Mmp9.pdf"), width = 10, height = 16)
print(p)
dev.off()



p <- plot_gsea(GSEA.Gpnmb@result, title = "Macrophages Gpnmb+. KO vs WT")

pdf(paste0(outdir,"/GSEA/Gpnmb.pdf"), width = 16, height = 12)
print(p)
dev.off()


p <- plot_gsea(GSEA.Arg1.Mrc1@result, title = "Macrophages Arg1.Mrc1+. KO vs WT")

pdf(paste0(outdir,"/GSEA/Arg1.Mrc1.pdf"), width = 16, height = 18)
print(p)
dev.off()



p <- plot_gsea(GSEA.Ciita.Mrc1.MHCII@result, title = "Macrophages Ciita.Mrc1.MHCII+. KO vs WT")

pdf(paste0(outdir,"/GSEA/Ciita.Mrc1.MHCII.pdf"), width = 10, height = 18)
print(p)
dev.off()


p <- plot_gsea(GSEA.Monocytes.Ly6c2Hi@result, title = "Macrophages Monocytes.Ly6c2Hi+. KO vs WT")

pdf(paste0(outdir,"/GSEA/Monocytes.Ly6c2Hi.pdf"), width = 10, height = 16)
print(p)
dev.off()


p <- plot_gsea(GSEA.IFN@result, title = "Macrophages IFN. KO vs WT")

pdf(paste0(outdir,"/GSEA/IFN.pdf"), width = 10, height = 16)
print(p)
dev.off()



p <- plot_gsea(GSEA.Neutrophils@result, title = "Neutrophils. KO vs WT")

pdf(paste0(outdir,"/GSEA/Neutrophils.pdf"), width = 10, height = 16)
print(p)
dev.off()



p <- plot_gsea(GSEA.Macros@result, title = "Macrophages All clusters. KO vs WT")

pdf(paste0(outdir,"/GSEA/All.pdf"), width = 10, height = 16)
print(p)
dev.off()






### Enrichment for functional analysis of clusters


clusters<- read.table(paste0(rsdir,"table.clusters.markers.filt.Clusterized.R2.tsv"), sep='\t', header=T)


# Función para extraer filas de un cluster específico
extraer_cluster <- function(df, cluster_name) {
  subset(df, cluster == cluster_name)
}

# Ejemplo de uso:
# Extraer todas las filas del cluster "Ly6c|Ms4ac Monocytes"
cluster_ly6c <- extraer_cluster(clusters, "Ly6c|Ms4ac Monocytes")

# Ver las primeras filas
dim(cluster_ly6c)


# 1️⃣ Filtrar solo los genes UP y diferencialmente expresados
up_genes <- cluster_ly6c %>%
  filter(avg_log2FC > 0, p_val_adj < 0.05) %>%
  pull(gene)

# 2️⃣ Mapear los símbolos de gen a ENTREZID (requerido por enrichGO)
gene_entrez <- bitr(up_genes, fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = org.Mm.eg.db)

# 3️⃣ Hacer enrichGO
ego <- enrichGO(gene          = gene_entrez$ENTREZID,
                OrgDb         = org.Mm.eg.db,
                keyType       = "ENTREZID",
                ont           = "BP",        # Biological Process
                pAdjustMethod = "BH",
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.2,
                readable      = TRUE)       # convierte de vuelta a SYMBOL





# seleccionar los 100 términos más significativos, pero mantener el objeto enrichResult
ego_top <- ego
ego_top@result <- ego@result[order(ego@result$p.adjust), ][1:100, ]

# ahora sí podemos simplificar
ego_simplified <- simplify(ego_top, cutoff = 0.7, by = "p.adjust", select_fun = min)


# 4️⃣ Revisar resultados
head(ego_simplified @result, n=10)



# Función completa para análisis GO de un cluster específico
enrichGO_cluster <- function(df, cluster_name,
                             adjp_cutoff = 0.05,
                             top_terms = 100,
                             simplify_cutoff = 0.7) {
  
  # 1️⃣ Extraer filas del cluster
  cluster_df <- df %>% filter(cluster == cluster_name)
  
  # 2️⃣ Filtrar genes UP y diferencialmente expresados
  up_genes <- cluster_df %>%
    filter(avg_log2FC > 0, p_val_adj < adjp_cutoff) %>%
    pull(gene)
  
  if(length(up_genes) == 0){
    warning(paste("No hay genes UP significativos para el cluster", cluster_name))
    return(NULL)
  }
  
  # 3️⃣ Mapear símbolos de gen a ENTREZID
  gene_entrez <- bitr(up_genes,
                      fromType = "SYMBOL",
                      toType   = "ENTREZID",
                      OrgDb    = org.Mm.eg.db)
  
  if(nrow(gene_entrez) == 0){
    warning(paste("No se pudieron mapear genes a ENTREZID para", cluster_name))
    return(NULL)
  }
  
  # 4️⃣ EnrichGO
  ego <- enrichGO(gene          = gene_entrez$ENTREZID,
                  OrgDb         = org.Mm.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff  = adjp_cutoff,
                  qvalueCutoff  = 0.2,
                  readable      = TRUE)
  
  # 5️⃣ Limitar a los top_terms más significativos
  ego_top <- ego
  ego_top@result <- ego@result %>%
    arrange(p.adjust) %>%
    head(top_terms)
  
  # 6️⃣ Simplificar GO para eliminar redundancias
  ego_simplified <- simplify(ego_top,
                             cutoff = simplify_cutoff,
                             by = "p.adjust",
                             select_fun = min)
  
  return(ego_simplified)
}


ego_ly6c <- enrichGO_cluster(clusters, "Ly6c|Ms4ac Monocytes")
ego_ly6c_sorted <- ego_ly6c
ego_ly6c_sorted@result <- ego_ly6c@result[order(ego_ly6c@result$GeneRatio, decreasing = TRUE), ]
head(ego_ly6c_sorted@result, 20)


ego_Arg1.Spp1 <- enrichGO_cluster(clusters, "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac")
head(ego_Arg1.Spp1@result, 20)

ego_Ciita.Ccl12<- enrichGO_cluster(clusters, "Ciita|Ccl12 Mac")
head(ego_Ciita.Ccl12@result, 20)


ego_Ciita.Siglec <- enrichGO_cluster(clusters, "Ciita|Siglec Mac")
head(ego_Ciita.Siglec@result, 20)


ego_Fn1.Vegfa <- enrichGO_cluster(clusters, "Fn1|Vegfa Mac")
head(ego_Fn1.Vegfa@result, 20)


ego_IFN <- enrichGO_cluster(clusters, "IFN Mac")
head(ego_IFN@result, 20)

ego_Mmp9.Ctsk  <- enrichGO_cluster(clusters, "Mmp9|Ctsk Mac")
head(ego_Mmp9.Ctsk@result, 20)

ego_Mrc1.C1qc <- enrichGO_cluster(clusters, "Mrc1|C1qc|Cbr2|Gas6 Mac")
head(ego_Mrc1.C1qc@result, 20)


ego_Neutrophils <- enrichGO_cluster(clusters, "Neutrophils")
head(ego_Neutrophils@result, 20)


ego_Npr2.Actn1 <- enrichGO_cluster(clusters, "Npr2|Actn1 Mac")
head(ego_Npr2.Actn1@result, 20)

ego_Slamf <- enrichGO_cluster(clusters, "Slamf Monocytes")
head(ego_Slamf@result, 20)

ego_Trem1 <- enrichGO_cluster(clusters, "Trem1|Ptgs2|Plaur|Celc4e Mac")
head(ego_Trem1@result, 20)



library(dplyr)
library(clusterProfiler)

# Lista de clusters
cluster_list <- c(
  "Ly6c|Ms4ac Monocytes",
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac",
  "Ciita|Ccl12 Mac",
  "Ciita|Siglec Mac",
  "Fn1|Vegfa Mac",
  "IFN Mac",
  "Mmp9|Ctsk Mac",
  "Mrc1|C1qc|Cbr2|Gas6 Mac",
  "Neutrophils",
  "Npr2|Actn1 Mac",
  "Slamf Monocytes",
  "Trem1|Ptgs2|Plaur|Celc4e Mac"
)

# Crear un dataframe vacío
ego_df <- data.frame()

# Iterar sobre los clusters
for (cl in cluster_list) {
  
  # Enrichment
  ego <- enrichGO_cluster(clusters, cl)
  
  # Ordenar por p.adjust ascendente y seleccionar top 20
  ego_sorted <- ego
  ego_sorted@result <- ego@result %>%
    arrange(p.adjust) %>%
    slice_head(n = 20)
  
  # Añadir columna de cluster
  ego_cluster_df <- ego_sorted@result %>%
    mutate(cluster = cl)
  
  # Combinar con el dataframe principal
  ego_df <- bind_rows(ego_df, ego_cluster_df)
}

# Revisar el resultado
head(ego_df, 20)
dim(ego_df)  # número total de filas


# Definir el orden deseado
cluster_order <- c(
  "Ly6c|Ms4ac Monocytes",
  "Trem1|Ptgs2|Plaur|Celc4e Mac",
  "Mrc1|C1qc|Cbr2|Gas6 Mac",
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac",
  "Npr2|Actn1 Mac",
  "Slamf Monocytes",
  "Mmp9|Ctsk Mac",
  "IFN Mac",
  "Fn1|Vegfa Mac",
  "Ciita|Siglec Mac",
  "Ciita|Ccl12 Mac",
  "Neutrophils"
)

# Convertir la columna cluster en factor con el orden deseado
ego_df$cluster <- factor(ego_df$cluster, levels = cluster_order)

# Ordenar el data frame por cluster y por p.adjust (o GeneRatio si quieres)
ego_df <- ego_df %>%
  arrange(cluster, p.adjust)

# Revisar
head(ego_df, 20)




library(dplyr)

# Suponiendo que ego_df$geneID tenga los genes separados por "/"
ego_df_clean <- ego_df %>%
  mutate(gene_list = strsplit(as.character(geneID), "/")) %>%
  rowwise() %>%
  mutate(n_genes = length(gene_list)) %>%
  ungroup()

# Comparación manual por pares para filtrar redundancias
# Ejemplo simple: eliminar filas que tengan >70% genes en común con otra más significativa
threshold <- 0.9
keep <- rep(TRUE, nrow(ego_df_clean))

for(i in 1:(nrow(ego_df_clean)-1)){
  for(j in (i+1):nrow(ego_df_clean)){
    overlap <- length(intersect(ego_df_clean$gene_list[[i]], ego_df_clean$gene_list[[j]])) /
               min(length(ego_df_clean$gene_list[[i]]), length(ego_df_clean$gene_list[[j]]))
    if(overlap > threshold & ego_df_clean$p.adjust[j] > ego_df_clean$p.adjust[i]){
      keep[j] <- FALSE
    }
  }
}

ego_df_nonredundant <- ego_df_clean[keep, ]



library(dplyr)
library(ggplot2)
library(forcats)




nora.colors <- c(
  "Ly6c|Ms4ac Monocytes"         = "#FF3B30",   # rojo coral vivo
  "Trem1|Ptgs2|Plaur|Celc4e Mac" = "#EE7942",   # naranja vivo
  "Mrc1|C1qc|Cbr2|Gas6 Mac"      = "#FFD92F",   # amarillo brillante
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac" = "#4DAF4A", # verde intenso
  "Npr2|Actn1 Mac"               = "#A6D854",   # verde lima
  "Cd8 T cells"                  = "#00BFC4",   # turquesa
  "Cd4 T cells"                  = "#FF69B4",   # rosa fuerte
  "Neutrophils"                  = "#4876FF",   # azul fuerte
  "DCs"                           = "#87CEEB",   # azul claro
  "NK"                            = "#AB82FF",   # violeta oscuro
  "Tgd"                           = "#D9B3FF",   # lila pastel
  "Mastocytes"                    = "#3CB371",   # verde saturado
  "B cells"                       = "#DC143C",   # rojo intenso
  "Cd8 Effector"                  = "#A52A2A",   # rojo ladrillo
  "Marco+ Mac"                    = "#1E3A8A",   # azul intermedio
  "Slamf Monocytes"             = "#66B2FF",   # azul medio claro
  "Mmp9|Ctsk Mac"                  = "#00723F",   # verde botella
  "IFN Mac"                        = "#C080FF",   # morado claro
  "Fn1|Vegfa Mac"                  = "#FFA500",   # naranja estándar
  "Ciita|Siglec Mac"               = "#1E90FF",   # azul cobalto vivo
  "Ciita|Ccl12 Mac"                = "#4682B4",   # azul acero
  "Activated B cells"              = "#FF1493"    # fucsia intenso
)





library(dplyr)
library(ggplot2)
library(forcats)

# Función para preparar datos y plotear dotplot de enrichGO
plot_ego_dot_horizontal <- function(ego_df_nonredundant, top_n = 10, colors, cluster_order) {
  
  # Seleccionar top pathways por cluster según GeneRatio
  ego_top <- ego_df_nonredundant %>%
    group_by(cluster) %>%
    slice_max(order_by = GeneRatio, n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  # Convertir cluster a factor para respetar el orden deseado
  ego_top$cluster <- factor(ego_top$cluster, levels = cluster_order)
  
  # Para que los pathways estén ordenados verticalmente por GeneRatio dentro de cada cluster
  ego_top <- ego_top %>%
    group_by(cluster) %>%
    arrange(GeneRatio) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
    ungroup()
  
  # Plot
  p <- ggplot(ego_top, aes(x = cluster, y = Description, size = Count, color = cluster)) +
    geom_point(alpha = 0.8) +
    scale_color_manual(values = colors) +
    theme_bw(base_size = 14) +
    labs(
      x = "Cluster",
      y = "GO Term",
      size = "Number of Genes",
      title = paste0("Top ", top_n, " enriched GO terms per cluster from Overexpressed Genes")
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey90", size = 0.2)
    )
  
  return(p)
}

# --- Usar la función ---
p <- plot_ego_dot_horizontal(
  ego_df_nonredundant, 
  top_n = 10, 
  colors = nora.colors, 
  cluster_order = cluster_order
)

# Guardar a PDF
pdf(paste0(outdir,"/GSEA/Top10_GO_terms_by_cluster_horizontal.pdf"), width = 26, height = 18)
print(p)
dev.off()



library(ggplot2)
library(dplyr)
library(forcats)

# Convertir GeneRatio de "x/y" a numérico
ego_df_nonredundant <- ego_df_nonredundant %>%
  mutate(GeneRatio_numeric = sapply(strsplit(as.character(GeneRatio), "/"),
                                    function(x) as.numeric(x[1]) / as.numeric(x[2])))

plot_ego_heatmap <- function(ego_df_nonredundant, top_n = 10, cluster_order, value_col = "GeneRatio_numeric") {
  
  # Seleccionar top pathways por cluster
  ego_top <- ego_df_nonredundant %>%
    group_by(cluster) %>%
    slice_max(order_by = !!sym(value_col), n = top_n, with_ties = FALSE) %>%
    ungroup()
  
  # Ordenar clusters
  ego_top$cluster <- factor(ego_top$cluster, levels = cluster_order)
  
  # Ordenar pathways verticalmente por valor dentro de cada cluster
  ego_top <- ego_top %>%
    group_by(cluster) %>%
    arrange(!!sym(value_col)) %>%
    mutate(Description = factor(Description, levels = rev(unique(Description)))) %>%
    ungroup()
  
  # Plot heatmap con gradiente azul-rojo
  p <- ggplot(ego_top, aes(x = cluster, y = Description, fill = !!sym(value_col))) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "blue", high = "red") +
    theme_bw(base_size = 14) +
    labs(
      x = "Cluster",
      y = "GO Term",
      fill = value_col,
      title = paste0("Top ", top_n, " enriched GO terms per cluster")
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )
  
  return(p)
}

# Generar heatmap
p_heatmap <- plot_ego_heatmap(
  ego_df_nonredundant,
  top_n = 10,
  cluster_order = cluster_order,
  value_col = "GeneRatio_numeric"
)

# Guardar a PDF
pdf(paste0(outdir, "/GSEA/Top10_GO_terms_by_cluster_heatmap.pdf"), width = 26, height = 18)
print(p_heatmap)
dev.off()

