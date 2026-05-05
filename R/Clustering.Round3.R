
library(Seurat)
library(dplyr)
library(cowplot)
library(scCustomize)
library(harmony)
library(ggplot2)
library(RPresto)
library("scProportionTest")
library(ggrepel)





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



# Colors



nora.colors <- c(
  # --- MONOCITOS & TAMs (ROJOS) ---
  "Ly6cHi Monocytes"     = "#FF0000",  # Rojo inflamatorio puro
  "Ly6cLo Monocytes"     = "#FF6A6A",  # Rojo salmón transición
  "Early IFN|MHCII-TAMs"       = "#B22222",  # Rojo vino (pre-TAM IFN)
  "Trem1|Ptgs2|Plaur|Celc4e Mac" = "#EE7942",   # naranja vivo (inflam. activados)
  "Mrc1|C1qc|Cbr2|Gas6 Mac"      = "#FFD92F",   # amarillo brillante (TAM residentes)
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac" = "#4DAF4A", # verde TAM reparadores
  "Npr2|Actn1 Mac"       = "#A6D854",   # verde lima (TAM estructurales)
  "Mmp9|Ctsk Mac"        = "#00723F",   # verde botella (remodelado matriz)
  "IFN Mac"              = "#C080FF",   # morado claro (TAM interferón maduros)
  "Fn1|Vegfa Mac"        = "#FFA500",   # naranja angiogénico
  "MHCII|Siglec Mac"     = "#1E90FF",   # azul cobalto (antigen presenting)
  "MHCII|Ccl12 Mac"      = "#4682B4",   # azul acero

  # --- CÉLULAS LINFOIDES ---
  "Cd4 Naive"          = "#00BFC4",   # turquesa
  "Cd4 Activated"          = "#FF69B4",   # rosa fuerte
  "Cd8 Effector"         = "#698B69",   # rojo ladrillo
  "Cd4 Th17"                  = "#D9B3FF",   # lila pastel
  "Cd4 Treg"                  = "#36648B",   # lila pastel
  "NK"                   = "#AB82FF",   
  "Mki67+ Tcell"         = "#43CD80",   
  "Activated B cells"    = "#FF1493",   # fucsia intenso
  "B cells"              = "#DC143C",   # rojo intenso

  # --- INNATAS ---
  "Neutrophils"          = "#4876FF",   # azul fuerte
  "DCs"                  = "#87CEEB",   # azul claro
  "Mastocytes"           = "#3CB371"    # verde saturado
)



### Adding Tcells subclustering

t.annotation <- read.csv(paste0(outdir,"/Tcells/subclustering_tcells_annotations.csv"))

# Asegura que el índice son los barcodes
rownames(t.annotation) <- t.annotation$barcode
t.annotation$barcode <- NULL

# Agregar metadata al objeto global
data <- AddMetaData(data, t.annotation)



# Crear la nueva columna como character para que permita nuevos valores
data$Clustering.Round3 <- as.character(data$Clustering.Round2)

# Reasignar según tus subclusters refinados
data$Clustering.Round3[data$Subclustering.T == "Cd4 Activated"] <- "Cd4 Activated"
data$Clustering.Round3[data$Subclustering.T == "Cd4 Naive"] <- "Cd4 Naive"
data$Clustering.Round3[data$Subclustering.T == "Cd4 Th17"] <- "Cd4 Th17"
data$Clustering.Round3[data$Subclustering.T == "Cd4 Treg"] <- "Cd4 Treg"
data$Clustering.Round3[data$Subclustering.T == "Cd8 Effector"] <- "Cd8 Effector"
data$Clustering.Round3[data$Subclustering.T == "Mki67+ Tcell"] <- "Mki67+ Tcell"
data$Clustering.Round3[data$Subclustering.T == "NK"] <- "NK"

# Convertimos de vuelta a factor, ya con todos los niveles presentes
data$Clustering.Round3 <- factor(data$Clustering.Round3)



data$type <- paste(data$group, data$DsRed)


saveRDS(data, paste0(rsdir,"objects/data.Clusterized.Round3.rds"))
data <- readRDS(paste0(rsdir,"objects/data.Clusterized.Round3.rds"))

#### Plots 

data <- SetIdent(data, value = "Clustering.Round3") 

pdf(paste0(outdir,"/Clustering.Round3/Umap.pdf"), width=18, height=8)
 DimPlot(data, reduction = "umap", label = FALSE, raster = FALSE, cols = nora.colors)

dev.off()





# Downsampling

set.seed(123)

n_cells <- 8000
cells_to_keep <- c()

for (group in unique(data$tag)) {
  # Subset de células por metadata
  group_cells <- rownames(data@meta.data)[data@meta.data$tag == group]
  
  if (length(group_cells) <= n_cells) {
    sampled_cells <- group_cells
  } else {
    sampled_cells <- sample(group_cells, n_cells)
  }
  
  cells_to_keep <- c(cells_to_keep, sampled_cells)
}

# Ahora sí crear el subset
data_downsampled <- subset(data, cells = cells_to_keep)

# Revisar número de células por grupo
table(data_downsampled$tag)



pdf(paste0(outdir,"/Clustering.Round3/Umap.split.downsampled.pdf"), width=22, height=8)
 DimPlot(data, reduction = "umap", split.by="tag", 
 label = FALSE, raster = FALSE, cols = nora.colors, ncol = 4)

dev.off()



pdf(paste0(outdir,"/Clustering.Round3/Umap.split.downsampled2.8kcells.pdf"), width=18, height=10)
 DimPlot(data, reduction = "umap", split.by="type", 
 label = FALSE, raster = FALSE, cols = nora.colors, ncol = 2)

dev.off()






# Quantification


# Proportions por tag
prop_df <- data_downsampled@meta.data %>%
  group_by(type, Clustering.Round3) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(type) %>%
  mutate(prop = n / sum(n))

# Order levels of Clustering.Round2
prop_df$Clustering.Round3 <- factor(
  prop_df$Clustering.Round3,
  levels = names(nora.colors)
)



pdf(paste0(outdir,"/Clustering.Round3/quantification.type.pdf"), width=24, height=12)
ggplot(prop_df, aes(x = type, y = prop, fill = Clustering.Round3)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_manual(values = nora.colors) +
  theme_minimal(base_size = 14) +
  labs(
    x = "Tag",
    y = "Cell proportion (Downsampled to 8k)",
    fill = "Clustering"   # ← aquí cambias la leyenda
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()






prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round3",
	sample_1 = "KO DsRedN", sample_2 = "KO DsRedP",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round3/Proportion.test.KO.pdf"), width=8, height=10)
plot1 <- permutation_plot(prop_test) + ggtitle("Proportion test.DsRedP vs DsRedN. KO samples")
plot1
dev.off()



prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round3",
	sample_1 = "WT DsRedN", sample_2 = "WT DsRedP",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round3/Proportion.test.WT.pdf"), width=8, height=10)
plot2 <- permutation_plot(prop_test) + ggtitle("Proportion test.DsRedP vs DsRedN. WT samples")
plot2
dev.off()




prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round3",
	sample_1 = "WT DsRedP", sample_2 = "KO DsRedP",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round3/Proportion.test.Positive.KOvsWT.pdf"), width=8, height=10)
plot3 <- permutation_plot(prop_test) + ggtitle("Proportion test.KO vs WT. DsRedP samples")
plot3
dev.off()



prop_test <- sc_utils(data_downsampled)

prop_test <- permutation_test(
	prop_test, cluster_identity = "Clustering.Round3",
	sample_1 = "WT DsRedN", sample_2 = "KO DsRedN",
	sample_identity = "type"
)


pdf(paste0(outdir,"/Clustering.Round3/Proportion.test.Negative.KOvsWT.pdf"), width=8, height=10)
plot4 <- permutation_plot(prop_test) + ggtitle("Proportion test.KO vs WT. DsRedN samples")
plot4
dev.off()



pdf(paste0(outdir,"/Clustering.Round3/Proportion.tests.all.clustR3.pdf"), width=30, height=8)
plot_grid(plot1, plot2, plot3, plot4, ncol=4)
dev.off() 



## Bubble plots DE genes




macro_clusters <- c(
  "Ly6cHi Monocytes",
  "Ly6cLo Monocytes",
  "Early IFN|MHCII-TAMs",
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

Idents(data) <- "Clustering.Round3"
mac <- subset(
  data,
  idents = macro_clusters
)


mac <- subset(
  mac,
  subset = tag %in% c("WT1-DsRedP", "DsRedP-KO2")
)



mac$cluster_condition <- paste(
  mac$Clustering.Round3,
  mac$tag,
  sep = "_"
)

Idents(mac) <- "cluster_condition"


Idents(mac) <- "Clustering.Round3"

Idents(data) <- "Clustering.Round3"

markers_round3 <- FindAllMarkers(
  mac,
  only.pos = TRUE,
  min.pct = 0.3,
  logfc.threshold = 0.4,
  test.use = "MAST",
  latent.vars = "nCount_RNA"
)


library(dplyr)

top30_markers <- markers_round3 %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  ungroup()

as.data.frame(top30_markers)

write.table(top30_markers, paste0(rsdir,"table.clusters.markers.filt.Clusterized.R3.MAC.tsv"), sep='\t')
write.table(markers_round3, paste0(rsdir,"table.clusters.markers.Clusterized.R3.MAC.tsv"), sep='\t')

markers_round3 <- read.table(paste0(rsdir,"table.clusters.markers.Clusterized.R3.MAC.tsv"), sep='\t', header=T)



library(dplyr)
library(tibble)
set.seed(123)

Idents(mac) <- "Clustering.Round3"

de_by_cluster <- lapply(
  levels(mac),
  function(cl){

    cells_cl <- WhichCells(mac, idents = cl)
    meta_cl  <- mac@meta.data[cells_cl, ]

    n1 <- sum(meta_cl$tag == "DsRedP-KO2")
    n2 <- sum(meta_cl$tag == "WT1-DsRedP")

    if(min(n1, n2) < 20) return(NULL)

    ratio <- max(n1, n2) / min(n1, n2)

    if(ratio > 2){
      n_use <- min(n1, n2)
      cells_use <- c(
        sample(rownames(meta_cl)[meta_cl$tag == "DsRedP-KO2"], n_use),
        sample(rownames(meta_cl)[meta_cl$tag == "WT1-DsRedP"], n_use)
      )
    } else {
      cells_use <- rownames(meta_cl)
    }

    # 🔑 crear objeto downsampleado
    mac_ds <- subset(mac, cells = cells_use)

    FindMarkers(
      mac_ds,
      group.by = "tag",
      ident.1 = "DsRedP-KO2",
      ident.2 = "WT1-DsRedP",
      assay = "RNA",
      test.use = "MAST",
      latent.vars = c("nCount_RNA", "percent.mt"),
      min.pct = 0,
      logfc.threshold = 0
    ) %>%
      rownames_to_column("gene") %>%
      mutate(
        cluster = cl,
        n_cells = length(cells_use),
        ratio = ratio
      )
  }
) |> bind_rows()






de_filtered <- de_by_cluster %>%
  filter(
    p_val_adj < 0.05,         # significancia
    pct.1 > 0.4 | pct.2 > 0.4,  # expresión mínima en al menos 20% de células de algún grupo
    abs(avg_log2FC) > 0.25    # fold change relevante
  )




top30_markers <- de_filtered %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  ungroup()

  as.data.frame(top30_markers)

write.table(top30_markers, paste0(rsdir,"table.DE.KOvsWT.by.cluster.Clusterized.R3.MAC.Batch.A.tsv"), sep='\t')
markers <-read.table(paste0(rsdir,"table.DE.KOvsWT.by.cluster.Clusterized.R3.MAC.Batch.A.tsv"), sep='\t', header=T) 



# Pre-filtrar genes para que no sean escasos
expr <- GetAssayData(mac, slot="data")
genes_use <- rownames(expr)[
  rowSums(expr[, mac$tag=="DsRedP-KO2"] > 0) / sum(mac$tag=="DsRedP-KO2") >= 0.3 |
  rowSums(expr[, mac$tag=="DsRedN-KO2"] > 0) / sum(mac$tag=="DsRedN-KO2") >= 0.3
]

# DE rápido por clúster
Idents(mac) <- "Clustering.Round3"

de_by_cluster <- lapply(
  levels(mac),
  function(cl){
    FindMarkers(
      mac,
      subset.ident = cl,
      group.by = "tag",
      ident.1 = "DsRedP-KO2",
      ident.2 = "DsRedN-KO2",
      test.use = "wilcox",
      features = genes_use
    ) %>%
      tibble::rownames_to_column("gene") %>%
      mutate(cluster = cl)
  }
) |> bind_rows()

head(de_by_cluster, n=100)


de_filtered <- de_by_cluster %>%
  filter(
    p_val_adj < 0.05,         # significancia
    pct.1 > 0.4 | pct.2 > 0.4,  # expresión mínima en al menos 20% de células de algún grupo
    abs(avg_log2FC) > 0.25    # fold change relevante
  )


top30_markers <- de_filtered  %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 20) %>%
  ungroup()

  as.data.frame(top30_markers)


##### Scratch 


DefaultAssay(data) <- "SCT"


citoquinas_macrofagos <- c(
  # ILs
  "Il1a", "Il1b", "Il6", "Il10",
  "Il12a", "Il12b", "Il18", "Il23a",
  "Il19", "Il20",
  # Receptores IL
  "Il1r1", "Il1r2", "Il6ra",
  "Il10ra", "Il10rb",
  "Il12rb1", "Il12rb2",
  "Il18r1", "Il18rap",
  # No IL – TNF family
  "Tnf", "Lta", "Tnfsf9", "Tnfsf12",
  # No IL – Quimiocinas
  "Ccl2", "Ccl3", "Ccl4", "Ccl5", "Ccl7",
  "Cxcl2", "Cxcl9", "Cxcl10", "Cxcl11",
  # Otros mediadores inflamatorios
  "Ifnb1", "Osm", "Lif", "Spp1"
)




genes_interes <- c(
  "Il1a",
  "Il1b",
  "Il12a",
  "Il6",
  "Tnf",
  "Cxcl2",
  "Nfkbia",
  "Il10",
  "Ccl2"
)



pdf(paste0(outdir,"/Clustering.Round3/cytokines.pdf"), width=22, height=12)

VlnPlot(
  data,
  features = genes_interes,
  split.by = "group",
  pt.size = 0,
  cols = c("red", "blue"),
  ncol = 4,
  raster = FALSE
)

dev.off()


p <- VlnPlot(
  data,
  features = genes_interes,
  split.by = "group",
  group.by = "Clustering.Round3",
  pt.size = 0,
  ncol = 4,
  raster = FALSE
) + 
  scale_fill_manual(values = c("WT" = "blue", "KO" = "red")) +
  labs(fill = "group")   # <- Esto agrega la leyenda


pdf(paste0(outdir,"/Clustering.Round3/cytokines.pdf"), width=22, height=10)
p
dev.off()



pdf(paste0(outdir,"/Clustering.Round3/cytokines.pdf"), width=22, height=12)
plots <- VlnPlot(object = data, 
features = genes_interes, split.by = "group", group.by = "Clustering.Round3", 
pt.size = 0, combine = T, split.plot=T, log=T,raster = FALSE)
dev.off()





data$group <- factor(data$group, levels = c("WT", "KO"))



library(patchwork)
library(ggplot2)

plots <- VlnPlot(
  data,
  features = genes_interes,
  split.by = "group",
  group.by = "Clustering.Round3",
  pt.size = 0,
  ncol = 4,         # Esto controla los violines dentro de cada facet
  raster = FALSE
)

# Añadir colores y leyenda FUERA del panel
plots <- lapply(
  plots,
  function(p) {
    p +
      scale_fill_manual(values = c("WT" = "blue", "KO" = "red")) +
      labs(fill = "Group") +
      guides(fill = guide_legend(override.aes = list(color = NA))) +
      theme(
        legend.position = "right",         # Leyenda a la derecha
        legend.key.size = unit(1.5, "lines"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 16, face = "bold")
      )
  }
)

# Guardar PDF con 3 columnas
pdf(paste0(outdir,"/Clustering.Round3/cytokines.pdf"), width=20, height=16)
wrap_plots(plots, ncol = 3, guides = "collect")  # ncol=3 aquí
dev.off()







pdf(paste0(outdir,"/Clustering.Round3/cytokines2.pdf"), width=42, height=12)

VlnPlot(
  data,
  features = "Il1a",
  split.by = "tag",
  pt.size = 0,
  raster = FALSE
)

dev.off()

##### Markers

png(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.png"), width=1600, height=1400)
FeaturePlot_scCustom(data, features= c("Ifit3", "Ifit2", "Irf7"), reduction = "umap",
split.by="group")
dev.off()


png(paste0(outdir,"/Clustering.Round3/Umap_Markers2.png"), width=1600, height=1400)
FeaturePlot_scCustom(data, features= c("Mrc1", "Siglec1", "Folr2"), reduction = "umap",
split.by="group")
dev.off()



png(paste0(outdir,"/Clustering.Round3/cytokines3.png"), width=1800, height=1200)

VlnPlot(
  data,
  features = c("Il1a", "Il1b", "Il6", "Tnf", "Arg1", "Cxcl2"),
  group.by = "type",
  pt.size = 0,
  raster = FALSE
)

dev.off()




# 1) Receptores de interferón
ifn_receptores <- c(
  "Ifnar1",
  "Ifnar2",
  "Ifngr1",
  "Ifngr2",
  "Ifnlr1",
  "Il10rb"
)

# 2) Genes productores de interferón
ifn_productores <- c(
  # IFN tipo I
  "Ifna1",
  "Ifna2",
  "Ifna4",
  "Ifna5",
  "Ifna6",
  "Ifna7",
  "Ifna8",
  "Ifna12",
  "Ifna13",
  "Ifnb1",
  "Ifne",
  
  # IFN tipo II
  "Ifng",
  
  # IFN tipo III
  "Ifnl2",
  "Ifnl3"
)

# 3) Genes de respuesta a interferón (ISGs)
ifn_signaling <- c(
  # Señalización / regulación
  "Stat1",
  "Stat2",
  "Irf1",
  "Irf7",
  "Irf9")
  
  # Efectores clásicos
  ifn_effector <- c(
  "Mx1",
  "Oas1a",
  "Oas2",
  "Isg15",
  "Ifit1",
  "Ifit2",
  "Ifit3",
  "Rsad2",
  "Usp18",
  "Eif2ak2")
  
  # Familia interferon-inducible
  ifn_inducible <- c(
  "Ifi202",
  "Ifi27",
  "Ifi44",
  "Ifi47")
  
  # Inmunomodulación / presentación de antígeno
  ifn_modulation <- c(
  "Cxcl10",
  "Icam1",
  "Tap1",
  "Tap2",
  "B2m"
)

Siglec_genes <- c("Siglec1", "Adgre1", "Timd4", "Lyve1", "Mrc1")


maria <- c("Ifnar1",
"Ifnar2",
"Stat1",
"Stat2",
"Irf7",
"Irf9",
"Ifit1",
"Ifit2",
"Ifit3",
"Isg15",
"Siglec1",
"Cxcl10")




pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.R.pdf"), width=18, height=22)
FeaturePlot_scCustom(data, features= ifn_receptores, reduction = "umap",
split.by="type")
dev.off()


pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.R.group.pdf"), width=12, height=22)
FeaturePlot_scCustom(data, features= ifn_receptores, reduction = "umap",
split.by="group")
dev.off()

data$group <- factor(data$group, levels = c("WT", "KO"))


pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.genes.group.pdf"), width=12, height=42)
FeaturePlot_scCustom(data, features= maria, reduction = "umap",
split.by="group")
dev.off()







pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.P.pdf"), width=18, height=22)
FeaturePlot_scCustom(data, features= ifn_productores, reduction = "umap",
split.by="type")
dev.off()


png(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.RES.png"), width=2400, height=7800)
FeaturePlot_scCustom(data, features= ifn_respuesta, reduction = "umap",
split.by="type")
dev.off()

pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.signaling.pdf"), width=24, height=28)
FeaturePlot_scCustom(data, features= ifn_signaling, reduction = "umap",
split.by="type")
dev.off()


pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.signaling.group.pdf"), width=12, height=28)
FeaturePlot_scCustom(data, features= ifn_signaling, reduction = "umap",
split.by="group")
dev.off()

pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.effector.type.pdf"), width=22, height=36)
FeaturePlot_scCustom(data, features= ifn_effector, reduction = "umap",
split.by="type")
dev.off()


pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.effector.group.pdf"), width=11, height=36)
FeaturePlot_scCustom(data, features= ifn_effector, reduction = "umap",
split.by="group")
dev.off()



pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.inducible.type.pdf"), width=22, height=16)
FeaturePlot_scCustom(data, features= ifn_inducible, reduction = "umap",
split.by="type")
dev.off()


pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.inducible.group.pdf"), width=11, height=16)
FeaturePlot_scCustom(data, features= ifn_inducible, reduction = "umap",
split.by="group")
dev.off()

pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.modulation.type.pdf"), width=22, height=16)
FeaturePlot_scCustom(data, features= ifn_modulation, reduction = "umap",
split.by="type")
dev.off()


pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.modulation.group.pdf"), width=11, height=16)
FeaturePlot_scCustom(data, features= ifn_modulation, reduction = "umap",
split.by="group")
dev.off()

pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.siglec.type.pdf"), width=22, height=20)
FeaturePlot_scCustom(data, features= Siglec_genes, reduction = "umap",
split.by="type")
dev.off()


pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.siglec.group.pdf"), width=11, height=20)
FeaturePlot_scCustom(data, features= Siglec_genes, reduction = "umap",
split.by="group")
dev.off()



pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.siglec.group.pdf"), width=11, height=20)
FeaturePlot_scCustom(data, features= Siglec_genes, reduction = "umap",
split.by="group")
dev.off()


pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.Tfam.group.pdf"), width=18, height=10)
FeaturePlot_scCustom(data, features= "Tfam", reduction = "umap",
split.by="group")
dev.off()


pdf(paste0(outdir,"/Clustering.Round3/Umap_Markers.IFN.Tfam.type.pdf"), width=32, height=10)
FeaturePlot_scCustom(data, features= "Tfam", reduction = "umap",
split.by="type")
dev.off()


DefaultAssay(data) <- "RNA"

Idents(data) <- "Clustering.Round3"

pdf(
  paste0(outdir, "/Clustering.Round3/Violin_Markers.IFN.Tfams.pdf"),
  width = 18,
  height = 28
)

VlnPlot(
  data,
  features = c("Tfam","mt-Co1","mt-Nd1"),
  split.by = "type",
  raster=FALSE,
  ncol=1,
  pt.size=0.05
)

dev.off()





pdf(
  paste0(outdir, "/Clustering.Round3/Violin_Markers.IFN.Receptors.pdf"),
  width = 18,
  height = 10
)

VlnPlot(
  data,
  features = ifn_receptores,
  pt.size = 0,
  stack = TRUE,
  flip = TRUE,
  split.by = "type"
)

dev.off()





pdf(
  paste0(outdir, "/Clustering.Round3/Violin_Markers.IFN.signaling.pdf"),
  width = 18,
  height = 10
)

VlnPlot(
  data,
  features = ifn_signaling,
  pt.size = 0,
  stack = TRUE,
  flip = TRUE,
  split.by = "type"
)

dev.off()




pdf(
  paste0(outdir, "/Clustering.Round3/Violin_Markers.IFN.effector.pdf"),
  width = 18,
  height = 14
)

VlnPlot(
  data,
  features = ifn_effector,
  pt.size = 0,
  stack = TRUE,
  flip = TRUE,
  split.by = "type"
)

dev.off()



pdf(
  paste0(outdir, "/Clustering.Round3/Violin_Markers.IFN.inducible.pdf"),
  width = 18,
  height = 10
)

VlnPlot(
  data,
  features = ifn_inducible,
  pt.size = 0,
  stack = TRUE,
  flip = TRUE,
  split.by = "type"
)

dev.off()




pdf(
  paste0(outdir, "/Clustering.Round3/Violin_Markers.IFN.siglec.pdf"),
  width = 18,
  height = 10
)

VlnPlot(
  data,
  features = Siglec_genes,
  pt.size = 0,
  stack = TRUE,
  flip = TRUE,
  split.by = "type"
)

dev.off()





genes_ifnI <- c("Irf7", "Stat1", "Stat2", "Isg15", "Ifit1")


pdf(
  paste0(outdir, "/Clustering.Round3/Violin_IFN_I_Core.pdf"),
  width = 18,
  height = 10
)

VlnPlot(
  data,
  features = genes_ifnI,
  pt.size = 0,
  stack = TRUE,
  flip = TRUE,
  same.y.lims = TRUE,
  split.by = "type"
)

dev.off()







### Bubble plots for clustyering round 3

## DEG


# Clusters of interest
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
    "MHCII|Ccl12 Mac"
)

macros <- subset(data, subset = Clustering.Round3 %in% macro_clusters)


Idents(macros) <- "Clustering.Round3"


markers <- FindAllMarkers(
  object = macros,
  only.pos = TRUE,
  min.pct = 0.2,
  logfc.threshold = 0.25
)

head(markers)


top30_per_cluster <- markers %>%
  group_by(cluster) %>%
  slice_min(order_by = p_val_adj, n = 30, with_ties = FALSE) %>%
  ungroup()


as.data.frame(top30_per_cluster)


write.table(markers, paste0(rsdir,"table.clusters.markers.Clusterized.R3.tsv"), sep='\t')


write.table(top30_per_cluster, paste0(rsdir,"table.clusters.markers.Clusterized.top30.R3.tsv"), sep='\t')



### Bubble plot for markers


clusters<- read.table(paste0(rsdir,"table.clusters.markers.Clusterized.top30.R3.tsv"), sep='\t', header=T)


# Epsilon to avoid zero multiplication
epsilon <- 1e-8

# New column 'score'
clusters <- clusters %>%
  mutate(score = avg_log2FC * abs(pct.1 - pct.2 + epsilon))



clusters_macro <- clusters %>%
  filter(cluster %in% macro_clusters)

clusters_otros <- clusters %>%
  filter(!cluster %in% macro_clusters)


macros <- subset(data, subset = Clustering.Round3 %in% macro_clusters)

others <- subset(data, subset = Clustering.Round2 %in% macro_clusters, invert = TRUE)





# For 5 genes per cluster
top_markers_subset <- clusters_macro %>%
  group_by(cluster) %>%
  slice_max(order_by = score, n = 4) %>%
  pull(gene)


top_markers_subset <- unique(top_markers_subset)

genes <- c(
  "F10","Slc7a2","Mmp12","Inhba","Spp1","Ms4a4c","Spon1","Oas3","Irf7",
  "Isg15","Emp1","Arg1","Fn1","Vegfa","F13a1","Slc27a1","Ifit3","Ifit2",
  "Oas2","Ifi44","Usp18","Plac8","Slfn4","Ly6c2","Sell","Rsad2","Gsr",
  "Hp","Cd244a","Vcan","Adgre5","Ciita","Tent5c","Slamf7","H2-DMb1",
  "Nr4a3","H2-DMa","Cp","Siglec1","Siglece","Ighm","Icosl","Mmp9",
  "Atp6v0d2","Ctsk","Actn1","Coq8a","Gas6","Mrc1","Cbr2","Fcrls",
  "C1qc","Npr2","Stmn1","Igf2r","Sdc1","Trem1","Ptgs2","Cd80",
  "Plaur","Clec4e"
)



# For 5 genes per cluster
top_markers_subset2 <- clusters_otros  %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  pull(gene)


cluster_order <- c(
  "Ly6cHi Monocytes", 
  "Ly6cLo Monocytes", 
  "Early IFN|MHCII-TAMs",
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


macros$Clustering.Round3 <- factor(macros$Clustering.Round3)
macros$Clustering.Round3 <- droplevels(macros$Clustering.Round3)
table(macros$Clustering.Round3)


macros$Clustering.Round3 <- factor(macros$Clustering.Round3, levels = cluster_order)
Idents(macros) <- "Clustering.Round3"


# --- Plots---

pdf(paste0(outdir, "/Clustering.Round3/bubble.macros.pdf"), width = 12, height = 14)
p <- Clustered_DotPlot(
  seurat_object = macros,
  features = genes,
  colors_use_idents = nora.colors2,
  k=5
)

print(p[[1]])
dev.off()



# ===============================
# Libraries
# ===============================
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# ===============================
# Colores de clusters (FIJOS)
# ===============================
nora.colors2 <- c(
  "Ly6cHi Monocytes"              = "#FF0000",
  "Ly6cLo Monocytes"              = "#FF6A6A",
  "Early IFN|MHCII-TAMs"           = "#B22222",
  "Trem1|Ptgs2|Plaur|Celc4e Mac"   = "#EE7942",
  "Mrc1|C1qc|Cbr2|Gas6 Mac"        = "#FFD92F",
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac" = "#4DAF4A",
  "Npr2|Actn1 Mac"                = "#A6D854",
  "Neutrophils"                   = "#4876FF",
  "Mmp9|Ctsk Mac"                 = "#00723F",
  "IFN Mac"                       = "#C080FF",
  "Fn1|Vegfa Mac"                 = "#FFA500",
  "MHCII|Siglec Mac"              = "#1E90FF",
  "MHCII|Ccl12 Mac"               = "#4682B4"
)

# ===============================
# Genes
# ===============================
genes <- unique(c(
  "F10","Slc7a2","Mmp12","Inhba","Spp1","Ms4a4c","Spon1","Oas3","Irf7",
  "Isg15","Emp1","Arg1","Fn1","Vegfa","F13a1","Slc27a1","Ifit3","Ifit2",
  "Oas3","Ifi44","Usp18","Plac8","Slfn4","Ly6c2","Sell","Rsad2","Gsr",
  "Hp","Cd244a","Vcan","Adgre5","Ciita","Tent5c","Slamf7","H2-DMb1",
  "Nr4a3","H2-DMa","Cp","Siglec1","Siglece","Ighm","Icosl","Mmp9",
  "Atp6v0d2","Ctsk","Actn1","Coq8a","Gas6","Mrc1","Cbr2","Fcrls",
  "C1qc","Nrp2","Stmn1","Igf2r","Sdc1","Trem1","Ptgs2","Cd80",
  "Plaur","Clec4e"
))




# ===============================
# Generar Clustered_DotPlot
# ===============================
pdf(paste0(outdir, "/Clustering.Round3/bubble_macros_clustered.pdf"),
    width = 12, height = 14)

p <- Clustered_DotPlot(
  seurat_object = macros,
  features = genes,
  colors_use_idents=nora.colors2,
  k = 5  # número de clusters jerárquicos de genes
)

# Imprimir plots
print(p[[1]])  # bubble plot principal
print(p[[2]])  # heatmap lateral
dev.off()





clusters_data <- c(
  "Activated B cells",
  "B cells",
  "Cd4 Activated",
  "Cd4 Naive",
  "Cd4 Th17",
  "Cd4 Treg",
  "Cd8 Effector",
  "DCs",
  "Mastocytes",
  "Neutrophils",
  "NK"
)




# --- Plots ---
pdf(paste0(outdir, "/Analysis/bubble.macros.k5.pdf"), width=12, height=14)

p <- Clustered_DotPlot(seurat_object = macros, features = top_markers_subset,
colors_use_idents = nora.colors2)


print(p[[1]])
dev.off()



nora.colors3 <- nora.colors[!names(nora.colors) %in% names(nora.colors2)]


# 1️⃣ Identificar los clusters con >0 células
clusters_present <- names(table(others$Cluster))[table(others$Cluster) > 0]

# 2️⃣ Crear nora.colors3 con solo esos clusters
nora.colors3 <- nora.colors[clusters_present]

# 3️⃣ Revisar
nora.colors3




Idents(others) <- factor(Idents(others), levels = names(nora.colors3))

# --- Graficar y exportar ---
pdf(paste0(outdir, "/Analysis/bubble.other.pdf"), width=10, height=10)

p <- Clustered_DotPlot(seurat_object = others, features = top_markers_subset2,
colors_use_idents = nora.colors3)


print(p[[1]])
dev.off()


# --- Graficar y exportar ---
pdf(paste0(outdir, "/Analysis/bubble.other.k4.pdf"), width=12, height=14)

p <- Clustered_DotPlot(seurat_object = others, features = top_markers_subset2,
colors_use_idents = nora.colors3, k=4)


print(p[[1]])
dev.off()




# Filtrar nora.colors para que solo tenga los clusters presentes en 'macros'
nora.colors.filtered <- nora.colors[names(nora.colors) %in% levels(Idents(macros))]

# Asegurarnos de que los niveles del objeto Seurat estén en el mismo orden
Idents(macros) <- factor(Idents(macros), levels = names(nora.colors.filtered))


pdf(paste0(outdir, "/Analysis/bubble.macros.pdf"), width=12, height=14)

p <- Clustered_DotPlot(seurat_object = macros, features = top_markers_subset)

# Alinear colores con los niveles actuales del objeto
p[[1]] <- p[[1]] +
  scale_color_manual(values = nora.colors[levels(macros)]) +
  scale_fill_manual(values = nora.colors[levels(macros)])

print(p[[1]])
dev.off()




# 2️⃣ Seleccionar los 8 principales genes por cluster
top8_markers <- clusters %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 6) %>%
  ungroup()

# Eliminar duplicados de la lista de genes
top8_genes <- unique(top8_markers$gene)

# 3️⃣ Crear el dotplot sin duplicados
dotplot <- DotPlot(
  object = macros,
  features = top_markers_subset
) +
  RotatedAxis() +
  scale_color_gradient(low = "lightgrey", high = "red") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))




pdf(paste0(outdir,"/Analysis/bubble1.pdf"), width=36, height=12)
print(dotplot)
dev.off()
