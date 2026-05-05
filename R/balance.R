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



imb_dir <- file.path(outdir, "DsRed_imbalance_analysis")
dir.create(imb_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(imb_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(imb_dir, "plots"), recursive = TRUE, showWarnings = FALSE)



# Read clusterized 1

data <- readRDS(paste0(rsdir,"/objects/data.Clustering.Round1.rds"))





# -----------------------------
# 1. Comprobaciones iniciales
# -----------------------------

required_cols <- c(
  "sample_id", "batch", "group", "DsRed", "mouse_id",
  "tet2_status", "cell_context", "percent.mt",
  "nCount_RNA", "nFeature_RNA",
  "Clustering.Round3", "Clustering.wide"
)

missing_cols <- setdiff(required_cols, colnames(data@meta.data))
if (length(missing_cols) > 0) {
  stop("Faltan columnas en metadata: ", paste(missing_cols, collapse = ", "))
}

# Fijar nombres finales estables
data$celltype_final <- as.character(data$Clustering.Round3)
data$celltype_wide <- as.character(data$Clustering.wide)

data$celltype_final <- factor(data$celltype_final)
data$celltype_wide <- factor(data$celltype_wide)

data$sample_id <- factor(data$sample_id)
data$batch <- factor(data$batch, levels = c("A", "B", "C"))
data$group <- factor(data$group, levels = c("WT", "KO"))
data$DsRed <- factor(data$DsRed, levels = c("DsRedN", "DsRedP"))
data$tet2_status <- factor(data$tet2_status, levels = c("Tet2_WT", "Tet2_KO"))
data$cell_context <- factor(
  data$cell_context,
  levels = c(
    "WT_DsRedN_control",
    "WT_DsRedP_control",
    "KO_DsRedN_WT_bystander",
    "KO_DsRedP_true_Tet2KO"
  )
)

Idents(data) <- "celltype_final"

# Diseño final
sample_design <- data@meta.data %>%
  as.data.frame() %>%
  distinct(sample_id, batch, group, DsRed, mouse_id, tet2_status, cell_context) %>%
  arrange(batch, group, mouse_id, DsRed)

write.csv(
  sample_design,
  file.path(imb_dir, "tables", "sample_design.csv"),
  row.names = FALSE
)






# -----------------------------
# 2. QC de profundidad por muestra/contexto/DsRed
# -----------------------------

qc_by_sample <- data@meta.data %>%
  as.data.frame() %>%
  group_by(sample_id, batch, group, DsRed, mouse_id, tet2_status, cell_context) %>%
  summarise(
    n_cells = n(),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    mean_percent_mt = mean(percent.mt, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(group, mouse_id, DsRed)

  as.data.frame(qc_by_sample)

qc_by_dsred <- data@meta.data %>%
  as.data.frame() %>%
  group_by(DsRed) %>%
  summarise(
    n_cells = n(),
    median_nCount_RNA = median(nCount_RNA, na.rm = TRUE),
    mean_nCount_RNA = mean(nCount_RNA, na.rm = TRUE),
    median_nFeature_RNA = median(nFeature_RNA, na.rm = TRUE),
    mean_nFeature_RNA = mean(nFeature_RNA, na.rm = TRUE),
    median_percent_mt = median(percent.mt, na.rm = TRUE),
    .groups = "drop"
  )

  as.data.frame(qc_by_dsred)

write.csv(qc_by_sample, file.path(imb_dir, "tables", "qc_depth_by_sample.csv"), row.names = FALSE)
write.csv(qc_by_dsred, file.path(imb_dir, "tables", "qc_depth_by_DsRed.csv"), row.names = FALSE)

## Apreciamos el desbalance de counts en Dsred, pero no es uniforme entre muestras
## Puede tener un componente biológico, pero también puede ser un artefacto técnico. Por eso es importante analizarlo por muestra y no solo por DsRed



comp_wide <- data@meta.data %>%
  as.data.frame() %>%
  dplyr::group_by(
    sample_id,
    batch,
    group,
    DsRed,
    mouse_id,
    tet2_status,
    cell_context,
    Clustering.wide
  ) %>%
  dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(freq_total = 100 * n / sum(n)) %>%
  dplyr::ungroup()

as.data.frame(comp_wide)



# La fracción DsRedP muestra una depleción reproducible de células T,
#especialmente CD8, tanto en muestras WT como KO. 
#Por tanto, este efecto no puede atribuirse directamente a Tet2 KO. 
#Podría reflejar un sesgo técnico/experimental asociado a DsRedP, 
#al sorting o al propio sistema de marcado, o bien un fenómeno biológico 
#inducido por DsRedP que afecte a la entrada, persistencia o recuperación de 
#células CD8 en el TME.


# Miramos las Tcells 

tcell_counts_Tonly <- tcell_counts %>%
  dplyr::filter(!Tcells %in% c("NK", "ILC")) %>%
  dplyr::group_by(sample_id) %>%
  dplyr::mutate(
    n_Tonly = sum(n),
    freq_within_Tonly = 100 * n / n_Tonly
  ) %>%
  dplyr::ungroup()

as.data.frame(tcell_counts_Tonly)



tcell_summary_Tonly <- tcell_counts_Tonly %>%
  dplyr::group_by(cell_context, Tcells) %>%
  dplyr::summarise(
    mean_freq_total = mean(freq_total),
    sd_freq_total = sd(freq_total),
    mean_freq_within_Tonly = mean(freq_within_Tonly),
    sd_freq_within_Tonly = sd(freq_within_Tonly),
    .groups = "drop"
  )

as.data.frame(tcell_summary_Tonly)


#  La fracción DsRedP muestra una depleción marcada de células CD8, especialmente Cd8 Exhausted,
# tanto en WT como en KO. Este efecto se observa tanto como porcentaje del total celular como dentro 
# del compartimento T-only, por lo que no refleja solo una pérdida global de T cells, sino un sesgo 
# específico contra estados CD8 efector/exhausted. Dado que ocurre también en WT_DsRedP, se interpreta 
# como efecto asociado a DsRedP/sorting/recuperación experimental, no como evidencia primaria de efecto Tet2.