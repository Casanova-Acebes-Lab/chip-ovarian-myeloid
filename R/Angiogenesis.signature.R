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


# Read clusterized 1

data <- readRDS(paste0(rsdir,"/objects/data.Clustering.Round1.rds"))



## Isolating Macrophages and Monocytes

macmono_obj <- subset(
  data,
  subset = Clustering.wide %in% c("Macrophages", "IFN Mac", "Monocytes")
)

table(macmono_obj$Clustering.wide)
table(macmono_obj$Clustering.Round3)
table(macmono_obj$sample_id)


## Angiogenesis signature

wt_functional_angio_genes <- c(
  "Adm",
  "Hbegf",
  "Pgf",
  "Vegfa",
  "Flt1",
  "Nrp2",
  "Igf1",
  "Areg",
  "Nos2",
  "Hif1a",
  "Jag1",
  "Jag2",
  "Mmp9",
  "Mmp12",
  "Mmp19",
  "Ctsk",
  "Timp1",
  "Timp3",
  "Plau",
  "Plaur",
  "Vcam1",
  "Icam1",
  "Il6",
  "Il1a"
)

ko_imperfect_ecm_genes <- c(
  "Osm",
  "Tgfb3",
  "Tgfbr1",
  "Col4a1",
  "Col4a2",
  "Col5a1",
  "Col5a2",
  "Col6a1",
  "Col6a2",
  "Col6a3",
  "Col18a1",
  "Col3a1",
  "Tnc",
  "Lama4",
  "Lamb1",
  "Nid1",
  "Pdgfrb",
  "Flt4",
  "Kdr",
  "Aplnr",
  "Bgn",
  "Adamts2",
  "Cxcl12"
)


## Check for genes

wt_present <- intersect(wt_functional_angio_genes, rownames(macmono_obj))
ko_present <- intersect(ko_imperfect_ecm_genes, rownames(macmono_obj))

setdiff(wt_functional_angio_genes, wt_present)
setdiff(ko_imperfect_ecm_genes, ko_present)

length(wt_present)
length(ko_present)


## Add module scores


DefaultAssay(macmono_obj) <- "RNA"

macmono_obj <- AddModuleScore(
  object = macmono_obj,
  features = list(wt_present),
  name = "WT_Functional_Angio",
  assay = "RNA"
)

macmono_obj <- AddModuleScore(
  object = macmono_obj,
  features = list(ko_present),
  name = "KO_Imperfect_ECM",
  assay = "RNA"
)

head(macmono_obj@meta.data[, c("WT_Functional_Angio1", "KO_Imperfect_ECM1")])


## Plotting


library(scCustomize)
library(ggplot2)

# Asegurar orden correcto
macmono_obj$group <- factor(macmono_obj$group, levels = c("WT", "KO"))

p_wt_angio <- FeaturePlot_scCustom(
  seurat_object = macmono_obj,
  features = "WT_Functional_Angio1",
  reduction = "umap",
  split.by = "group",
  order = TRUE,
  colors_use = viridis_plasma_dark_high,
  na_color = "lightgray"
) 

p_ko_ecm <- FeaturePlot_scCustom(
  seurat_object = macmono_obj,
  features = "KO_Imperfect_ECM1",
  reduction = "umap",
  split.by = "group",
  order = TRUE,
  colors_use = viridis_plasma_dark_high,
  na_color = "lightgray"
)

p_wt_angio
p_ko_ecm

pb_dir <- file.path(outdir, "Pseudobulk_Macrophages_Monocytes_Round3")
dir.create(pb_dir, recursive = TRUE, showWarnings = FALSE)

score_dir <- file.path(pb_dir, "plots", "SingleCell_ModuleScores")
dir.create(score_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
  filename = file.path(score_dir, "UMAP_WT_functional_angiogenesis_score_split_WT_vs_KO_scCustomize.pdf"),
  plot = p_wt_angio,
  width = 12,
  height = 5
)

ggsave(
  filename = file.path(score_dir, "UMAP_KO_imperfect_ECM_score_split_WT_vs_KO_scCustomize.pdf"),
  plot = p_ko_ecm,
  width = 12,
  height = 5
)