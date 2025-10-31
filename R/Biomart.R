 

# Enviroment R.4.3.3
library(Seurat)
library(dplyr)
library(ggplot2)
library(org.Mm.eg.db)
library(biomaRt)
library(org.Hs.eg.db)  # Humanos
library(AnnotationDbi)



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


### Read data

signature.Bul.Tet2 <- read.table(paste0('data/Signatures/Signature.Huerga_Encabo.H.2023.Bulk.tsv'),sep='\t', header=T)


