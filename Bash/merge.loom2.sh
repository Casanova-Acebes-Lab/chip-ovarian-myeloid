#!/bin/sh
#SBATCH --job-name=merge_velocity
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/home/dcaceres/Tet2sc/Bash/merge_velocity_%A.out
#SBATCH --error=/home/dcaceres/Tet2sc/Bash/merge_velocity_%A.err

# ============================================================
# 🚀 MERGE DE VELOCYTO LOOM FILES + ADATA BASE
# Script todo en uno listo para cluster con SLURM
# ============================================================

echo "------------------------------------------------------------"
echo "🚀 INICIO DEL JOB MERGE VELOCYTO - $(date)"
echo "------------------------------------------------------------"

# --- Activar entorno con scvelo ---
source /home/dcaceres/miniconda3/etc/profile.d/conda.sh
conda activate scvelo_env

# --- Directorio de salida ---
OUTDIR="/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/processed_h5ad"
mkdir -p "$OUTDIR"

# --- Ejecutar merge directamente con Python ---
python3 - << 'END_PYTHON'
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, sys
from datetime import datetime
import scanpy as sc
import anndata as ad
import numpy as np

# -------------------------------
# Configuración
# -------------------------------
OUT_DIR = "/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/processed_h5ad"
os.makedirs(OUT_DIR, exist_ok=True)
LOG_FILE = os.path.join(OUT_DIR, "merge_velocity.log")
sys.stdout = open(LOG_FILE, "w")
sys.stderr = sys.stdout

print("------------------------------------------------------------")
print("🚀 INICIO DEL PROCESO MERGE -", datetime.now())
print("------------------------------------------------------------\n")

# -------------------------------
# Cargar adata base
# -------------------------------
adata_path = "/home/dcaceres/Tet2sc/Results/plots/scVelo/my_data.h5ad"
print(f"Cargando adata base desde: {adata_path}")
adata = sc.read_h5ad(adata_path)
adata.var_names_make_unique()
print(f"✅ adata original: {adata.n_obs} células x {adata.n_vars} genes\n")

if "barcode" in adata.obs.columns:
    adata.obs = adata.obs.drop(columns=["barcode"])

# -------------------------------
# Rutas de loom
# -------------------------------
loom_paths = [
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedN-KO2.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedN-KO3.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedP-KO2.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedP-KO3.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT1-DsRedn.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT1-DsRedp.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT2-DsRedn.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT2-DsRedp.loom',
]
batch_keys = ["KO2n", "KO3n", "KO2p", "KO3p", "WT1n", "WT1p", "WT2n", "WT2p"]

# -------------------------------
# Cargar y alinear looms
# -------------------------------
ldata_list = []
for path in loom_paths:
    print(f"Leyendo: {path}")
    ld = sc.read_loom(path)
    ld.var_names_make_unique()
    ldata_list.append(ld)
print("\n✅ Todos los loom cargados correctamente")

# -------------------------------
# Genes comunes
# -------------------------------
common_genes = set(adata.var_names)
for ld in ldata_list:
    common_genes &= set(ld.var_names)
common_genes = sorted(common_genes)
print(f"✅ Genes comunes entre adata y loom: {len(common_genes)}")

# Subset de genes y asignación de batch
for i, ld in enumerate(ldata_list):
    ld = ld[:, common_genes].copy()
    ld.obs["batch"] = batch_keys[i]
    ldata_list[i] = ld

ldata_merged = ad.concat(ldata_list, join="inner", label="batch", keys=batch_keys, index_unique=None)
ldata_merged.obs_names_make_unique()
print(f"✅ Loom merge completo: {ldata_merged.n_obs} células x {ldata_merged.n_vars} genes")

# -------------------------------
# Limpiar barcodes (solo secuencia)
# -------------------------------
_re_seq = re.compile(r'([ACGTacgt]{8,})')

def extract_sequence_only(bc):
    if not isinstance(bc, str):
        return bc
    matches = _re_seq.findall(bc.strip())
    return matches[-1].upper() if matches else bc

adata.obs_names = adata.obs_names.to_series().apply(extract_sequence_only)
ldata_merged.obs_names = ldata_merged.obs_names.to_series().apply(extract_sequence_only)

# -------------------------------
# Hacer obs_names únicos
# -------------------------------
adata.obs_names_make_unique()
ldata_merged.obs_names_make_unique()

# -------------------------------
# Intersección de celdas
# -------------------------------
common_cells = adata.obs_names.intersection(ldata_merged.obs_names)
print(f"\n✅ Celdas con capas de velocity disponibles: {len(common_cells)}")

# -------------------------------
# Inicializar capas
# -------------------------------
for layer in ["spliced", "unspliced", "ambiguous"]:
    adata.layers[layer] = np.zeros(adata.shape, dtype=np.float32)

# -------------------------------
# Transferencia de capas correctamente alineada
# -------------------------------
if len(common_cells) > 0:
    print("🔗 Transfiriendo capas a celdas comunes...")
    loom_sub = ldata_merged[common_cells, :].copy()
    
    # Alinear genes
    common_gene_set = adata.var_names.intersection(loom_sub.var_names)
    adata = adata[:, list(common_gene_set)].copy()
    loom_sub = loom_sub[:, list(common_gene_set)].copy()
    print(f"✅ Genes alineados para transferencia: {len(common_gene_set)}")
    
    # Índices exactos
    adata_idx = adata.obs_names.get_indexer(common_cells)
    loom_idx = loom_sub.obs_names.get_indexer(common_cells)
    
    for layer in ["spliced", "unspliced", "ambiguous"]:
        if layer in loom_sub.layers:
            layer_data = loom_sub.layers[layer]
            if not isinstance(layer_data, np.ndarray):
                layer_data = layer_data.toarray()
            adata.layers[layer][adata_idx, :] = layer_data[loom_idx, :]
else:
    print("⚠️ No se encontraron celdas comunes, capas permanecerán en cero")

# -------------------------------
# Guardar h5ad final
# -------------------------------
out_file = os.path.join(OUT_DIR, "adata_with_velocity_layers_subsetted.h5ad")
adata.write(out_file)
print(f"\n✅ Archivo final guardado: {out_file}")

print("\n------------------------------------------------------------")
print("🎉 PROCESO FINALIZADO CORRECTAMENTE -", datetime.now())
print("------------------------------------------------------------")
sys.stdout.close()
END_PYTHON

echo "------------------------------------------------------------"
echo "✅ JOB FINALIZADO $(date)"
echo "------------------------------------------------------------"
