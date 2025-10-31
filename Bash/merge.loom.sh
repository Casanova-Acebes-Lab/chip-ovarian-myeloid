#!/bin/sh
#SBATCH --job-name=merge_velocity
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/home/dcaceres/Tet2sc/Bash/merge_velocity_%A.out
#SBATCH --error=/home/dcaceres/Tet2sc/Bash/merge_velocity_%A.err

# ============================================================
# 🚀 MERGE DE VELOCYTO LOOM FILES + ADATA BASE
# Limpia barcodes dejando SOLO la secuencia de nucleótidos (ACGT),
# asigna capas spliced/unspliced/ambiguous a las celdas del adata subset
# ============================================================

echo "------------------------------------------------------------"
echo "🚀 INICIO DEL JOB MERGE VELOCYTO - $(date)"
echo "------------------------------------------------------------"

# --- Activar entorno con scvelo ---
source /home/dcaceres/miniconda3/etc/profile.d/conda.sh
conda activate scvelo_env

# --- Directorios principales ---
OUTDIR="/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/processed_h5ad"
mkdir -p "$OUTDIR"

PYTMP="$OUTDIR/merge_velocity_tmp.py"

cat << 'EOF' > "$PYTMP"
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Merge automático de loom files sobre adata original (de Seurat),
limpiando los barcodes para quedarse SOLO con la secuencia de nucleótidos (A/C/G/T).
"""

import os, sys, re
from datetime import datetime
import scanpy as sc
import anndata as ad
import numpy as np

OUT_DIR = "/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/processed_h5ad"
os.makedirs(OUT_DIR, exist_ok=True)
LOG_FILE = os.path.join(OUT_DIR, "merge_velocity.log")
sys.stdout = open(LOG_FILE, "w")
sys.stderr = sys.stdout

print("------------------------------------------------------------")
print("🚀 INICIO DEL PROCESO MERGE -", datetime.now())
print("------------------------------------------------------------\n")

# --- Cargar adata base ---
adata_path = "/home/dcaceres/Tet2sc/Results/plots/scVelo/my_data.h5ad"
print(f"Cargando adata base desde: {adata_path}")
adata = sc.read_h5ad(adata_path)
adata.var_names_make_unique()
print(f"✅ adata original: {adata.n_obs} células x {adata.n_vars} genes\n")

# Evitar conflicto 'barcode' si existe como columna
if "barcode" in adata.obs.columns:
    print("🧹 Eliminando columna duplicada 'barcode' de adata.obs para evitar conflictos.")
    adata.obs = adata.obs.drop(columns=["barcode"])

# --- Rutas de loom ---
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

# --- Cargar looms ---
ldata_list = []
for path in loom_paths:
    print(f"Leyendo: {path}")
    ld = sc.read_loom(path)
    if "matrix" in ld.layers and ld.layers["matrix"].shape != ld.X.shape:
        del ld.layers["matrix"]
    ld.var_names_make_unique()
    ldata_list.append(ld)
print("\n✅ Todos los loom cargados correctamente")

# --- Alinear genes comunes ---
common_genes = set(adata.var_names)
for ld in ldata_list:
    common_genes &= set(ld.var_names)
common_genes = sorted(common_genes)
print(f"✅ Genes comunes entre adata y loom: {len(common_genes)}")

for i, ld in enumerate(ldata_list):
    ld = ld[:, common_genes].copy()
    ld.obs["batch"] = batch_keys[i]
    ldata_list[i] = ld

ldata_merged = ad.concat(ldata_list, join="inner", label="batch", keys=batch_keys, index_unique=None)
ldata_merged.obs_names_make_unique()
print(f"✅ Loom merge completo: {ldata_merged.n_obs} células x {ldata_merged.n_vars} genes")

# --- Función para limpiar barcodes ---
_re_seq = re.compile(r'([ACGTacgt]{8,})')

def extract_sequence_only(bc):
    if not isinstance(bc, str):
        return bc
    bc = bc.strip()
    matches = _re_seq.findall(bc)
    if not matches:
        return bc
    return matches[-1].upper()

# Aplicar limpieza
print("\n🔬 Extrayendo solo la secuencia de nucleótidos de los barcodes...")
adata_seq = adata.obs_names.to_series().apply(extract_sequence_only)
ldata_seq = ldata_merged.obs_names.to_series().apply(extract_sequence_only)

def make_unique(series):
    vals = list(series)
    seen = {}
    out = []
    for v in vals:
        if v not in seen:
            seen[v] = 0
            out.append(v)
        else:
            seen[v] += 1
            newv = f"{v}_{seen[v]}"
            while newv in seen:
                seen[v] += 1
                newv = f"{v}_{seen[v]}"
            seen[newv] = 0
            out.append(newv)
    return out, sum(1 for c in vals if vals.count(c) > 1)

adata_seq_unique, ndup_adata = make_unique(adata_seq)
ldata_seq_unique, ndup_loom = make_unique(ldata_seq)

if ndup_adata > 0:
    print(f"⚠️ Atención: {ndup_adata} duplicados detectados tras limpiar adata. Se han renombrado añadiendo sufijos.")
if ndup_loom > 0:
    print(f"⚠️ Atención: {ndup_loom} duplicados detectados tras limpiar loom. Se han renombrado añadiendo sufijos.")

adata.obs_names = adata_seq_unique
ldata_merged.obs_names = ldata_seq_unique

print("\n🧩 Ejemplo de barcodes después de extraer secuencia (adata):")
print(list(adata.obs_names[:5]))
print("🧩 Ejemplo de barcodes después de extraer secuencia (loom):")
print(list(ldata_merged.obs_names[:5]))

# --- Intersección de celdas ---
common_cells = adata.obs_names.intersection(ldata_merged.obs_names)
print(f"\n✅ Celdas con capas de velocity disponibles: {len(common_cells)}")

# --- Subset y transferencia de capas ---
print(f"\n🧬 Subseteando adata (mantiene {adata.n_obs} células del análisis original).")
adata = adata.copy()
adata.layers["spliced"] = np.zeros(adata.shape, dtype=np.float32)
adata.layers["unspliced"] = np.zeros(adata.shape, dtype=np.float32)
adata.layers["ambiguous"] = np.zeros(adata.shape, dtype=np.float32)

if len(common_cells) > 0:
    print(f"🔗 Transfiriendo capas de {len(common_cells)} células comunes...")
    common_cells_sorted = sorted(common_cells)
    loom_sub = ldata_merged[common_cells_sorted, :].copy()

    # 💡 NUEVO: Alinear genes antes de usar adata.var_names
    common_gene_set = adata.var_names.intersection(loom_sub.var_names)
    adata = adata[:, list(common_gene_set)].copy()
    loom_sub = loom_sub[:, list(common_gene_set)].copy()
    print(f"✅ Genes alineados para transferencia: {len(common_gene_set)}")

    idx_adata = adata.obs_names.isin(common_cells_sorted)
    for layer in ["spliced", "unspliced", "ambiguous"]:
        if layer in loom_sub.layers:
            layer_data = loom_sub.layers[layer]
            if not isinstance(layer_data, np.ndarray):
                try:
                    layer_data = layer_data.A
                except Exception:
                    layer_data = layer_data.toarray()
            adata.layers[layer][idx_adata, :] = layer_data
else:
    print("⚠️ No se encontraron celdas comunes, capas vacías permanecerán en cero.")

out_file = os.path.join(OUT_DIR, "adata_with_velocity_layers_subsetted.h5ad")
adata.write(out_file)
print(f"\n✅ Archivo final guardado: {out_file}")

print("\n------------------------------------------------------------")
print("🎉 PROCESO FINALIZADO CORRECTAMENTE -", datetime.now())
print("------------------------------------------------------------")
sys.stdout.close()
EOF

# --- Ejecutar script Python ---
echo "🐍 Ejecutando merge_velocity_tmp.py dentro del job..."
python "$PYTMP"

echo "------------------------------------------------------------"
echo "✅ JOB FINALIZADO $(date)"
echo "------------------------------------------------------------"
