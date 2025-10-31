#!/bin/sh
#SBATCH -J scvelo_script
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 08:00:00
#SBATCH -o /home/dcaceres/Tet2sc/Bash/script.read.out
#SBATCH -e /home/dcaceres/Tet2sc/Bash/script.read.err

# ============================================================
# CONFIGURACIÓN DEL ENTORNO
# ============================================================
source /home/dcaceres/miniconda3/etc/profile.d/conda.sh
conda activate scvelo_env

python << 'EOF'

import os
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scv
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)

# --- Directorios ---
outdir = "/home/dcaceres/Tet2sc/Results/plots"
scvelo_outdir = os.path.join(outdir, "scVelo")
os.makedirs(scvelo_outdir, exist_ok=True)
sc.settings.figdir = scvelo_outdir

# --- Cargar AnnData ya procesado desde Seurat ---
adata = sc.read_h5ad(os.path.join(scvelo_outdir, "my_data.h5ad"))

# --- Loom files de Velocyto ---
loom_paths = [
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedN-KO2.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedN-KO3.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedP-KO2.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedP-KO3.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT1-DsRedn.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT1-DsRedp.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT2-DsRedn.loom',
    '/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT2-DsRedp.loom'
]

ldata_list = []

for i, path in enumerate(loom_paths, start=1):
    print(f"🔹 Cargando {path}")
    ld = sc.read_loom(path, obs_names='CellID', var_names='Gene')
    ld.var_names_make_unique()
    
    # --- Asegurar que spliced y unspliced estén en layers ---
    if 'spliced' not in ld.layers:
        ld.layers['spliced'] = ld.X.copy()
    if 'unspliced' not in ld.layers:
        # Intentar cargar unspliced si existe en X_name original
        ld.layers['unspliced'] = ld.X.copy()
    
    # Etiqueta de batch
    ld.obs['batch'] = f'sample_{i}'
    
    # Corregir barcodes para merge
    barcodes = [bc.split(':')[-1] for bc in ld.obs.index]
    ld.obs.index = [f"{bc}_{i}" for bc in barcodes]
    
    ldata_list.append(ld)

# --- Concatenar todos los loom files ---
ldata_merged = ad.concat(ldata_list, join='outer', label='batch')
ldata_merged.var_names_make_unique()
print(f"✅ Merge completo: {ldata_merged.n_obs} células, {ldata_merged.n_vars} genes")

# --- Guardar el objeto combinado ---
processed_h5ad = os.path.join(scvelo_outdir, 'processed_h5ad')
os.makedirs(processed_h5ad, exist_ok=True)
ldata_merged.write_h5ad(os.path.join(processed_h5ad, 'merged_all.h5ad'))
print("✅ Archivo guardado en processed_h5ad/merged_all.h5ad")

# --- Merge con el AnnData original (metadata, PCA, UMAP) ---
adata = scv.utils.merge(adata, ldata_merged)

# --- Preprocesamiento de scVelo ---
sc.pp.neighbors(adata, n_neighbors=15, use_rep='X_pca')  # calcular vecinos antes
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# --- Calcular RNA velocity ---
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

# --- Visualización ---
scv.pl.velocity_embedding(adata, basis='umap', frameon=False, save='embedding.pdf')
scv.pl.velocity_embedding_stream(adata, basis='umap', color=['celltype', 'condition'], save='embedding_stream.pdf')
scv.pl.proportions(adata, groupby='celltype_full')

# --- Guardar AnnData final ---
adata.write(os.path.join(scvelo_outdir, "adata_velocity.h5ad"))

EOF

