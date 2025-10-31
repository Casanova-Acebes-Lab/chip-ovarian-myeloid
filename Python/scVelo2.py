
import os
import scanpy as sc
import anndata
import pandas as pd
import numpy as np
import scvelo as scv

# ============================================================
# 1️⃣ Cargar adata base (desde Seurat convertido)
# ============================================================
adata_path = "/home/dcaceres/Tet2sc/Results/plots/scVelo/my_data.h5ad"
adata = sc.read_h5ad(adata_path)
print("adata cargado:", adata)

# ============================================================
# 2️⃣ Cargar archivos loom generados por velocyto
# ============================================================
loom_dir = "/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files"
loom_files = [
    "smpR065_DsRedN-KO2.loom",
    "smpR065_DsRedN-KO3.loom",
    "smpR065_DsRedP-KO2.loom",
    "smpR065_DsRedP-KO3.loom",
    "WT1-DsRedn.loom",
    "WT1-DsRedp.loom",
    "WT2-DsRedn.loom",
    "WT2-DsRedp.loom"
]

ldata = []
for i, lf in enumerate(loom_files, start=1):
    path = os.path.join(loom_dir, lf)
    print(f"🔹 Cargando {path}")
    ld = scv.read_loom(path, X_name='spliced', obs_names='CellID', var_names='Gene')
    ld.var_names_make_unique()
    ld.obs['batch'] = f"sample_{i}"
    # ajustar barcodes para que coincidan con tu adata
    ld.obs.index = [bc.split(":")[-1] + f"_{i}" for bc in ld.obs.index]
    ldata.append(ld)

# ============================================================
# 3️⃣ Concatenar los loom
# ============================================================
ldata_merged = ldata[0].concatenate(ldata[1:], batch_key='batch', index_unique=None)
print("ldata_merged:", ldata_merged)

# ============================================================
# 4️⃣ Merge con adata base
# ============================================================
adata = scv.utils.merge(adata, ldata_merged)
print("adata final:", adata)
print("Capas disponibles:", adata.layers.keys())

# ============================================================
# 5️⃣ Preprocesamiento scVelo
# ============================================================
scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)

scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata)

# ============================================================
# 6️⃣ Calcular RNA velocity
# ============================================================
scv.tl.velocity(adata, mode='stochastic')  # o 'dynamical' si prefieres el modelo 2020
scv.tl.velocity_graph(adata)

# ============================================================
# 7️⃣ Visualización
# ============================================================
# Proporciones de spliced/unspliced por cluster
scv.pl.proportions(adata, groupby='celltype_full', save="_proportions.png")

# UMAP con RNA velocity en stream
scv.pl.velocity_embedding_stream(
    adata,
    basis='umap',
    color=['celltype', 'condition'],
    save="_embedding_stream.png"
)

# ============================================================
# 8️⃣ Guardar adata final
# ============================================================
output_file = os.path.join(loom_dir, "processed_h5ad/adata_scvelo_merged.h5ad")
os.makedirs(os.path.dirname(output_file), exist_ok=True)
adata.write(output_file)
print(f"✅ Archivo guardado en {output_file}")