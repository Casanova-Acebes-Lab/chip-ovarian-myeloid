import pandas as pd

# Leer archivo TSV
config = pd.read_csv("config.tsv", sep="\t")

# Convertir a diccionario clave → valor
config_dict = dict(zip(config['key'], config['value']))

# Acceder a las variables
rsdir = config_dict["rsdir"]
datadir = config_dict["datadir"]
outdir = config_dict["outdir"]

print(rsdir, datadir, outdir)

# Read data


import scanpy as sc
import anndata
from scipy import io
import numpy as np
import os
import subprocess


# --- Directorios ---
scvelo_dir = os.path.join(rsdir, "objects/scVelo")
scvelo_outdir = os.path.join(outdir, "scVelo")
os.makedirs(scvelo_outdir, exist_ok=True)

# --- Cargar matriz ---
X = io.mmread(os.path.join(scvelo_dir, "counts.mtx"))
adata = anndata.AnnData(X=X.transpose().tocsr())

# --- Cargar metadatos ---
cell_meta = pd.read_csv(os.path.join(scvelo_dir, "metadata.csv"))
cell_meta.columns = cell_meta.columns.str.strip()
adata.obs = cell_meta
adata.obs.index = adata.obs['barcode']

# --- Cargar nombres de genes ---
with open(os.path.join(scvelo_dir, "gene_names.csv"), 'r') as f:
    adata.var.index = f.read().splitlines()

# --- Cargar PCA/UMAP ---
pca = pd.read_csv(os.path.join(scvelo_dir, "pca.csv"))
pca.index = adata.obs.index
adata.obsm['X_pca'] = pca.to_numpy()
adata.obsm['X_umap'] = np.vstack((adata.obs['UMAP_1'], adata.obs['UMAP_2'])).T


# --- Paleta de R como diccionario ---
nora_colors = {
    "Ly6c|Ms4ac Monocytes": "#FF3B30",
    "Trem1|Ptgs2|Plaur|Celc4e Mac": "#EE7942",
    "Mrc1|C1qc|Cbr2|Gas6 Mac": "#FFD92F",
    "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac": "#4DAF4A",
    "Npr2|Actn1 Mac": "#A6D854",
    "Cd8 T cells": "#00BFC4",
    "Cd4 T cells": "#FF69B4",
    "Neutrophils": "#4876FF",
    "DCs": "#87CEEB",
    "NK": "#AB82FF",
    "Tgd": "#D9B3FF",
    "Mastocytes": "#3CB371",
    "B cells": "#DC143C",
    "Cd8 Effector": "#A52A2A",
    "Marco+ Mac": "#1E3A8A",
    "Slamf Monocytes": "#66B2FF",
    "Mmp9|Ctsk Mac": "#00723F",
    "IFN Mac": "#C080FF",
    "Fn1|Vegfa Mac": "#FFA500",
    "Ciita|Siglec Mac": "#1E90FF",
    "Ciita|Ccl12 Mac": "#4682B4",
    "Activated B cells": "#FF1493"
}

# --- Columna de interés ---
col = 'Cluster'

# --- Filtrar solo categorías presentes en tu subset ---
cats_present = [cat for cat in nora_colors.keys() if cat in adata.obs[col].unique()]

# Convertir la columna a categórica usando solo las categorías presentes
adata.obs[col] = pd.Categorical(
    adata.obs[col],
    categories=cats_present,
    ordered=True
)

# Crear lista de colores para estas categorías
palette = [nora_colors[cat] for cat in cats_present]

# Configurar directorio de salida
scvelo_outdir = os.path.join(outdir, "scVelo")
os.makedirs(scvelo_outdir, exist_ok=True)
sc.settings.figdir = scvelo_outdir

# Graficar UMAP con colores exactos del subset
sc.pl.umap(
    adata,
    color=col,
    palette=palette,
    frameon=False,
    save='_Cell_type2.png'
)



# --- Guardar AnnData ---
adata.write(os.path.join(scvelo_outdir, "my_data.h5ad"))





# Step 0: Constructing spliced and unspliced counts matrices
# See Bash/Velocyto.sh for this step

# Step 1: Load data

import scvelo as scv
import scanpy as sc
import cellrank as cr
import numpy as np
import pandas as pd
import anndata as ad

scv.settings.verbosity = 3
scv.settings.set_figure_params('scvelo', facecolor='white', dpi=100, frameon=False)
cr.settings.verbosity = 2


adata = sc.read_h5ad(f'{outdir}/scVelo/my_data.h5ad')




ldata1 = sc.read_loom('/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedN-KO2.loom')
ldata2 = sc.read_loom('/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedN-KO3.loom')
ldata3 = sc.read_loom('/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedP-KO2.loom')
ldata4 = sc.read_loom('/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/smpR065_DsRedP-KO3.loom')
ldata5 = sc.read_loom('/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT1-DsRedn.loom')
ldata6 = sc.read_loom('/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT1-DsRedp.loom')
ldata7 = sc.read_loom('/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT2-DsRedn.loom')
ldata8 = sc.read_loom('/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/WT2-DsRedp.loom')


# rename barcodes in order to merge:
barcodes = [bc.split(':')[-1] for bc in ldata1.obs.index]
barcodes = [bc + "_1" for bc in barcodes]
ldata1.obs.index = barcodes


barcodes = [bc.split(':')[-1] for bc in ldata2.obs.index]
barcodes = [bc + "_1" for bc in barcodes]
ldata2.obs.index = barcodes


barcodes = [bc.split(':')[-1] for bc in ldata3.obs.index]
barcodes = [bc + "_1" for bc in barcodes]
ldata3.obs.index = barcodes

barcodes = [bc.split(':')[-1] for bc in ldata4.obs.index]
barcodes = [bc + "_1" for bc in barcodes]
ldata4.obs.index = barcodes

barcodes = [bc.split(':')[-1] for bc in ldata5.obs.index]
barcodes = [bc + "_1" for bc in barcodes]
ldata5.obs.index = barcodes


 barcodes = [bc.split(':')[-1] for bc in ldata6.obs.index]
barcodes = [bc + "_1" for bc in barcodes]
ldata6.obs.index = barcodes

barcodes = [bc.split(':')[-1] for bc in ldata7.obs.index]
barcodes = [bc + "_1" for bc in barcodes]
ldata7.obs.index = barcodes


barcodes = [bc.split(':')[-1] for bc in ldata8.obs.index]
barcodes = [bc + "_1" for bc in barcodes]
ldata8.obs.index = barcodes


# make variable names unique
ldata1.var_names_make_unique()
ldata2.var_names_make_unique()
ldata3.var_names_make_unique()
ldata4.var_names_make_unique()
ldata5.var_names_make_unique()
ldata6.var_names_make_unique()
ldata7.var_names_make_unique()
ldata8.var_names_make_unique()




# concatenate the three loom
ldata = ldata1.concatenate([ldata2, ldata3,ldata4,ldata5, ldata6, ldata7, ldata8])



# Rutas a tus archivos loom
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



import anndata as ad

# 🔧 1️⃣ Alinear genes
common_genes = ldata1.var_names
for l in [ldata2, ldata3, ldata4, ldata5, ldata6, ldata7, ldata8]:
    common_genes = common_genes.intersection(l.var_names)
print(f"✅ Genes comunes: {len(common_genes)}")

# 🔧 2️⃣ Filtrar genes y convertir a copias reales
ldata_list = []
for i, ld in enumerate([ldata1, ldata2, ldata3, ldata4, ldata5, ldata6, ldata7, ldata8], start=1):
    print(f"Procesando ldata{i} ...")
    ld = ld[:, common_genes].copy()  # aseguramos vista -> copia
    # 🔹 eliminar 'matrix' de forma segura (sin activar validación)
    if hasattr(ld, "_layers") and "matrix" in ld._layers:
        print(f"⚠️ Eliminando capa 'matrix' (inconsistente) en ldata{i}")
        del ld._layers["matrix"]
    ldata_list.append(ld)

# 🔧 3️⃣ Asignar etiquetas de batch
batch_keys = ["KO2n", "KO3n", "KO2p", "KO3p", "WT1n", "WT1p", "WT2n", "WT2p"]
for key, adata in zip(batch_keys, ldata_list):
    adata.obs["batch"] = key

# 🔧 4️⃣ Concatenar
print("🧩 Concatenando todos los AnnData...")
ldata_merged = ad.concat(
    ldata_list,
    join="inner",
    label="batch",
    keys=batch_keys,
    index_unique=None
)

print(f"✅ Merge completo: {ldata_merged.n_obs} células, {ldata_merged.n_vars} genes")

# 🔧 5️⃣ Guardar resultado
out_path = "/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/processed_h5ad/merged_all.h5ad"
ldata_merged.write_h5ad(out_path)
print(f"✅ Archivo guardado en: {out_path}")



print(adata)





import scanpy as sc
import scvelo as scv

# Ruta al archivo que generaste
adata = sc.read_h5ad("/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/processed_h5ad/merged_with_layers.h5ad")

# Revisar contenido
print(adata)





# plot umap to check
sc.pl.umap(adata, 
frameon=False, legend_loc='on data',      color=col,
 title='', save='_celltypes.3.pdf')



sc.pl.umap(
    adata,
    color=col,
    palette=palette,
    frameon=False,
    save='_Cell_type2.png'
)


scv.pl.proportions(adata, groupby='celltype_full')
















import os
import scanpy as sc
import scvelo as scv
import anndata as ad

# --- Rutas de tus loom ---
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

# --- Cargar todos los loom ---
ldata_list = []
for path in loom_paths:
    ld = scv.read(path, cache=True)
    # limpiar capa 'matrix' inconsistente
    if "matrix" in ld.layers and ld.layers["matrix"].shape != ld.X.shape:
        del ld.layers["matrix"]
    # hacer nombres de genes únicos
    ld.var_names_make_unique()
    ldata_list.append(ld)

# --- Alinear genes: solo los comunes entre todos los loom ---
common_genes = set(ldata_list[0].var_names)
for ld in ldata_list[1:]:
    common_genes &= set(ld.var_names)
common_genes = sorted(common_genes)
print(f"✅ Genes comunes: {len(common_genes)}")

# Filtrar cada AnnData por los genes comunes y hacer barcodes únicos
for i, ld in enumerate(ldata_list):
    ld = ld[:, common_genes].copy()
    ld.obs_names = [f"{batch_keys[i]}_{c}" for c in ld.obs_names]  # barcodes únicos
    ld.obs["batch"] = batch_keys[i]
    ldata_list[i] = ld

# --- Concatenar todos los loom ---
ldata_merged = ad.concat(
    ldata_list,
    join='inner',  # mantener solo genes comunes
    label='batch',
    keys=batch_keys,
    index_unique=None
)
ldata_merged.obs_names_make_unique()
print(f"✅ Loom merge completo: {ldata_merged.n_obs} células, {ldata_merged.n_vars} genes")

# --- Sincronizar con tu adata original ---
common_cells = adata.obs_names.intersection(ldata_merged.obs_names)
print(f"✅ Celdas comunes con adata original: {len(common_cells)}")

adata_sub = adata[common_cells].copy()
ldata_sub = ldata_merged[common_cells].copy()

# --- Transferir layers spliced/unspliced/ambiguous ---
adata_sub.layers["spliced"] = ldata_sub.layers["spliced"]
adata_sub.layers["unspliced"] = ldata_sub.layers["unspliced"]
if "ambiguous" in ldata_sub.layers:
    adata_sub.layers["ambiguous"] = ldata_sub.layers["ambiguous"]

# Reemplazar adata original con el subdataset listo para velocity
adata = adata_sub
adata.obs_names_make_unique()

print(f"✅ Merge final listo: {adata.n_obs} células, {adata.n_vars} genes")





adata = sc.read_h5ad("/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/processed_h5ad/merged_clean.h5ad")
print(adata)





# --- Recalcular PCA y UMAP ---
sc.pp.normalize_total(ldata_merged)
sc.pp.log1p(ldata_merged)
sc.pp.highly_variable_genes(ldata_merged, n_top_genes=2000)
sc.pp.scale(ldata_merged)
sc.tl.pca(ldata_merged)
sc.pp.neighbors(ldata_merged)
sc.tl.umap(ldata_merged)

# --- Guardar resultado final ---
out_path = "/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files/processed_h5ad/merged_clean.h5ad"
ldata_merged.write(out_path)
print(f"✅ Archivo guardado en: {out_path}")

# --- Plot UMAP de prueba ---
col = "batch"
sc.pl.umap(adata, color=col, frameon=False, save="_batch_test.png")

scv.pl.proportions(adata, groupby='celltype_full',save="_proportion.png")