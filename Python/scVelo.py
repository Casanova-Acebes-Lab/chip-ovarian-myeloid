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
import scvelo as scv


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
    # --- MONOCITOS & TAMs (ROJOS) ---
  "Ly6cHi Monocytes"     = "#FF0000",  # Rojo inflamatorio puro
  "Ly6cLo Monocytes"     = "#FF6A6A",  # Rojo salmón transición
  "Early IFN|MHCII-TAMs"       = "#B22222",  # Rojo vino (pre-TAM IFN)
  "Trem1|Ptgs2|Plaur|Celc4e Mac" = "#EE7942",   # naranja vivo (inflam. activados)
  "Mrc1|C1qc|Cbr2|Gas6 Mac"      = "#FFD92F",   # amarillo brillante (TAM residentes)
  "Arg1|Spp1|Mmp12|Mmp19|Il1a Mac" = "#4DAF4A", # verde TAM reparadores
  "Npr2|Actn1 Mac"       = "#A6D854",   # verde lima (TAM estructurales)
  "Marco+ Mac"           = "#1E3A8A",   # azul intermedio (scavenger)
  "Mmp9|Ctsk Mac"        = "#00723F",   # verde botella (remodelado matriz)
  "IFN Mac"              = "#C080FF",   # morado claro (TAM interferón maduros)
  "Fn1|Vegfa Mac"        = "#FFA500",   # naranja angiogénico
  "MHCII|Siglec Mac"     = "#1E90FF",   # azul cobalto (antigen presenting)
  "MHCII|Ccl12 Mac"      = "#4682B4",   # azul acero

  # --- CÉLULAS LINFOIDES ---
  "Cd8 T cells"          = "#00BFC4",   # turquesa
  "Cd4 T cells"          = "#FF69B4",   # rosa fuerte
  "Cd8 Effector"         = "#A52A2A",   # rojo ladrillo
  "Tgd"                  = "#D9B3FF",   # lila pastel
  "NK"                   = "#AB82FF",   # violeta oscuro
  "Activated B cells"    = "#FF1493",   # fucsia intenso
  "B cells"              = "#DC143C",   # rojo intenso

  # --- INNATAS ---
  "Neutrophils"          = "#4876FF",   # azul fuerte
  "DCs"                  = "#87CEEB",   # azul claro
  "Mastocytes"           = "#3CB371"    # verde saturado
}

# --- Columna de interés ---
col = 'Clustering.Round2'

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

# Step 1: merge loom files and integrate with existing adata
# See Bash/merge.loom.sh for this step
import scvelo as scv
import scanpy as sc

scv.set_figure_params(style='white', dpi=150)

# ==========
# LOAD DATA
# ==========
adata = sc.read_h5ad(f"{outdir}/scVelo/velocyto_loom_files/processed_h5ad/adata_with_velocity_layers_subsetted.h5ad")

print("Genes totales antes del filtrado:", adata.n_vars)

# ==========
# 1) FILTRADO + NORMALIZACIÓN (pipeline recomendado)
# ==========
scv.pp.filter_and_normalize(
    adata,
    min_shared_counts=20,   # ← robusto, mantiene genes expresados en suficientes células
    n_top_genes=2000        # ← selecciona las más informativas
)

# ==========
# 2) PCA + Vecinos
# ==========
sc.pp.pca(adata, n_comps=30)
sc.pp.neighbors(adata, n_neighbors=30, n_pcs=30)

# ==========
# 3) MOMENTS
# ==========
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)

# ==========
# 4) VELOCITY + GRAPH
# ==========
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)

print("Genes tras filtrado:", adata.n_vars)

# ==========
# 5) VISUALIZACIÓN GLOBAL
# ==========
scv.pl.velocity_embedding_stream(
    adata, basis='umap', color='Cluster',
    legend_loc='right margin', save='velocity_stream_all.pdf'
)

scv.pl.velocity_embedding_grid(
    adata, basis='umap', color='Cluster',
    arrow_size=2.0, arrow_length=4.0,
    legend_loc='right margin', save='velocity_grid_all.pdf'
)

# ==================================
# WT

# ==================================
adata_WT = adata[adata.obs['group'] == 'WT'].copy()
scv.tl.velocity_embedding(adata_WT, basis='umap')

scv.pl.velocity_embedding_stream(
    adata_WT, basis='umap', color='Cluster',
    legend_loc='right margin', save='WT_velocity_stream.pdf'
)

scv.pl.velocity_embedding_grid(
    adata_WT, basis='umap', color='Cluster',
    arrow_size=2.0, arrow_length=4.0,
    legend_loc='right margin', save='WT_velocity_grid.pdf'
)

# ==================================
# KO
# ==================================
adata_KO = adata[adata.obs['group'] == 'KO'].copy()
scv.tl.velocity_embedding(adata_KO, basis='umap')

scv.pl.velocity_embedding_stream(
    adata_KO, basis='umap', color='Cluster',
    legend_loc='right margin', save='KO_velocity_stream.pdf'
)

scv.pl.velocity_embedding_grid(
    adata_KO, basis='umap', color='Cluster',
    arrow_size=2.0, arrow_length=4.0,
    legend_loc='right margin', save='KO_velocity_grid.pdf'
)






scv.set_figure_params(
    style='white',
    dpi=150,
    frameon=True,             # ← importante para fondo blanco
    facecolor='white'         # ← fuerza fondo blanco
)



scv.pl.velocity(
    adata_WT,
    var_names=['Ifi206'],
    color='Cluster',
    dpi=150,
    save='Ifi206_WT_velocity_WHITE.png'
)




### Part 3: Downstream analysis


import pandas as pd

df = pd.DataFrame(adata_WT.uns['rank_velocity_genes']['names'])

# Mostrar tabla completa por pantalla
print(df)


scv.tl.velocity_confidence(adata_WT)
keys = 'velocity_length', 'velocity_confidence'
scv.pl.scatter(adata_WT, c=keys, cmap='coolwarm', perc=[5, 95],
save='WT_velocity.lenght.confidence.pdf')


scv.pl.velocity_graph(adata_WT, threshold=.1, color='Cluster',
save='WT_velocity.graph.pdf')




scv.tl.velocity_pseudotime(adata_WT)
scv.pl.scatter(adata_WT, color='velocity_pseudotime', cmap='gnuplot',
save='WT_velocity.pseudotime.png')



scv.tl.velocity_pseudotime(adata_KO)
scv.pl.scatter(adata_KO, color='velocity_pseudotime', cmap='gnuplot',
save='KO_velocity.pseudotime.png')




import scvelo as scv

# 1) Completely remove old graph info
for key in ["neighbors", "paga"]:
    if key in adata_WT.uns:
        del adata_WT.uns[key]

for key in ["distances", "connectivities"]:
    if key in adata_WT.obsp:
        del adata_WT.obsp[key]

# 2) Recompute PCA if missing (won't overwrite existing X_umap or clustering)
scv.pp.filter_and_normalize(adata_WT)
scv.pp.pca(adata_WT)

# 3) ***REBUILD NEIGHBOR GRAPH FOR VELOCITY***
scv.pp.moments(adata_WT, n_pcs=30, n_neighbors=30)

# 4) Recompute velocities
scv.tl.velocity(adata_WT, mode="stochastic")
scv.tl.velocity_graph(adata_WT)

# 5) Now PAGA will run successfully
import scanpy as sc
sc.tl.paga(adata_WT, groups="Cluster")






scv.pl.paga(adata_WT, basis='umap', size=50, alpha=.1,
             min_edge_width=2, node_size_scale=1.5,
             save='wt_velocity.graph.png')


scv.pl.paga(adata_WT, basis='umap',
             size=50, alpha=0.1,
             min_edge_width=2,
             node_size_scale=1.5,
             save='wt_paga.pdf')







import scvelo as scv

# 1) Completely remove old graph info
for key in ["neighbors", "paga"]:
    if key in adata_KO.uns:
        del adata_KO.uns[key]

for key in ["distances", "connectivities"]:
    if key in adata_KO.obsp:
        del adata_KO.obsp[key]

# 2) Recompute PCA if missing (won't overwrite existing X_umap or clustering)
scv.pp.filter_and_normalize(adata_KO)
scv.pp.pca(adata_KO)

# 3) ***REBUILD NEIGHBOR GRAPH FOR VELOCITY***
scv.pp.moments(adata_KO, n_pcs=30, n_neighbors=30)

# 4) Recompute velocities
scv.tl.velocity(adata_KO, mode="stochastic")
scv.tl.velocity_graph(adata_KO)

# 5) Now PAGA will run successfully
import scanpy as sc
sc.tl.paga(adata_KO, groups="Cluster")




scv.pl.paga(adata_KO, basis='umap',
             size=50, alpha=0.1,
             min_edge_width=2,
             node_size_scale=1.5,
             save='wt_paga.KO.pdf')





# Part 4: Analyzing a specific cell population



df = adata_WT.var
df = df[(df['fit_likelihood'] > .1) & df['velocity_genes'] == True]

kwargs = dict(xscale='log', fontsize=16)
with scv.GridSpec(ncols=3) as pl:
    pl.hist(df['fit_alpha'], xlabel='transcription rate', **kwargs)
    pl.hist(df['fit_beta'] * df['fit_scaling'], xlabel='splicing rate', xticks=[.1, .4, 1], **kwargs)
    pl.hist(df['fit_gamma'], xlabel='degradation rate', xticks=[.1, .4, 1], **kwargs)

scv.get_df(adata_WT, 'fit*', dropna=True).head()


scv.tl.latent_time(adata_WT)
scv.pl.scatter(adata_WT, color='latent_time', color_map='gnuplot', size=80,
save='WT_velocity.latent.time.pdf')


import scvelo as scv
import pandas as pd

# ----------------------------
# 1️⃣ Selección de clusters
# ----------------------------
cluster_A = "Slamf Monocytes"
cluster_B = "IFN Mac"

subset = adata_WT[adata_WT.obs["Cluster"].isin([cluster_A, cluster_B])].copy()

# ----------------------------
# 2️⃣ Rankear genes dinámicos (usar min_corr bajo para asegurar resultados)
# ----------------------------
scv.tl.rank_velocity_genes(subset, groupby="Cluster", min_corr=0.05)

# ----------------------------
# 3️⃣ Extraer top genes y exportar
# ----------------------------
top_genes = subset.uns['rank_velocity_genes']['names'][:30]  # top 30 genes
pd.DataFrame(top_genes, columns=['gene']).to_csv("WT_top_genes_Slamf_IFN.csv", index=False)

# ----------------------------
# 4️⃣ Visualizar los top genes en UMAP
# ----------------------------
for gene in top_genes[:10]:  # top 10 para visualizar
    scv.pl.scatter(subset, color=gene, basis='umap', show=True)





gene = "Trem1"
import numpy as np

if gene not in adata_WT.var_names:
    print("El gen NO está en adata_WT.var_names")
else:
    idx = adata_WT.var_names.get_loc(gene)

    total_X = np.sum(adata_WT.X[:, idx].data)
    total_spliced = np.sum(adata_WT.layers["spliced"][:, idx].data) if "spliced" in adata_WT.layers else None
    total_unspliced = np.sum(adata_WT.layers["unspliced"][:, idx].data) if "unspliced" in adata_WT.layers else None

    print(f"\nGen: {gene}")
    print("  Total X:", total_X)
    print("  Total spliced:", total_spliced)
    print("  Total unspliced:", total_unspliced)






gene = "Trem1"
if gene in adata_WT.var_names:
    idx = adata_WT.var_names.get_loc(gene)
    total_X = float(adata_WT.X[:, idx].sum())
    total_spliced = float(adata_WT.layers["spliced"][:, idx].sum())
    total_unspliced = float(adata_WT.layers["unspliced"][:, idx].sum())
    print("Total counts in X:", total_X)
    print("Total spliced:", total_spliced)
    print("Total unspliced:", total_unspliced)
else:
    print("El gen NO está en adata_WT.var_names")






import numpy as np

for layer in ['spliced', 'unspliced', 'ambiguous']:
    if layer in adata.layers:
        X = adata.layers[layer]
        # toma una muestra pequeña de 1000 células (o menos si hay menos)
        n = min(1000, X.shape[0])
        rows = np.random.choice(X.shape[0], n, replace=False)

        # suma solo esas filas
        try:
            partial_sum = X[rows, :].sum()
        except Exception:
            partial_sum = np.array(X[rows, :].sum()).sum()

        mean_per_cell = partial_sum / n
        print(f"{layer}: shape={X.shape}, mean counts/cell≈{mean_per_cell:.2f}")
    else:
        print(f"{layer}: not found")




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