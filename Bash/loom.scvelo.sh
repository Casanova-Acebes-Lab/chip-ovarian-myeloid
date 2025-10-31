#!/bin/sh
#SBATCH -J scvelo_load
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=128G
#SBATCH -t 08:00:00
#SBATCH -o /home/dcaceres/Tet2sc/Bash/loom.read.out
#SBATCH -e /home/dcaceres/Tet2sc/Bash/loom.read.err

# ============================================================
# CONFIGURACIÓN DEL ENTORNO
# ============================================================
source /home/dcaceres/miniconda3/etc/profile.d/conda.sh
conda activate scvelo_env

# ============================================================
# RUTAS DE ENTRADA / SALIDA
# ============================================================
LOOM_DIR=/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files
OUT_DIR=$LOOM_DIR/processed_h5ad
mkdir -p "$OUT_DIR"

# ============================================================
# BLOQUE PYTHON
# ============================================================
python <<PYCODE
import os
import glob
import scvelo as scv

loom_dir = "$LOOM_DIR"
out_dir = "$OUT_DIR"

# Tomar todos los archivos .loom de la carpeta
loom_files = sorted(glob.glob(os.path.join(loom_dir, "*.loom")))

for idx, loom_path in enumerate(loom_files, start=1):
    print("========================================")
    print(f"Procesando archivo: {loom_path}")
    print(f"ID de muestra: {idx}")
    print(f"Fecha: {os.popen('date').read().strip()}")
    print(f"Salida final: {out_dir}")
    print("========================================")

    # Leer archivo loom (scVelo 0.3.3)
    ldata = scv.read(loom_path, cache=True)

    # Renombrar barcodes
    barcodes = [bc.split(':')[1] for bc in ldata.obs.index]
    barcodes = [bc[:-1] + f"_{idx}" for bc in barcodes]
    ldata.obs.index = barcodes

    # Guardar en subdirectorio como .h5ad
    out_path = os.path.join(out_dir, os.path.splitext(os.path.basename(loom_path))[0] + f"_proc_{idx}.h5ad")
    ldata.write(out_path)
    print(f"[{idx}] Guardado en {out_path}")

PYCODE





