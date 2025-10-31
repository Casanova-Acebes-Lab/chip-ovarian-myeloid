#!/bin/sh
#SBATCH --job-name=velocyto_run_multi
#SBATCH --array=0-7
#SBATCH --cpus-per-task=8
#SBATCH --mem=96G
#SBATCH --time=24:00:00
#SBATCH --output=/home/dcaceres/Tet2sc/Bash/velocyto_run_multi_%A_%a.out
#SBATCH --error=/home/dcaceres/Tet2sc/Bash/velocyto_run_multi_%A_%a.err

# ============================================================
# VELOCYTO RUN SCRIPT PARA CELLRANGER MULTI
# Optimizado para 3 librerías y 8 muestras
# ============================================================

echo "------------------------------------------------------------"
echo "🚀 INICIO DEL PROCESO VELOCYTO - $(date)"
echo "------------------------------------------------------------"

# --- Activar entorno con scvelo y velocyto ---
source /home/dcaceres/miniconda3/etc/profile.d/conda.sh
conda activate velocyto_env

# --- Rutas principales ---
DATADIR="/storage/scratch01/groups/ci/Trem1_GEICAM/Sarai/Single.Cell.Tet2/outs"
OUTDIR="/home/dcaceres/Tet2sc/Results/plots/scVelo/velocyto_loom_files"
mkdir -p "$OUTDIR"

REPEATS="/storage/scratch01/groups/ci/Trem1_GEICAM/Sarai/Annonation_Files/rmsk_velocyto.txt"
TRANSCRIPTOME="/storage/scratch01/groups/ci/Trem1_GEICAM/Sarai/Annonation_Files/gencode.vM25.annotation.gtf"

# --- Archivo de resumen global ---
SUMMARY_FILE="$OUTDIR/velocyto_summary.log"
touch "$SUMMARY_FILE"

# --- Lista de las 8 muestras ---
samples=(
  "out.smR039a.9/outs/per_sample_outs/WT1-DsRedn"
  "out.smR039a.9/outs/per_sample_outs/WT1-DsRedp"
  "out.smR039b.9/outs/per_sample_outs/WT2-DsRedn"
  "out.smR039b.9/outs/per_sample_outs/WT2-DsRedp"
  "out.smR065.9/outs/per_sample_outs/smpR065_DsRedN-KO2"
  "out.smR065.9/outs/per_sample_outs/smpR065_DsRedN-KO3"
  "out.smR065.9/outs/per_sample_outs/smpR065_DsRedP-KO2"
  "out.smR065.9/outs/per_sample_outs/smpR065_DsRedP-KO3"
)

SAMPLE=${samples[$SLURM_ARRAY_TASK_ID]}
SAMPLE_NAME=$(basename $SAMPLE)
COUNT_DIR="$DATADIR/$SAMPLE/count"
BAM="$COUNT_DIR/sample_alignments.bam"
BARCODES="$COUNT_DIR/sample_filtered_feature_bc_matrix/barcodes.tsv.gz"
OUTPUT_LOOM="$OUTDIR/${SAMPLE_NAME}.loom"

echo "=== Procesando muestra: $SAMPLE_NAME ==="
echo "Directorio de entrada: $COUNT_DIR"
echo "Nodo: $SLURMD_NODENAME | JobID: $SLURM_JOB_ID | TaskID: $SLURM_ARRAY_TASK_ID"
echo "------------------------------------------------------------"

# --- Verificar archivos necesarios ---
if [ ! -f "$BAM" ]; then
    echo "❌ ERROR: No se encontró BAM: $BAM"
    echo -e "$(date)\t$SAMPLE_NAME\tERROR\tBAM no encontrado" >> "$SUMMARY_FILE"
    exit 1
fi

if [ ! -f "$BARCODES" ]; then
    echo "❌ ERROR: No se encontró barcodes.tsv.gz: $BARCODES"
    echo -e "$(date)\t$SAMPLE_NAME\tERROR\tBarcodes no encontrado" >> "$SUMMARY_FILE"
    exit 1
fi

# --- Ejecutar velocyto ---
echo "🧬 Ejecutando Velocyto..."
velocyto run \
    -b "$BARCODES" \
    -e "$SAMPLE_NAME" \
    -o "$OUTDIR" \
    -m "$REPEATS" \
    "$BAM" \
    "$TRANSCRIPTOME"

# --- Comprobación de salida ---
if [ -f "$OUTPUT_LOOM" ]; then
    echo "✅ ÉXITO: Loom creado correctamente: $OUTPUT_LOOM"
    STATUS="OK"
else
    echo "⚠ ADVERTENCIA: No se encontró el archivo Loom esperado: $OUTPUT_LOOM"
    STATUS="ERROR"
fi

# --- Registro en resumen global ---
echo -e "$(date)\t$SAMPLE_NAME\t$STATUS" >> "$SUMMARY_FILE"

# --- Mensaje final ---
echo "------------------------------------------------------------"
if [ "$STATUS" = "OK" ]; then
    echo "🎉 Proceso finalizado correctamente para $SAMPLE_NAME"
else
    echo "❌ Proceso finalizado con errores para $SAMPLE_NAME"
fi
echo "Finalizado en: $(date)"
echo "------------------------------------------------------------"
