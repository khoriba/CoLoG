#!/usr/bin/env bash
set -euo pipefail

##############################################
# Batch execution of CoLoG pipeline Ver.2
# Usage:
# Usage:
#   bash run_colog_all.sh <input_dir> [primer.bed]
##############################################

INPUT_DIR="$1"
PRIMER_BED="${2:-}"

echo "Input directory: ${INPUT_DIR}"
echo "Primer BED: ${PRIMER_BED:-<none>}"

cd "${INPUT_DIR}"

echo "=== Running mkchank.sh ==="
mkchank.sh
echo "=== mkchank.sh completed ==="

for SAMPLE_DIR in */; do
    SAMPLE_DIR="${SAMPLE_DIR%/}"

    # Validate input folder
    if [ ! -d "${SAMPLE_DIR}" ]; then
        continue
    fi
    if ! compgen -G "${SAMPLE_DIR}"/*_R1_*.fastq.gz > /dev/null; then
        echo ">>> Skipping ${SAMPLE_DIR} (no *_R1_*.fastq.gz found)"
        continue
    fi

    echo "=== Running COLOG.sh on: ${SAMPLE_DIR} ==="

    if [ -n "${PRIMER_BED}" ] && [ -f "${PRIMER_BED}" ]; then
        COLOG.sh --primer "${PRIMER_BED}" --threads 16 "${SAMPLE_DIR}"
    else
        COLOG.sh --threads 16 "${SAMPLE_DIR}"
    fi

    echo "=== Finished: ${SAMPLE_DIR} ==="
    echo
done

echo "All samples completed successfully."
