#!/usr/bin/env bash
set -euo pipefail

###################################################
## Batch execution of CoLoG pipeline
## Usage:
##   bash run_colog_all.sh /path/to/input_dir
## If no input directory is provided, the current
## working directory will be used as input.
## If your fastq files are in current directory
## (colog_env)$ ./run_colog_all.sh
##
## If your fastq files are in another directory
## (colog_env)$ ./run_colog_all.sh /path/to/inputs
###################################################

INPUT_DIR="${1:-.}"

echo "Input directory: ${INPUT_DIR}"

# Change to the input directory
cd "${INPUT_DIR}"

# Step 1: Split fastq files into sample folders
echo "=== Running mkchank.sh ==="
mkchank.sh
echo "=== mkchank.sh completed ==="

# Step 2: Run COLOG.sh on each sample directory sequentially
for SAMPLE_DIR in */; do
    SAMPLE_DIR="${SAMPLE_DIR%/}"  # remove trailing slash

    # Check the item is a directory
    if [ ! -d "${SAMPLE_DIR}" ]; then
        continue
    fi

    # Check if *_R1_*.fastq.gz exists in the folder (sample validation)
    if ! compgen -G "${SAMPLE_DIR}"/*_R1_*.fastq.gz > /dev/null; then
        echo ">>> Skipping ${SAMPLE_DIR} (no *_R1_*.fastq.gz found)"
        continue
    fi

    echo "=== Running COLOG.sh on: ${SAMPLE_DIR} ==="
    COLOG.sh --threads 16 "${SAMPLE_DIR}"
    echo "=== Finished: ${SAMPLE_DIR} ==="
    echo
done

echo "All samples completed successfully."
