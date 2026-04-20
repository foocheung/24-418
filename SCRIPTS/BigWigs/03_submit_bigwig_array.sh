#!/usr/bin/bash
set -euo pipefail

WORKROOT="/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe"
MANIFEST="${WORKROOT}/bigwig_manifest.tsv"

NJOBS=$(($(wc -l < "${MANIFEST}") - 1))

if [[ "${NJOBS}" -le 0 ]]; then
    echo "No jobs found in manifest: ${MANIFEST}"
    exit 1
fi

echo "Submitting ${NJOBS} array jobs"

sbatch \
  --job-name=bigwig_ct \
  --array=1-"${NJOBS}"%32 \
  --cpus-per-task=8 \
  --mem=32G \
  --time=24:00:00 \
  --output="${WORKROOT}/slurm_bigwig_%A_%a.out" \
  --error="${WORKROOT}/slurm_bigwig_%A_%a.err" \
  /data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/04_sbatch_wrapper.sh
