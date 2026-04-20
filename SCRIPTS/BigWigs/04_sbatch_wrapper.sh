#!/usr/bin/bash
#SBATCH --export=ALL

set -euo pipefail

WORKROOT="/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe"
MANIFEST="${WORKROOT}/bigwig_manifest.tsv"

# Ensure your environment is active
#source ~/.bashrc
module load samtools

# Optional but safer (forces correct path)
export PATH="$HOME/.local/bin:$PATH"

export THREADS="${SLURM_CPUS_PER_TASK:-8}"
export BINSIZE=25

bash "${WORKROOT}/02_bigwig_worker.sh" "${MANIFEST}"

