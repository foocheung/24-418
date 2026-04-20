#!/usr/bin/bash
set -euo pipefail

MANIFEST="${1:-/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/bigwig_manifest.tsv}"
THREADS="${THREADS:-8}"
BINSIZE="${BINSIZE:-25}"

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID is not set"
    exit 1
fi

LINE_NUM=$((SLURM_ARRAY_TASK_ID + 1))

ROW=$(sed -n "${LINE_NUM}p" "${MANIFEST}")
if [[ -z "${ROW}" ]]; then
    echo "No manifest row for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

IFS=$'\t' read -r GEX LANE CELLTYPE BCFILE BAM OUTDIR TMPDIR LOGDIR BW SORTBAM <<< "${ROW}"

mkdir -p "${OUTDIR}" "${TMPDIR}" "${LOGDIR}"

echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "GEX: ${GEX}"
echo "Cell type: ${CELLTYPE}"
echo "BAM: ${BAM}"
echo "Barcode file: ${BCFILE}"
echo "Output bigWig: ${BW}"

if [[ ! -f "${BAM}" ]]; then
    echo "Missing BAM: ${BAM}"
    exit 1
fi

if [[ ! -f "${BCFILE}" ]]; then
    echo "Missing barcode file: ${BCFILE}"
    exit 1
fi

if [[ -f "${BW}" ]]; then
    echo "bigWig already exists, skipping: ${BW}"
    exit 0
fi

RAWBAM="${TMPDIR}/${CELLTYPE}.bam"

if [[ ! -f "${BAM}.bai" ]]; then
    samtools index -@ "${THREADS}" "${BAM}"
fi

samtools view -@ "${THREADS}" -h "${BAM}" | \
awk -v BCFILE="${BCFILE}" '
    BEGIN {
        while ((getline < BCFILE) > 0) {
            keep[$1] = 1
        }
    }
    /^@/ { print; next }
    {
        cb = ""
        for (i = 12; i <= NF; i++) {
            if ($i ~ /^CB:Z:/) {
                cb = substr($i, 6)
                break
            }
        }
        if (cb in keep) print
    }
' | samtools view -@ "${THREADS}" -b -o "${RAWBAM}" -

NREADS=$(samtools view -@ "${THREADS}" -c "${RAWBAM}")
if [[ "${NREADS}" -eq 0 ]]; then
    echo "No reads for ${GEX} ${CELLTYPE}, skipping" | tee -a "${LOGDIR}/empty_groups.log"
    rm -f "${RAWBAM}"
    exit 0
fi

samtools sort -@ "${THREADS}" -o "${SORTBAM}" "${RAWBAM}"
samtools index -@ "${THREADS}" "${SORTBAM}"
rm -f "${RAWBAM}"

~/.local/bin/bamCoverage \
    --bam "${SORTBAM}" \
    --outFileName "${BW}" \
    --outFileFormat bigwig \
    --binSize "${BINSIZE}" \
    --normalizeUsing CPM \
    --numberOfProcessors "${THREADS}" \
    --minMappingQuality 30 \
    --exactScaling \
    --skipNonCoveredRegions \
    > "${LOGDIR}/${CELLTYPE}.bamCoverage.stdout.log" \
    2> "${LOGDIR}/${CELLTYPE}.bamCoverage.stderr.log"

echo "Finished ${GEX} ${CELLTYPE}"
