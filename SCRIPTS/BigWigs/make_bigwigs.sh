#!/usr/bin/bash
set -euo pipefail

THREADS=8
BINSIZE=25

RUNROOT="/data/chi/PROJECTS/24-418_Fraser_HIV_Reactivation/Bioinformatics/Runs/22MN7CLT3/RUN2"
BARCODE_ROOT="barcode_lists_by_GEX"

declare -A BAMS
BAMS[GEX1]="${RUNROOT}/3config_24_418_GEX1/outs/possorted_genome_bam.bam"
BAMS[GEX2]="${RUNROOT}/3config_24_418_GEX2/outs/possorted_genome_bam.bam"
BAMS[GEX3]="${RUNROOT}/3config_24_418_GEX3/outs/possorted_genome_bam.bam"
BAMS[GEX4]="${RUNROOT}/3config_24_418_GEX4/outs/possorted_genome_bam.bam"
BAMS[GEX5]="${RUNROOT}/3config_24_418_GEX5/outs/possorted_genome_bam.bam"
BAMS[GEX6]="${RUNROOT}/3config_24_418_GEX6/outs/possorted_genome_bam.bam"
BAMS[GEX7]="${RUNROOT}/3config_24_418_GEX7/outs/possorted_genome_bam.bam"
BAMS[GEX8]="${RUNROOT}/3config_24_418_GEX8/outs/possorted_genome_bam.bam"

for GEX in GEX1 GEX2 GEX3 GEX4 GEX5 GEX6 GEX7 GEX8; do

    BAM="${BAMS[$GEX]}"
    BC_DIR="${BARCODE_ROOT}/${GEX}"

    if [[ ! -f "${BAM}" ]]; then
        echo "Missing BAM for ${GEX}: ${BAM}"
        continue
    fi

    if [[ ! -d "${BC_DIR}" ]]; then
        echo "Missing barcode directory for ${GEX}: ${BC_DIR}"
        continue
    fi

    OUTDIR="${RUNROOT}/${GEX}"
    TMPDIR="${OUTDIR}/tmp_bams"
    LOGDIR="${OUTDIR}/logs"

    mkdir -p "${OUTDIR}" "${TMPDIR}" "${LOGDIR}"

    if [[ ! -f "${BAM}.bai" ]]; then
        samtools index -@ ${THREADS} "${BAM}"
    fi

    for LIST in "${BC_DIR}"/*.txt; do
        [[ -e "$LIST" ]] || continue

        CELLTYPE=$(basename "${LIST}" .txt)
        RAWBAM="${TMPDIR}/${CELLTYPE}.bam"
        SORTBAM="${TMPDIR}/${CELLTYPE}.sorted.bam"
        BW="${OUTDIR}/${CELLTYPE}.bw"

        echo "Processing ${GEX} ${CELLTYPE}"

        samtools view -@ ${THREADS} -h "${BAM}" | \
        awk -v BCFILE="${LIST}" '
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
        ' | samtools view -@ ${THREADS} -b -o "${RAWBAM}" -

        NREADS=$(samtools view -@ ${THREADS} -c "${RAWBAM}")
        if [[ "${NREADS}" -eq 0 ]]; then
            echo "No reads for ${GEX} ${CELLTYPE}, skipping" | tee -a "${LOGDIR}/empty_groups.log"
            rm -f "${RAWBAM}"
            continue
        fi

        samtools sort -@ ${THREADS} -o "${SORTBAM}" "${RAWBAM}"
        samtools index -@ ${THREADS} "${SORTBAM}"
        rm -f "${RAWBAM}"

        bamCoverage \
            --bam "${SORTBAM}" \
            --outFileName "${BW}" \
            --outFileFormat bigwig \
            --binSize ${BINSIZE} \
            --normalizeUsing CPM \
            --numberOfProcessors ${THREADS} \
            --minMappingQuality 30 \
            --exactScaling \
            --skipNonCoveredRegions \
            > "${LOGDIR}/${CELLTYPE}.bamCoverage.stdout.log" \
            2> "${LOGDIR}/${CELLTYPE}.bamCoverage.stderr.log"

    done
done

