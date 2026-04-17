#!/bin/bash
#SBATCH --job-name=group_bw_all
#SBATCH --output=group_bw_all.out
#SBATCH --error=group_bw_all.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

set -euo pipefail

module load r-xml/3.99-0.14-ey4x3fd

Rscript run_group_bw.R "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/work/a7/961ca94603d713bbcbc755966180e8/proj_scrnaseq_unconstrained.rds" 8
Rscript run_group_bw.R "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/work/d2/3855bf25254f5e25ece914412c904c/proj_scrnaseq_unconstrained.rds" 8
Rscript run_group_bw.R "/data/chi/TEMP/FOO/TEMP1/TEMP/scATACpipe/work/a7/05f2f92f55e67bc78318f4c3065114/proj_scrnaseq_unconstrained.rds" 8
