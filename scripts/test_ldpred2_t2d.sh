#!/usr/bin/env bash
#
# Test run of LDpred2 for T2D
#
# Author: Nik Baya (2021-08-27)
#
#$ -N test_ldpred2
#$ -wd /well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline
#$ -o logs/ldpred2/
#$ -e logs/ldpred2/
#$ -P lindgren.prjc
#$ -pe shmem 24
#$ -q short.qe
#$ -t 22

set -o errexit
set -o nounset
module purge

source /well/lindgren/UKBIOBANK/nbaya/resources/ukb_utils/bash/qsub_utils.sh

readonly chr=$( get_chr ${SGE_TASK_ID} )

# R script
readonly script="utils/test_ldpred2.R"

job_set_up() {
  module load R/3.6.2-foss-2019b
}

job_set_up

SECONDS=0
Rscript ${script}
print_update "Finished LDpred2 for chrom22" "$SECONDS"
