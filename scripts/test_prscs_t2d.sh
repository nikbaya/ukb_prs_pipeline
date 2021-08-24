#!/usr/bin/env bash
#
# Test run of PRS-CS for T2D
#
# Author: Nik Baya (2021-08-23)
#
#$ -N test_prscs
#$ -wd /well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline
#$ -o logs/prs_cs.log
#$ -e logs/prs_cs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 4
#$ -q short.qe
#$ -t 1-22

set -o errexit
set -o nounset
module purge

for util_type in {qsub,hail}; do
  source /well/lindgren/UKBIOBANK/nbaya/resources/ukb_utils/bash/${util_type}_utils.sh
done

# options
readonly chr=$( get_chr ${SGE_TASK_ID} )

# Python script
readonly script="ukb_prs_pipeline/main.py"

job_set_up() {
  module load Anaconda3/2020.07
  set_up_conda
  conda activate prs_cs
  export PYTHONPATH="/well/lindgren/UKBIOBANK/nbaya/resources/ukb_utils"
  export PYTHONPATH="$PYTHONPATH:/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline"
  export PYTHONPATH="$PYTHONPATH:/well/lindgren/UKBIOBANK/nbaya/resources/prs"
}

job_set_up

SECONDS=0
python3 ${script} \
  --chrom ${chr}
print_update "Finished PRS-CS for chrom${chr}" "$SECONDS"
