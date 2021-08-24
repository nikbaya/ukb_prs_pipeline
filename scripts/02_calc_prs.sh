#!/usr/bin/env bash
#
# Calculate PRS
#
# Author: Nik Baya (2021-08-24)
#
#$ -N calc_prs
#$ -wd /well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline
#$ -o logs/calc_prs.log
#$ -e logs/calc_prs.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 8
#$ -q short.qe

set -o errexit
set -o nounset
module purge

for util_type in {qsub,hail}; do
  source /well/lindgren/UKBIOBANK/nbaya/resources/ukb_utils/bash/${util_type}_utils.sh
done

# options
readonly chrom=$( get_chr ${SGE_TASK_ID} )
readonly prs_method="prs_cs"

# Hail script
readonly hail_script="ukb_prs_pipeline/calc_prs.py"

job_set_up() {
  module load Anaconda3/2020.07
  set_up_hail
  export PYTHONPATH="/well/lindgren/UKBIOBANK/nbaya/resources/ukb_utils"
  export PYTHONPATH="$PYTHONPATH:/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline"
  export PYTHONPATH="$PYTHONPATH:/well/lindgren/UKBIOBANK/nbaya/resources/prs"
}

job_set_up

SECONDS=0

if [ "${prs_method}" = "prs_cs" ]; then
  readonly prs_cs_sst_prefix="data/prs_cs/out/t2d_pst_eff_a1_b0.5_phiauto"
  readonly out_dir="data/prs_cs/out/t2d"
  readonly out="${out_dir}/prs.ukb_imputed_v3.t2d.prs_cs.tsv.gz"

  mkdir -p ${out_dir}

  python3 ${hail_script} \
    --get_prs \
    --prs_method ${prs_method} \
    --beta_est_path ${prs_cs_sst_prefix} \
    --dataset "ukb_imputed_v3" \
    --out ${out} \
    && print_update "Finished writing PRS: ${out}" "${SECONDS}"

fi

