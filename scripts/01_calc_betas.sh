#!/usr/bin/env bash
#
# Calculate estimated effect sizes for PRS
#
# Author: Nik Baya (2021-08-24)
#
#$ -N calc_betas
#$ -wd /well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline
#$ -o logs/calc_betas.log
#$ -e logs/calc_betas.errors.log
#$ -P lindgren.prjc
#$ -pe shmem 8
#$ -q short.qe
#$ -t 1-22

set -o errexit
set -o nounset
module purge

for util_type in {qsub,hail}; do
  source /well/lindgren/UKBIOBANK/nbaya/resources/ukb_utils/bash/${util_type}_utils.sh
done

# options
readonly chrom=$( get_chr ${SGE_TASK_ID} )
readonly method="prs_cs"

job_set_up() {
  module load Anaconda3/2020.07
  set_up_conda
  if [ "${method}" = "prs_cs" ]; then
    conda activate prs_cs
  fi
  export PYTHONPATH="/well/lindgren/UKBIOBANK/nbaya/resources/ukb_utils"
  export PYTHONPATH="$PYTHONPATH:/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline"
  export PYTHONPATH="$PYTHONPATH:/well/lindgren/UKBIOBANK/nbaya/resources/prs"
}

job_set_up

SECONDS=0

if [ "${method}" = "prs_cs" ]; then
  readonly prs_cs_dir="/well/lindgren/UKBIOBANK/nbaya/resources/prs/PRScs"
  readonly prs_cs_script="${prs_cs_dir}/PRScs.py"
  readonly sst_file="data/prs_cs/input_sumstats/t2d.tsv"
  # readonly tmp_sst_file="data/tmp/prs_cs/t2d_chr${chr}.tsv"
  readonly out_dir="data/prs_cs/out/t2d"

  mkdir -p ${out_dir}

  # TODO: Consider modifying PRS-CS code to read in compressed sumstats
  python3 ${prs_cs_script} \
    --ref_dir="data/prs_cs/ld_ref_panels/ldblk_ukbb_eur" \
    --bim_prefix="data/prs_cs/bim/ukb_imputed_v3" \
    --sst_file="data/prs_cs/input_sumstats/t2d.tsv" \
    --n_gwas=100000 \
    --out_dir=${out_dir} \
    --chrom=${chrom} \
  && print_update "Finished PRS-CS for chrom${chrom}" "$SECONDS"
fi
