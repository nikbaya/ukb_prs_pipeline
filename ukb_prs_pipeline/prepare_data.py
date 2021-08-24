#!/usr/bin/env python3
"""
Prepare data for use with PRS method

Author: Nik Baya (2021-08-24)
"""

import argparse
import hail as hl

from ukb_utils.hail import hail_bmrc_init

DATA_DIR = "/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline/data"

# NOTE: Consider using Pandas instead of Hail to make the dependencies easier


def munge_sumstats(ss_path: str, prs_method: str, out: str):
    """Prepare summary statistics for use with PRS method

    :param ss_path: Path to summary statistics to munge
    :param prs_method: The PRS method to be used
    :param out: Path to write output
    """

    hail_bmrc_init(
        log="logs/munge_sumstats-hail.log",
        default_reference="GRCh37"
    )

    if prs_method == "prs_cs":
        # NOTE: This block will change when we generalise to other traits

        ss = hl.import_table(
            paths=ss_path,
            types={'beta': hl.tfloat, 'pval': hl.tfloat}
        )
        variant_manifest = hl.import_table(
            paths=f"{DATA_DIR}/neale_lab_round2/variants.tsv.bgz",
            key="variant",
        )

        ss = ss.select(
            SNP=variant_manifest[ss.variant].rsid,
            A1=variant_manifest[ss.variant].ref,
            A2=variant_manifest[ss.variant].alt,
            BETA=ss.beta,
            P=ss.pval
        )

        ss.export(out)


def get_bim_prefix_path(dataset):
    """Get PLINK bim file prefix path for a dataset
    """
    if dataset == "ukb_imputed_v3":
        return f'{DATA_DIR}/prs_cs/bim/ukb_imputed_v3'


def get_bim(dataset):
    """Get bim file for use with PRS-CS
    """
    from ukb_utils.genotypes import get_ukb_imputed_v3_bgen, AUTOSOMES

    bim_prefix = get_bim_prefix_path(dataset)

    if not hl.hadoop_is_file(bim_prefix+'.bim'):
        mt = get_ukb_imputed_v3_bgen(chroms=AUTOSOMES)
        ht = mt.rows().key_by()
        ht = ht.select(
            CHR=ht.locus.contig,
            VARID=ht.rsid,
            CM=0,
            POS=ht.locus.position,
            A1=ht.alleles[0],
            A2=ht.alleles[1]
        )
        ht.export(bim_prefix+'.bim', header=False)


def main(args):
    if args.munge_sumstats:
        munge_sumstats(
            ss_path=args.ss_path,
            prs_method=args.prs_method,
            out=args.out
        )
   if args.get_bim:
       get_bim(
           dataset=args.dataset
       )

if __name__=="__main__":
   parser = argparse.ArgumentParser()
   parser.add_argument('--munge_sumstats', help='')
   parser.add_argument('--ss_path', help='')
   parser.add_argument('--prs_method', help='')
   parser.add_argument('--out', help='')
   parser.add_argument('--get_bim', help='')
   parser.add_argument('--dataset', help='')
   args = parser.parse_args()
   main(args)
