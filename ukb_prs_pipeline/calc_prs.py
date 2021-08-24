#!/usr/bin/env python3
"""
Calculate PRS

@author: nbaya (2021-08-20)
"""

import argparse
import hail as hl


def calc_scores(mt, beta_field: str, gt_field: str = 'GT', standardize_genotypes=False):
    """Calculate polygenic scores using effect sizes `beta`, a row expression in `mt`.
    """

    assert standardize_genotypes == False, "Error: standardize_genotypes=True not implemented yet."
    assert mt[
        gt_field].dtype == hl.tcall, f"Error: gt_field={gt_field} is not a CallExpression (type: {mt[gt_field].dtype})"
    mt = mt.annotate_cols(score=hl.agg.sum(
        mt[gt_field].n_alt_alleles()*mt[beta_field]))
    return mt


def get_prs(prs_method, beta_est_path, dataset, out):
    """Write out PRS for a dataset

    :param prs_method: The PRS method used to get the estimated effect sizes.
    :param beta_est_path: Path to table of estimated effect sizes to use
    :param dataset: The dataset for which to calculate PRS
    :param out: Path to which the PRS will be written
    """
    from ukb_utils.genotypes import get_ukb_imputed_v3_bgen_path, get_ukb_imputed_v3_bgen, AUTOSOMES

    # TODO: Generalise this later
    if prs_method == "prs_cs":
        beta_est = hl.import_table(
            paths=[f'{beta_est_path}_chr{chrom}.txt' for chrom in AUTOSOMES],
	    no_header=True,
	    types={"f5":hl.tfloat}
        )
        beta_est = beta_est.rename(
            {'f0': 'chr',
             'f1': 'rsid',
             'f2': 'pos',
             'f3': 'a1',
             'f4': 'a2',
             'f5': 'beta'
             }
        )
        beta_est = beta_est.key_by('rsid')

    if dataset == "ukb_imputed_v3":
        mt = get_ukb_imputed_v3_bgen()

    mt = mt.annotate_rows(
        beta_est=beta_est[mt.varid].beta
    )

    mt = mt.filter_rows(hl.is_defined(beta_est[mt.varid]))
    row_ct = mt.count_rows()
    print(f'\nNumber of variants shared between MatrixTable and sumstats: {row_ct}')

    mt = calc_scores(
        mt=mt,
        beta_field="beta_est"
    )

    mt.cols().select("score").export(out)


def main(args):
    from ukb_utils.hail import hail_bmrc_init

    hail_bmrc_init(
        log="logs/calc_prs-hail.log",
        default_reference="GRCh37"
    )

    if args.get_prs:
        get_prs(
            prs_method=args.prs_method,
            beta_est_path=args.beta_est_path,
            dataset=args.dataset,
            out=args.out
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--get_prs', default=False,
                        action="store_true", help='Whether to get PRS')
    parser.add_argument('--prs_method', default=None,
                        help='PRS method used to get estimated effect sizes')
    parser.add_argument('--beta_est_path', default=None,
                        help='Path to estimated effect sizes')
    parser.add_argument('--dataset', default=None,
                        help='Which dataset to use when calculating PRS')
    parser.add_argument('--out', default=None,
                        help='Path to write output file')
    args = parser.parse_args()

    main(args)
