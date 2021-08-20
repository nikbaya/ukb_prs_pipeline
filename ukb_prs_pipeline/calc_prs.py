#!/usr/bin/env python3
"""
Calculate PRS

@author: nbaya (2021-08-20)
"""

import argparse


def run_prscs(chrom_list, ref_dir, bim_prefix, sst_file, n_gwas, a, b, phi, n_iter,
        n_burnin, thin, out_dir, beta_std, seed):

    from PRScs.parse_genet import parse_ref, parse_bim, parse_sumstats, parse_ldblk
    from PRScs.mcmc_gtb import mcmc

    for chrom in chrom_list:
        print('##### process chromosome %d #####' % int(chrom))

        ref_basename = ref_dir.split("/")[-1]

        if '1kg' in ref_basename:
            ref_dict = parse_ref(ref_dir + '/snpinfo_1kg_hm3', int(chrom))
        elif 'ukbb' in ref_basename):
            ref_dict=parse_ref(ref_dir + '/snpinfo_ukbb_hm3', int(chrom))

        vld_dict=parse_bim(
                bim_prefix, int(chrom)
                )

        sst_dict=parse_sumstats(
                ref_dict, vld_dict, sst_file, n_gwas
                )

        ld_blk, blk_size=parse_ldblk(
                ref_dir, sst_dict, int(chrom)
                )

        mcmc(
                a = a, b = b, phi = phi, sst_dict = sst_dict, n = n_gwas,
                ld_blk = ld_blk, blk_size = blk_size, n_iter = n_iter,
                n_burnin = n_burnin, thin = thin, chrom = chrom,
                out_dir = out_dir, beta_std = beta_std, seed = seed
        )


def main(args):
    pass


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output_path', default = None,
                        help = 'Output file path')
    args=parser.parse_args()

    main(args)
