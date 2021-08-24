#!/usr/bin/env python3
"""
Run full PRS pipeline

@author: nbaya (2021-08-20)
"""

import argparse
#  import get_genotypes
#  import calc_betas
#  import calc_prs
#  from get_phens import get_phenotype
#  import hail as hl

DATA_DIR = "/well/lindgren/UKBIOBANK/nbaya/ukb_prs_pipeline/data"


def t2d_test(chrom):
    import os
    from ukb_utils.genotypes import get_ukb_imputed_v3_bgen, AUTOSOMES
    from ukb_prs_pipeline.calc_betas import run_prscs

    chrom_list=[int(chrom)]

    #  hl.init(
    #      log="logs/t2d_hail.log",
    #      default_reference="GRCh37",
    #      master='local[%s]' % os.environ.get('NSLOTS', 1)
    #  )

    ss_prscs_path = f'{DATA_DIR}/prs_cs/input_sumstats/t2d.tsv'
    #  if not hl.hadoop_is_file(ss_prscs_path):
    #      ss = hl.import_table(
    #          paths=f'{DATA_DIR}/neale_lab_round2/2443.gwas.imputed_v3.both_sexes.tsv.bgz',
    #          types={'beta':hl.tfloat, 'pval':hl.tfloat}
    #      )

    #      # TODO: Remove after testing
    #      # ss = ss.annotate(**hl.parse_variant(ss.variant, reference_genome="GRCh37"))
    #      # ss = ss.filter(hl.literal(list(map(str, chroms))).contains(ss.locus.contig))

    #      variant_manifest = hl.import_table(
    #          paths=f"{DATA_DIR}/neale_lab_round2/variants.tsv.bgz",
    #          key="variant",
    #      )

    #      ss = ss.select(
    #          SNP = variant_manifest[ss.variant].rsid,
    #          A1 = variant_manifest[ss.variant].ref,
    #          A2 = variant_manifest[ss.variant].alt,
    #          BETA = ss.beta,
    #          P = ss.pval
    #      )

    #      ss.export(ss_prscs_path)
    #
    bim_prefix = f'{DATA_DIR}/prs_cs/bim/ukb_imputed_v3'
    #  if not hl.hadoop_is_file(bim_prefix+'.bim'):
    #      mt = get_ukb_imputed_v3_bgen(chroms=AUTOSOMES)
    #      ht = mt.rows().key_by()
    #      ht = ht.select(
    #          CHR = ht.locus.contig,
    #          VARID = ht.rsid,
    #          CM = 0,
    #          POS = ht.locus.position,
    #          A1 = ht.alleles[0],
    #          A2 = ht.alleles[1]
    #      )
    #      ht.export(bim_prefix+'.bim', header=False)

    # path for table of estimated betas from PRS-CS
    prscs_out_dir = f'{DATA_DIR}/prs_cs/output/t2d'
    n_gwas = 100000  # GWAS sample size (approximate for now)
    run_prscs(
        chrom_list=chrom_list,
        ref_dir=f'{DATA_DIR}/prs_cs/ld_ref_panels/ldblk_ukbb_eur',
        bim_prefix=bim_prefix,
        sst_file=ss_prscs_path,
        n_gwas=n_gwas,
        out_dir=prscs_out_dir
    )

    #  prscs_df = pd.from_dict({k:v for k,v in sst_dict if k in {'SNP','BP','A1','A2','beta_est'}})
    #  prscs_ht = hl.Table.from_pandas(df)
    #  prscs_ht.write(prscs_out_path)

    pass

    mt = get_ukb_imputed_v3_bgen(chrom=chroms)

    case_ids = hl.import_table(
        paths="/well/lindgren/UKBIOBANK/for_nik/t2d_case_ids.txt",
        key='f0',
        no_header=True,
    )
    case_ids = case_ids.rename({'f0': 's'})


def main(args):
    #hl.init(log="logs/main-hail.log", default_reference="GRCh38")
    t2d_test(chrom=args.chrom)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--chrom', help='Chromosome to run')
    args = parser.parse_args()

    main(args)
