#!/usr/bin/env python3
"""
Get genotypes

@author: nbaya (2021-08-20)
"""

import hail as hl

AUTOSOMES = list(range(1, 23))


def get_ukb_imputed_v3_bgen(chroms: list = AUTOSOMES, entry_fields: list = ['GT', 'GP', 'dosage']):
    # Use these bgens because I have write access to create .idx2 directorie
    # original link: /well/ukbb-wtchg/v3/imputation/ukb_imp_chr*_v3.bgen
    def bgen_path_fn(chrom): return f"/well/lindgren/UKBIOBANK/nbaya/resources/ref/ukb_wes_200k/ukb_ref_panel/data/imputed_v3_bgen/ukb_imp_chr{chrom}_v3.bgen"

    bgen_path_list = list(map(bgen_path_fn, chroms))

    for chrom, bgen_path in zip(chroms, bgen_path_list):
        if not hl.hadoop_is_dir(bgen_path+".idx2"):
            hl.index_bgen(
	        path=bgen_path,
		reference_genome="GRCh37",
		contig_recoding={str(chrom).zfill(2): str(chrom)},
	    )

    mt = hl.import_bgen(
        path=bgen_path_list,
        entry_fields=entry_fields,
        sample_file="/well/lindgren/UKBIOBANK/DATA/SAMPLE_FAM/ukb11867_imp_chr1_v3_s487395.sample"

    )
    return mt
