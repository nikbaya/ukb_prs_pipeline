#!/usr/bin/env python3
"""
Calculate PRS

@author: nbaya (2021-08-20)
"""

def calc_scores(mt, beta_field: str, gt_field: str = 'GT', standardize_genotypes=False):
	"""
	Calculate polygenic scores using effect sizes `beta`, a row expression in `mt`.
	"""

    assert standardize_genotypes==False, "Error: standardize_genotypes=True not implemented yet."	
	mt = mt.annotate_cols(score = hl.agg.sum(mt[gt_field]*mt[beta_field]))
	return mt
