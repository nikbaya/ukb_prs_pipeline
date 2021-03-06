#!/usr/bin/env python3
"""
Get UK Biobank phenotype tables

@author: nbaya (2021-08-13)
"""

import argparse
import pandas as pd


def get_phenotype(input_path, sample_id_col, phenotype_map: dict, read_csv_kwargs):
    """Return reformatted phenotype data as a Pandas DataFrame

    :param input_path: Path to input file
    :param sample_id_col: Name of column with sample IDs
    :param phenotype_map: Values of dict can either be the name of the phenotype column to use, or
    they can be a list of phenotype column names to reduce into a single phenotype.
    :param read_csv_kwargs: kwargs for Pandas read_csv
    """

    df = pd.read_csv(input_path, **read_csv_kwargs)
    for final_col_name, col_names in phenotype_map:
        if type(col_names) is list:
            # Take union of boolean fields
            df[final_col_name] = reduce(lambda x, y: (x == 1) | (
                y == 1), [df[field] for field in col_names]).sum()
        else:
            df.rename(columns={col_names: final_col_name}, inplace=True)
    return df
