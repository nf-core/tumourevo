#!/usr/bin/env python3

import os
import pandas as pd



def input_processing(data):
    new_columns = {'Project': "dataset_id", 'Genome': '$reference_genome', 'Type': "SOMATIC", 'mut_type': "SNP"}
    df = data.assign(**new_columns)
    df['chr'] = df['chr'].astype(str).str[3:]
    df = df.rename(columns={'Indiv': 'Sample', 'chr': 'chrom', 'from': 'pos_start', 'to': 'pos_end'})
    df["ID"] = df["Sample"]
    df = df.loc[:, ['Project', 'Sample', 'ID', 'Genome', 'mut_type', 'chrom', 'pos_start', 'pos_end', 'ref', 'alt', 'Type']]
    return df



def process_tsv_join(tsv_join):
    patients_tsv = tsv_join.split()
    # Read each file into a pandas DataFrame and ensure all columns are of type 'string'
    tables = []
    for p_table in patients_tsv:
        df = pd.read_csv(p_table, sep='\\t', dtype=str)
        tables.append(df)
    multisample_table = pd.concat(tables, ignore_index=True)
    return multisample_table
