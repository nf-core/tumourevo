#!/usr/bin/env python3

import sys
import pandas as pd
import argparse
## Read command line args

parser = argparse.ArgumentParser(description='Process input and output files.')
parser.add_argument('--joint', help='Joint table path')
parser.add_argument('--best_fit', help='Pyclone-vi best fit path')
parser.add_argument('--ctree_input', help='Table for ctree path')
#parser.add_argument('--pyclone_joint', help='Cluster-annotated joint table')
args = parser.parse_args()


## Read pyclone best fit table
#best_fit_file=sys.argv[2]
best_fit_file = args.best_fit
df_output = pd.read_csv(best_fit_file, sep='\t')

## Read pyclone input table
joint_table = args.joint
#df_input = pd.read_csv(sys.argv[1], sep = '\t')
df_input = pd.read_csv(joint_table, sep = '\t')
## Caluclate number of mutations per cluster and add to the subset orginal df_output dataframe

df_output_small = df_output[['mutation_id','sample_id','cluster_id','cellular_prevalence']]
df_output_small['nMuts'] = df_output_small.groupby('cluster_id')['mutation_id'].transform('nunique')

## Find the clonal cluster and add the colum 'is.clonal' to identify it

samples = df_output_small['sample_id'].unique()
top_clusters = 0
for s in samples:
  ind = df_output_small[df_output_small['sample_id']==s]['cellular_prevalence'].idxmax()
  top_cluster = df_output_small.loc[ind,'cluster_id']
  if (top_clusters != top_cluster):
    top_clusters = top_cluster

i = df_output_small[df_output_small['cluster_id']==top_clusters].index

df_output_small['is.clonal'] = "F"
df_output_small.loc[i,'is.clonal'] = "T"

## Merge the input information about driver genes to the clusters
if 'driver_label' not in df_input.columns or 'is_driver' not in df_input.columns:
    df_input['driver_label'] = "NA"
    df_input['is_driver'] = False
    df_input.loc[1,'driver_label'] = ""
    df_input.loc[1,'is_driver'] = True

df_merged = pd.merge(df_output_small, df_input[['mutation_id','sample_id','patient_id','driver_label','is_driver']], on = ["mutation_id","sample_id"], how="inner")
df_final = df_merged.drop_duplicates(subset=['sample_id', 'driver_label','cluster_id','is.clonal','is_driver'])
df_final = df_final.rename(columns={'cellular_prevalence':'CCF','cluster_id': 'cluster','driver_label':'variantID','patient_id':'patientID','is_driver':'is.driver'})
df_final['is.driver'] = df_final['is.driver'].replace({False: 'F', True: 'T'})

######### OLD SCRIPT #########
## Find clusters with driver genes and find cluster wihtout driver genes
## For the ones without drivers add the 'NA' to variantID and False to is.driver
#driver_genes_df = df_final[df_final['is_driver'] == True]

#non_driver_genes_df = df_final[df_final['is_driver'] == False].drop_duplicates(subset=['sample_id','cluster_id']).copy()
#non_driver_genes_df['gene'] = 'NA'
#non_driver_genes_df['is_driver'] = False
#result_df = pd.concat([driver_genes_df, non_driver_genes_df], ignore_index=True)

## Identify clusters for which is.driver is both True and False
#conflicting_clusters = result_df.groupby(['sample_id', 'cluster_id'])['is_driver'].nunique() == 2
#conflicting_clusters = conflicting_clusters[conflicting_clusters].index

## Create a new DataFrame with only True values for conflicting clusters
#true_values_df = result_df[(result_df['is_driver'] == True) & result_df.set_index(['sample_id', 'cluster_id']).index.isin(conflicting_clusters)]

## Create a new DataFrame with all values for non-conflicting clusters
#non_conflicting_df = result_df[~result_df.set_index(['sample_id', 'cluster_id']).index.isin(conflicting_clusters)]

## Concatenate the DataFrames
#result_df = pd.concat([true_values_df, non_conflicting_df], ignore_index=True)
#result_df['is_driver'] = result_df['is_driver'].replace({False: 'F', True: 'T'})

#result_df = result_df.rename(columns={'cellular_prevalence':'CCF','cluster_id': 'cluster'})
#result_df['tool'] = pd.Series(["pyclonevi" for x in range(len(result_df.index))])
#ctree_input = result_df[['patient_id','gene','is_driver','is.clonal','cluster','nMuts','sample_id','CCF','tool']]

## Final ctree input dataframe
final_jt_file = args.ctree_input
df_final.to_csv(final_jt_file, sep="\t",index=False, header=True)

## Annotate input joint table with pyclone-vi assignments
#df2 = pd.merge(df_input, df_output[['mutation_id','sample_id','cluster_id','cellular_prevalence','cluster_assignment_prob']], on=['mutation_id','sample_id'])
#df_joint_new=df2.rename(columns={"cluster_id": "pyclonevi_cluster_id", "cellular_prevalence": "pyclonevi_cellular_prevalence", "cluster_assignment_prob": "pyclone_vi_cluster_assignment_prob"})
#annotated_joint=args.pyclone_joint
#df_joint_new.to_csv(annotated_joint, sep="\t",index=False, header=True)
