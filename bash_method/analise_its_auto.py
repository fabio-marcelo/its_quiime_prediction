#!/usr/bin/env python3

# Import libraries
#!conda install -c conda-forge pandas -y
#!conda install -c conda-forge matplotlib -y

import pandas as pd
import re
import numpy as np
import itertools as IT
import glob
import os
import sys

path = sys.argv[1] 
file_final = path.split('/')[-1]
#final_table = 'final_table.csv'
final_table = sys.argv[2]

def planilha(x):
    
    #name = x
    file = x+'-final_output.csv'

    

    df = pd.read_csv(f'{path}/{file}', sep = '\t', 
                    skiprows=1, header = 0)
    df_filter1 = df.query("taxonomy != 'Unassigned'")
    df_filter1 = df.loc[df['taxonomy'] != 'Unassigned']
    

    df_filter1_indexed = df_filter1.set_index('taxonomy')
    df_filter1_indexed  

    df_transposed = df_filter1_indexed.transpose()
    df_transposed

    taxons = list(df_transposed.columns)
    count_taxons = pd.melt(df_transposed, value_vars=taxons,value_name='Count', ignore_index=False)
    count_taxons = count_taxons.query("Count != 0")
    count_taxons = count_taxons.sort_index()
    count_taxons['Tool'] = x
    count_taxons

    count_taxons.to_csv(f'{path}/count-taxons_{x}.csv',
                                   index=True, sep= '\t')

    df_result = df_transposed.where(df_transposed == 0.0, 
                                           df_transposed.columns.to_series(), 
                                           axis=1)

    # remover colunas e linhas contendo somente zero
    df_result_final = df_result.loc[:, (df_result != 0).any(axis=0)]
    df_result_final  = df_result_final.loc[(df_result_final !=0).any(axis=1)]
    # remove zeros
    df_result_final.replace(to_replace = 0.0, value = 'Nan', inplace=True)
    # remove NA
    # df_result_final.replace(re.compile('NA|,| '), '', inplace=True)
    # remove linhas completamente vazias
    df_result_final.replace('', np.nan, inplace=True)
    df_result_final.dropna(how='all')
    df_result_final.stack()

    # Define a function to keep unique values in a row
    def keep_unique(row):
        unique_values = []
        for item in row:
            if isinstance(item, list):
                unique_values.extend([val for val in item if val not in unique_values])
            elif isinstance(item, str):
                if item not in unique_values:
                    unique_values.append(item)
            else:
                if str(item) not in unique_values:
                    unique_values.append(str(item))
        return unique_values

    # Apply the function using apply and axis=1
    df_result_final1 = df_result_final.apply(keep_unique, axis=1)

    df_result_final2 = pd.DataFrame(df_result_final1)
    df_result_final2.columns =['taxonomy']

    df_result_final2.reset_index(inplace=True)
    df_result_final2 = df_result_final2.rename(columns = {'index':'sample'})

    df_split = pd.DataFrame(df_result_final2['taxonomy'].tolist()).fillna('')
    df_result_final3 = pd.concat([df_result_final2, df_split], axis=1)

    df_indexed = df_result_final3.set_index('sample')
    df_final = df_indexed.loc[:, df_indexed.columns!='taxonomy']

    df_final = (df_final.stack()
    .reset_index(level=1, drop=True)
    .reset_index(name='Taxonomy')
    )

    df_final['Taxonomy'].replace('', np.nan, inplace=True)
    df_final['Taxonomy'].replace('Nan', np.nan, inplace=True)
    df_final.dropna(subset=['Taxonomy'], inplace=True)

    df_final.to_csv(f'{path}/planilha_{x}.csv',
                                   index=True, sep= '\t')


# define methods
metodos = ["blast", "vsearch", "sklearn"]

# run def planilha over methods
for i in metodos:
    planilha(i)

blast = pd.read_csv(f'{path}/planilha_blast.csv', sep = '\t', index_col = 'sample')
# blast.rename(columns = {'Taxonomy':'blast'}, inplace = True)
blast['tool'] = 'blast'

vsearch = pd.read_csv(f'{path}/planilha_vsearch.csv', sep = '\t', index_col = 'sample')
# blast.rename(columns = {'Taxonomy':'blast'}, inplace = True)
vsearch['tool'] = 'vsearch'

sklearn = pd.read_csv(f'{path}/planilha_sklearn.csv', sep = '\t', index_col = 'sample')
# blast.rename(columns = {'Taxonomy':'blast'}, inplace = True)
sklearn['tool'] = 'sklearn'

# Delete column
del sklearn['Unnamed: 0']
del blast['Unnamed: 0']
del vsearch['Unnamed: 0']

frames = [blast, vsearch, sklearn]

result = pd.concat(frames)
# result = blast.join(sklearn).join(vsearch)
np.unique(result[['Taxonomy', 'tool']].values)
result = result.sort_index()

result.to_csv(f'{path}/final-{file_final}.csv',
                                   index=True, sep= '\t')

files = glob.glob(f'{path}/count-taxons*')
dfs = [pd.read_csv(f, header=0, sep="\t") for f in files]

total_counts = pd.concat(dfs,ignore_index=False)
total_counts = total_counts.groupby(['Unnamed: 0', 'taxonomy', 'Tool']).sum('Count')

# save as csv
total_counts.to_csv(f'{path}/{final_table}',
                                   index=True, sep= '\t')
