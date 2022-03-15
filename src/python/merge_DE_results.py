'''
This file merges all of the DESeq2 results in a folder. For each file it will extract the log2foldchange and adjust pvalue information.
Steps to run:
1. put all result files in a folder, the files needs to be ended with .result.tsv
2. run python merge_DE_results.py -p path -o out_file -g gene_file
The gene_file is optional. If provided, each line should have one gene. Then this code would extract information for the genes in the file 
'''

import argparse
import pandas as pd
import glob

parser = argparse.ArgumentParser(description='merge DESeq2 results')
parser.add_argument('-p','--path',action='store',dest='path',help='path that has all the DESeq2 results')
parser.add_argument('-o','--out',action='store',dest='out_fn',help='out merged file')
parser.add_argument('-g','--gene',action='store',dest='gene_fn',help='gene file if you want to extarct merged results for specific genes',default='')

args   = parser.parse_args()
path   = args.path
out_fn = args.out_fn
gene_fn = args.gene_fn


fns = sorted(glob.glob(f'{path}/*.result.csv'))
fns = [f.replace('\\','/') for f in fns]
if gene_fn == '':
	genes = ''
else:
	genes = pd.read_csv(gene_fn, header=None)[0].tolist()

dfs = []
gene_name_dic = {}
for fn in fns:
    df = pd.read_csv(fn,header=0,index_col=0)
    df.index.name = 'geneid'
    name = '.'.join(fn.split('/')[-1].split('.')[:-2])
    df.rename(columns={'log2FoldChange':f'{name}_log2FoldChange','padj':f'{name}_padj'},
              inplace=True)
    dfs.append(df[[f'{name}_log2FoldChange',f'{name}_padj']])
    dic = df['symbol'].to_dict()
    gene_name_dic.update(dic)
res = pd.concat(dfs, axis=1, sort=False)
res['symbol'] = res.index.map(lambda x: gene_name_dic[x])
res.index.name = 'geneid'

columns = ['symbol'] + res.columns.tolist()[:-1]
res = res[columns]

if genes != '':
	res = res.query('symbol in @genes')
res[columns].to_csv(out_fn,sep='\t')