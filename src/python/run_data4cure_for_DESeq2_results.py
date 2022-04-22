'''Created by Shangzhong.Li@pfizer.com at 2022/04
'''
import pandas as pd
import argparse
import os
import logging.config

parser = argparse.ArgumentParser(description='transform the condition table for DESeq2 input')

parser.add_argument('-i','--input',action='store',dest='in_file',help='input DESeq2 results file')
parser.add_argument('-o','--output',action='store',dest='out_file',help='output tbl file')
parser.add_argument('--domain',action='store',dest='domain',help='domain which is used to describe the experiment')

# d4c specific parameters
parser.add_argument('--user',action='store',dest='user',help='user name for d4c')
parser.add_argument('--pw',action='store',dest='pw',help='password for d4c')

args = parser.parse_args()
in_fn = args.in_file
out_tbl = args.out_file
domain = args.domain

user = args.user
pw = args.pw


from d4client.common import AuthProvider    
from d4client.ops_manager import makeOpsManager
import d4client.d4capp as d4capp

server_url = 'https://data4cure.pfizer.com'
ap = AuthProvider(server_url, verify_ssl=False)
ap.login(user, pw)
token = ap.get_token({ 'user': True })


def format_ETSV2tbl(in_fn, out_fn):
    head = {}
    head['data_type.name'] = ['numeric'] * 6
    head['data_type.type'] = ['numeric'] * 6
    head['entity.type']    =['generic_entity'] * 6
    head['entity.domain.name'] = ['gene_level_results'] * 6
    head['entity.label'] = ['baseMean','foldChange','lfcSE','stat','pvalue','padj']
    head['entity.name'] = head['entity.label']
    head['label'] = head['entity.name']

    final_head = [['#ETSV Gene Statistics','data_type.name'] + head['data_type.name'],
                ['', 'data_type.type'] + head['data_type.type'],
                ['', 'entity.type'] + head['entity.type'],
                ['', 'entity.domain.name'] + head['entity.domain.name'],
                ['', 'entity.label'] + head['entity.label'],
                ['', 'entity.name'] + head['entity.name'],
                ['', 'label'] + head['label'],
                ['gene name'] + [''] * 7
    ]

    in_df = pd.read_csv(in_fn,header=0,index_col=0)
    in_df['foldChange'] = 2 ** in_df['log2FoldChange']
    in_df = in_df[head['entity.label']]
    if 'symbol' in in_df.columns:
        del in_df['symbol']

    with open(out_fn,'w') as f:
        for i in final_head:
            f.write('\t'.join(i) + '\n')
    in_df.insert(0,'insert','')
    in_df.to_csv(out_fn,sep='\t',header=None,mode='a')


def format_file2tbl(in_fn, out_tbl):
    in_df = pd.read_csv(in_fn,header=0,index_col=0)
    in_df.index.name = 'gene'
    in_df['foldChange'] = 2 ** in_df['log2FoldChange']
    columns = ['baseMean','foldChange','lfcSE','stat','pvalue','padj']
    in_df = in_df[columns]
    inter_fn = f'{in_fn}.inter'
    in_df.to_csv(inter_fn,sep='\t')
    cmd = f'import_data_table genestats --use-ensembl {inter_fn} {out_tbl}'
    os.system(cmd)
    os.remove(inter_fn)


def upload_file2_d4c(DESeq2_file, platform_path):
    # upload file to d4c
    table_loader = d4capp.FileLoader(
            server_url,
            token=token
        )
    try:
        table_loader.upload_file(
            DESeq2_file,
            platform_path,
        )
    except d4capp.UploadTsvTableError as e:
        print(e.error_details())
        # Above prints errors related to table import, if any. This is equivalent to the errors shown via the `info` button found adjacent to uploaded tables in DataHub.
        raise e

def run_pathway_expression(deseq2_tbl, platform_path, domain):
    m = makeOpsManager(server_url=server_url, auth_provider=ap)
    app = m.applications_ops.application_by_name(
    m.applications_ops.PATHWAY_EXPRESSION_ANALYSIS
    )
    base_name = os.path.basename(deseq2_tbl[:-4])
    context = {
        'entity_measurement':{
            'entity':{
                'type': 'health_condition',
                'name': base_name,
                'label': base_name,
                'domain': {'name':'health_condition'},
            },
            'label': base_name,
        }
    }

    spec = {
        'name_field': domain,
        'description_field': f'Pathway analysis of DESeq2 results for {deseq2_tbl}',
        'gene_stats_table': 'My Data/'+platform_path+ '/' + os.path.basename(deseq2_tbl), 
        'p_value_col': ['pvalue'],
        'q_value_col': ['padj'],
        'fold_change_col': ['foldChange'],
        'test_stat_col': ['stat'],
        'results_dir_selector': 'My Data/' +platform_path+'/DESeq2_pathway_analysis',
        'domains': ['REACTOME', 'KEGG', "GO Biological Process", "GO Cellular Component", "GO Molecular Function", "MSigDB Canonical Pathways", "MSigDB Hallmark Gene Sets", "TF Targets HC","Immune Signatures","Pfizer Gene Sets v5"],
        'input_data_type': 'Gene Expression Statistics',
        'condition_entity_measurement': [context]
    }
    submission_spec = app.application_ops().submission_spec_from_dict(spec)
    res = m.applications_ops.submit(submission_spec)

    job, job_uuid = None, None
    if res.has_error():
        print('Error:', res.error)
    else:
        job = res.job
        job_uuid = res.job_uuid


LOGGING_CONFIG = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'simple': {
            'format': '%(levelname)s %(asctime)s %(message)s',
        },
    },
    'handlers': {
        'stdout': {
            'class': 'logging.StreamHandler',
            'formatter': 'simple',
            'level': 'DEBUG',
            'stream': 'ext://sys.stdout',
        },
    },
    'loggers': {
        'd4client': {
            'handlers': [ 'stdout', ],
            'level': 'INFO',
        },
    },
}
logging.config.dictConfig(LOGGING_CONFIG)

if __name__ == '__main__':
    # 1. format the table
    format_file2tbl(in_fn, out_tbl)
    # 2.1 submit to the data4cure using api
    platform_path = 'bulkRNASeq/' + os.path.basename(out_tbl).split('.')[0] + '_DESeq2_result'
    upload_file2_d4c(out_tbl, platform_path)
    # 2.2 run data4cure
    run_pathway_expression(out_tbl, platform_path, domain)


