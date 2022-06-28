'''Created by Shangzhong.Li@pfizer.com at 2022/04
'''
import pandas as pd
import numpy as np
import os
import yaml
import argparse
import logging.config

parser = argparse.ArgumentParser(description='transform the condition table for DESeq2 input')

parser.add_argument('-i','--input_tsv',action='store',dest='in_file',help='input expression file, row is gene, column is sample')
parser.add_argument('-m','--meta',action='store',dest='meta',help='meta file')
parser.add_argument('-c','--cmp',action='store',dest='cmp',help='compare file, each row is a comparision, has two columns: control and test')
parser.add_argument('-v','--vari',action='store',dest='variables',help='columns in meta_fn separated by ,')
parser.add_argument('-d','--dt',action='store',dest='data_type',help='data type of expression file, can be tpm or count')
parser.add_argument('--domain',action='store',dest='domain',help='domain which is used to describe the experiment')

# d4c specific parameters
parser.add_argument('--user',action='store',dest='user',help='user name for d4c')
parser.add_argument('--pw',action='store',dest='pw',help='password for d4c')


args = parser.parse_args()
expr_fn = args.in_file
meta_fn = args.meta
cmp_fn = args.cmp
variables = args.variables
conditions = variables.split(',')
data_type = args.data_type
domain = args.domain

user = args.user
pw = args.pw


def log2_transform(df):
    df = np.log2(df + 1)
    return df


def sub_meta(meta_fn, sub_meta_fn, conditions, ctrl, test):
    meta_df = pd.read_csv(meta_fn,sep='\t|,',header=0,engine='python')
    meta_df['mergeCond'] = meta_df[conditions].apply(lambda x: ':'.join(x) ,axis=1)
    sub = [ctrl, test]
    meta_df = meta_df.query('mergeCond in @sub')
    meta_df = meta_df.reset_index(drop=True)
    meta_df['mergeCond'] = meta_df['mergeCond'].map(lambda x: 0 if x == ctrl else 1)
    meta_df.to_csv(sub_meta_fn, sep='\t', index=False)
    samples = meta_df['Sample'].tolist()
    return samples


def sub_expr(expr_fn, sub_expr_fn, samples, data_type):
    '''extract sub expression matrix of interested samples'''
    expr_df = pd.read_csv(expr_fn,sep='\t|,',header=0,index_col=0,engine='python')
    expr_df = expr_df.round().astype(int)
    expr_df.index.name = 'gene'
    if data_type == 'tpm':
        expr_df = log2_transform(expr_df)
    expr_df[samples].T.to_csv(sub_expr_fn,sep='\t')


def get_config(config_fn, ctrl, test):
    '''prepare configure file'''
    ctrl_name = '_'.join(ctrl.split(':'))
    test_name = '_'.join(test.split(':'))
    config = {'mergeCond': {
            'name': 'condition',
            'data_type_name': 'factor',
            'data_type_type': 'numeric',
            'factor_levels': [
                {'value': 0, 'label': ctrl_name},
                {'value': 1, 'label': test_name}
            ],
            'entity_type': 'generic_entity'
        }
    }
    with open(config_fn, 'w') as f:
        yaml.dump(config, f, sort_keys=False,default_flow_style=False, encoding='utf-8')
    

def format_files(expr_fn, meta_fn, conditions, data_type, domain, ctrl, test):
    '''
    * expr_fn: expression file, row is gene, column is sample
    * meta_fn: meta data file, row is sample, column is condition, First column needs to be Sample.
    * conditions: a list of columns in meat_fn to combine together to define case or control
    * data_type: data type for expr_fn, can be 'count', 'tpm'.

    '''
    path = os.path.dirname(expr_fn)
    config_fn = f'{path}/config.txt'

    ctrl_name = '_'.join(ctrl.split(':'))
    test_name = '_'.join(test.split(':'))

    # get sub meta fn
    sub_meta_fn = f'{path}/sub_meta.tsv'
    sub_samples = sub_meta(meta_fn, sub_meta_fn, conditions, ctrl, test)
    # get sub expression fn
    sub_expr_fn = f'{path}/sub_expr.tsv'
    sub_expr(expr_fn, sub_expr_fn, sub_samples, data_type)
    # transform expression data
    expr_tbl = f'{path}/{test_name}_VS_{ctrl_name}.expr.tbl'
    cmd = f'import_data_table omics --use-ensembl {sub_expr_fn} rna_expression {domain} {expr_tbl}'
    os.system(cmd)
    # transform meta data
    config_fn = f'{path}/config.yaml'
    get_config(config_fn, ctrl, test)
    meta_tbl = f'{path}/{test_name}_VS_{ctrl_name}.meta.tbl'
    cmd = f'import_data_table metad {sub_meta_fn} {config_fn} {domain} {meta_tbl}'
    os.system(cmd)
    # return transformed file name
    return expr_tbl, meta_tbl


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


from d4client.common import AuthProvider    
from d4client.ops_manager import makeOpsManager
import d4client.d4capp as d4capp

server_url = 'https://data4cure.pfizer.com'
ap = AuthProvider(server_url, verify_ssl=False)
ap.login(user, pw)
token = ap.get_token({ 'user': True })

def upload_file2_d4c(ExprFilePath, MetaFilePath, platform_path):
    # upload file to d4c
    table_loader = d4capp.FileLoader(
            server_url,
            token=token
        )
    try:
        table_loader.upload_file(
            MetaFilePath,
            platform_path,
        )
    except d4capp.UploadTsvTableError as e:
        print(e.error_details())
        # Above prints errors related to table import, if any. This is equivalent to the errors shown via the `info` button found adjacent to uploaded tables in DataHub.
        raise e
    try:
        table_loader.upload_file(
            ExprFilePath,
            platform_path,
        )
    except d4capp.UploadTsvTableError as e:
        print(e.error_details())
        # Above prints errors related to table import, if any. This is equivalent to the errors shown via the `info` button found adjacent to uploaded tables in DataHub.
        raise e


def run_pathway_expression(expr_tbl, meta_tbl, platform_path, domain, test_name, ctrl_name, data_type):
    m = makeOpsManager(server_url=server_url, auth_provider=ap)
    app = m.applications_ops.application_by_name(
    m.applications_ops.PATHWAY_EXPRESSION_ANALYSIS
    )
    if data_type == 'count':
        spec = {
            'name_field': domain,
            'description_field': f'DESeq2 analysis between {test_name} (case) vs {ctrl_name} (control)',
            'condition_table': 'My Data/'+platform_path+ '/' + os.path.basename(meta_tbl),
            'gene_counts_table': 'My Data/'+platform_path+ '/' + os.path.basename(expr_tbl),
            'condition_var': ['condition'],
            'results_dir_selector': 'My Data/' +platform_path+'/DESeq2_analysis',
            'domains': ['REACTOME', 'KEGG', "GO Biological Process", "GO Cellular Component", "GO Molecular Function", "MSigDB Canonical Pathways", "MSigDB Hallmark Gene Sets", "TF Targets HC","Immune Signatures","Pfizer Gene Sets v6"],
            'input_data_type': 'Gene Expression (counts)'
        }
    elif data_type == 'tpm':
        spec = {
            'name_field': domain,
            'description_field': f'limma analysis between {test_name} (case) vs {ctrl_name} (control)',
            'condition_table': 'My Data/'+platform_path+ '/' + os.path.basename(meta_tbl),
            'gene_counts_table': 'My Data/'+platform_path+ '/' + os.path.basename(expr_tbl),
            'condition_var': ['condition'],
            'results_dir_selector': 'My Data/' +platform_path+'/DESeq2_analysis',
            'domains': ['REACTOME', 'KEGG', "GO Biological Process", "GO Cellular Component", "GO Molecular Function", "MSigDB Canonical Pathways", "MSigDB Hallmark Gene Sets", "TF Targets HC","Immune Signatures","Pfizer Gene Sets v6"],
            'input_data_type': 'Gene Expression (continuous)'
        }
    submission_spec = app.application_ops().submission_spec_from_dict(spec)
    res = m.applications_ops.submit(submission_spec)
    job, job_uuid = None, None
    if res.has_error():
        print('Error:', res.error)
    else:
        job = res.job
        job_uuid = res.job_uuid




if __name__ == '__main__':
    cmp_df = pd.read_csv(cmp_fn,sep='\t|,',header=0,engine='python')
    for idx, row in cmp_df.iterrows():
        ctrl = row['control']
        test = row['test']
        ctrl_name = '_'.join(ctrl.split(':'))
        test_name = '_'.join(test.split(':'))
        # 1. format files
        expr_tbl, meta_tbl = format_files(expr_fn, meta_fn, conditions, data_type, domain, ctrl, test)
        # 2. upload to d4c
        platform_path = f'BulkRNASeq/{domain}/{test_name}_VS_{ctrl_name}'
        upload_file2_d4c(expr_tbl, meta_tbl, platform_path)
        # # 3. run pathway expression
        run_pathway_expression(expr_tbl, meta_tbl, platform_path, domain, test_name, ctrl_name, data_type)

