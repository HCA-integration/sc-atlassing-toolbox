import numpy as np
import pandas as pd
import logging
logger = logging.getLogger('Run metric')
logger.setLevel(logging.INFO)

try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    logger.info('no hardware acceleration for sklearn')

from metrics import metric_map
from metrics.utils import write_metrics
from utils.io import read_anndata


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params
threads = snakemake.threads

dataset = wildcards.dataset
file_id = wildcards.file_id
batch_key = wildcards.batch
label_key = wildcards.label
metric = wildcards.metric

metric_type = params.get('metric_type')
assert metric_type in ['batch_correction', 'bio_conservation'], f'Unknown metric_type: {metric_type}'
allowed_output_types = params.get('output_types')
input_type = params.get('input_type')
comparison = params.get('comparison', False)
cluster_key = params.get('cluster_key', 'leiden')
metric_function = metric_map[metric]

uns = read_anndata(input_file, uns='uns').uns
output_type = uns.get('output_type', 'full') # Same as in prepare.py

if output_type not in allowed_output_types:
    logging.info(
        f'Skip metric={metric} for data output type={output_type}\n'
        f'allowed output types={allowed_output_types}'
    )
    write_metrics(
        scores=[np.nan],
        output_types=[output_type],
        metric=metric,
        metric_type=metric_type,
        batch=batch_key,
        label=label_key,
        file_id=file_id,
        dataset=dataset,
        filename=output_file,
        **uns.get('wildcards', {}),
        )
    exit(0)

kwargs = {'obs': 'obs', 'uns': 'uns'}
if input_type == 'knn':
    kwargs |= {'obsp': 'obsp'}
if input_type == 'embed':
    kwargs |= {'obsm': 'obsm'}
if input_type == 'full':
    kwargs |= {'X': 'X', 'var': 'var'}

logger.info(f'Read {input_file} of input_type {input_type}...')
adata = read_anndata(input_file, **kwargs)
print(adata, flush=True)
adata_raw = None

if comparison:
    adata_raw = read_anndata(
        input_file,
        X='raw/X',
        obs='obs',
        obsm='raw/obsm',
        var='raw/var',
        uns='raw/uns',
        dask=True,
        backed=True,
    )
    print('adata_raw')
    print(adata_raw, flush=True)

# prepare clustering columns
def replace_last(source_string, replace_what, replace_with, cluster_key='leiden'):
    """
    adapted from
    https://stackoverflow.com/questions/3675318/how-to-replace-some-characters-from-the-end-of-a-string/3675423#3675423
    """
    if not source_string.startswith(cluster_key):
        return source_string
    if not source_string.endswith(replace_what):
        return source_string + '_orig'
    head, _sep, tail = source_string.rpartition(replace_what)
    return head + replace_with + tail

# subset obs columns for metrics
cluster_columns = [col for col in adata.obs.columns if col.startswith(cluster_key) and col.endswith('_1')]
adata.obs = adata.obs[cluster_columns+[batch_key, label_key]].copy()
adata.obs.rename(columns=lambda x: replace_last(x, '_1', ''), inplace=True)

logger.info(f'Run metric {metric} for {output_type}...')
adata.obs[batch_key] = adata.obs[batch_key].astype(str).fillna('NA').astype('category')
score = metric_function(
    adata,
    output_type,
    batch_key=batch_key,
    label_key=label_key,
    adata_raw=adata_raw,
    cluster_key=cluster_key,
    n_threads=threads,
)

write_metrics(
    scores=[score],
    output_types=[output_type],
    metric=metric,
    metric_type=metric_type,
    batch=batch_key,
    label=label_key,
    file_id=file_id,
    dataset=dataset,
    filename=output_file,
    **adata.uns.get('wildcards', {}),
)
