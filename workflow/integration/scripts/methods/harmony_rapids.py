import logging
logging.basicConfig(level=logging.INFO)
from scipy.sparse import issparse
import rapids_singlecell as rsc
import cupy as cp

from utils import add_metadata
from utils_pipeline.io import read_anndata
from utils_pipeline.accessors import select_layer
from utils_pipeline.processing import process

input_adata = snakemake.input[0]
output_adata = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

logging.info(f'Read {input_adata}...')
adata_raw = read_anndata(input_adata)
adata_raw.X = select_layer(adata_raw, params['norm_counts'])

# subset to HVGs
adata_raw = adata_raw[:, adata_raw.var['highly_variable']]

# prepare for integration
adata = adata_raw.copy()
adata.X = select_layer(adata, params['norm_counts'], force_dense=True)

# run method
logging.info('Run rapids_singlecell harmony...')
rsc.tl.harmony_integrate(
    adata,
    key=wildcards.batch,
    basis='X_pca',
    adjusted_basis='X_emb'
)

# prepare output adata
adata = process(adata=adata, adata_raw=adata_raw, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

adata.write_zarr(output_adata)