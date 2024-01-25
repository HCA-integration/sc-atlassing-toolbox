from pathlib import Path
import scanpy as sc
import logging
logging.basicConfig(level=logging.INFO)

from utils import add_metadata, remove_slots
from utils_pipeline.io import read_anndata, write_zarr_linked
from utils_pipeline.processing import assert_neighbors


input_file = snakemake.input[0]
output_file = snakemake.output[0]
wildcards = snakemake.wildcards
params = snakemake.params

adata = read_anndata(
    input_file,
    X='layers/norm_counts',
    obs='obs',
    var='var',
    obsp='obsp',
    obsm='obsm',
    varm='varm',
    uns='uns'
)

# prepare output adata
files_to_keep = ['obsm', 'uns']

if 'X_pca' not in adata.obsm:
    sc.pp.pca(adata, use_highly_variable=True, batch=wildcards.batch)
    files_to_keep.extend(['varm'])
adata.obsm['X_emb'] = adata.obsm['X_pca']

logging.info(adata.__str__())
logging.info(adata.uns.keys())
try:
    assert_neighbors(adata)
    logging.info(adata.uns['neighbors'].keys())
    files_to_keep.extend(['obsp', 'uns'])
except AssertionError:
    logging.info('Compute neighbors...')
    sc.pp.neighbors(adata)
    print(adata.uns['neighbors'])

adata = remove_slots(adata=adata, output_type=params['output_type'])
add_metadata(adata, wildcards, params)

# write file
logging.info(f'Write {output_file}...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    input_file,
    output_file,
    files_to_keep=['obsm', 'uns'],
    slot_map={'X': 'layers/norm_counts'},
)