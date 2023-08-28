"""
Compute UMAP
"""
from pathlib import Path
import logging
logging.basicConfig(level=logging.INFO)
import scanpy as sc

from utils.io import read_anndata, link_zarr


def check_and_update_neighbors_info(adata, neighbors_key, params):
    # check if representation is available
    try:
        assert neighbors_key in adata.uns
        assert adata.uns[neighbors_key]['connectivities_key'] in adata.obsp
        assert adata.uns[neighbors_key]['distances_key'] in adata.obsp
    except AssertionError as e:
        logging.info(f'neighbors_key: {neighbors_key}')
        logging.info(adata.__str__())
        raise e

    # check if representation is available in current anndata, otherwise read from file
    use_rep = adata.uns[neighbors_key]['params'].get('use_rep', None)
    if use_rep not in adata.obsm:
        logging.info(f'Read {input_rep}...')
        adata_rep = read_anndata(input_rep)
        if use_rep == 'X' or use_rep is None:
            adata.obsm[use_rep] = adata_rep.X
        else:
            adata.obsm[use_rep] = adata_rep.obsm[use_rep]


def compute_umap(adata, params):
    try:
        logging.info('Compute UMAP...')
        sc.tl.umap(adata, method='rapids', **params)
    except Exception as e:
        print('sc.tl.umap: Rapids failed, defaulting to UMAP implementation')
        print(e)
        sc.tl.umap(adata, **params)


input_file = snakemake.input[0]
input_rep = snakemake.input.rep
output_file = snakemake.output[0]
params = dict(snakemake.params.items())

logging.info(f'Read {input_file}...')
adata = read_anndata(input_file)

if adata.n_obs == 0:
    logging.info('No data, write empty file...')
    adata.write_zarr(output_file)
    exit(0)

neighbors_key = params.get('neighbors_key', 'neighbors')

# compute UMAP
if isinstance(neighbors_key, list):
    for key in neighbors_key:
        check_and_update_neighbors_info(adata, key, params)
        params |= dict(neighbors_key=key)
        compute_umap(adata, params)
        adata.obsm[f'X_umap_{key}'] = adata.obsm['X_umap']
else:
    check_and_update_neighbors_info(adata, neighbors_key, params)
    params |= dict(neighbors_key=neighbors_key)
    compute_umap(adata, params)

logging.info(f'Write to {output_file}...')
del adata.raw
del adata.X
del adata.layers
del adata.obsp
adata.write_zarr(output_file)

if input_file.endswith('.zarr'):
    input_files = [f.name for f in Path(input_file).iterdir()]
    files_to_keep = [f for f in input_files if f not in ['obsm', 'uns']]
    link_zarr(
        in_dir=input_file,
        out_dir=output_file,
        file_names=files_to_keep,
        overwrite=True,
)