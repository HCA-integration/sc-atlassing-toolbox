import logging
import warnings

import patient_representation as pr
import pandas as pd
import scanpy as sc

from utils.io import read_anndata, write_zarr_linked


warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)

sc.set_figure_params(dpi=100, frameon=False)
input_zarr = snakemake.input.zarr
prepare_zarr = snakemake.input.prepare
output_zarr = snakemake.output.zarr
sample_key = snakemake.params.get('sample_key')
cell_type_key = snakemake.params.get('cell_type_key')
use_rep = snakemake.params.get('use_rep')
n_epochs = snakemake.params.get('n_epochs')


logging.info(f'Read "{input_zarr}"...')
n_obs = read_anndata(input_zarr, obs='obs').n_obs
dask = n_obs > 2e6
adata = read_anndata(
    input_zarr,
    X=use_rep,
    obs='obs',
    var='var',
    backed=dask,
    dask=dask,
    stride=int(n_obs / 5),
)

adata.X = adata.X.astype('int32')
adata.obs[cell_type_key] = adata.obs[cell_type_key].astype(str).fillna('NA').astype('category')

logging.info(f'Calculating scPoli representation for "{cell_type_key}", using cell features from "{use_rep}"')
representation_method = pr.tl.SCPoli(
    sample_key=sample_key,
    cells_type_key=cell_type_key,
    layer='X',
    n_epochs=n_epochs,
    unknown_ct_names=['NA'],
)
representation_method.prepare_anndata(adata)

# compute distances
distances = representation_method.calculate_distance_matrix(force=True)

# create new AnnData object for patient representations
adata = sc.AnnData(
    obs=pd.DataFrame(index=representation_method.samples),
    obsm={'distances': distances}
)

# compute kNN graph
sc.pp.neighbors(
    adata,
    use_rep='distances',
    metric='precomputed'
)

logging.info(f'Write "{output_zarr}"...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=prepare_zarr,
    out_dir=output_zarr,
    files_to_keep=['obsm', 'obsp', 'uns']
)
