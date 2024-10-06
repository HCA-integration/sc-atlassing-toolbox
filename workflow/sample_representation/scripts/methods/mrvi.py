import logging
import warnings
import torch

import patient_representation as pr
import pandas as pd
import scanpy as sc

from utils.io import read_anndata, write_zarr_linked


warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)

sc.set_figure_params(dpi=100, frameon=False)
input_file = snakemake.input.zarr
prepare_zarr = snakemake.input.prepare
output_file = snakemake.output.zarr
sample_key = snakemake.params.get('sample_key')
cell_type_key = snakemake.params.get('cell_type_key')
use_rep = snakemake.params.get('use_rep')
max_epochs = snakemake.params.get('hyperparams', {}).get('max_epochs')
max_epochs = int(max_epochs) if max_epochs is not None else max_epochs

use_gpu = torch.cuda.is_available()
print(f'GPU available: {use_gpu}', flush=True)

logging.info(f'Read "{input_file}"...')
n_obs = read_anndata(input_file, obs='obs').n_obs
dask = n_obs > 2e6
adata = read_anndata(
    input_file,
    X=use_rep,
    obs='obs',
    backed=dask,
    dask=dask,
    stride=int(n_obs / 5),
)

logging.info(f'Calculating MrVI representation for "{cell_type_key}", using cell features from "{use_rep}"')
representation_method = pr.tl.MrVI(
    sample_key=sample_key,
    cells_type_key=cell_type_key,
    layer='X',
    max_epochs=max_epochs,
    accelerator='gpu' if use_gpu else 'auto',
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

logging.info(f'Write "{output_file}"...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=prepare_zarr,
    out_dir=output_file,
    files_to_keep=['obsm', 'obsp', 'uns']
)
