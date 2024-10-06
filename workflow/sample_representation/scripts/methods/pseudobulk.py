import logging
import warnings

import scanpy as sc

from utils.io import read_anndata, write_zarr_linked

warnings.simplefilter("ignore", UserWarning)
logging.basicConfig(level=logging.INFO)

sc.set_figure_params(dpi=100, frameon=False)
input_file = snakemake.input.zarr
prepare_zarr = snakemake.input.prepare
output_file = snakemake.output.zarr

logging.info(f'Read "{input_file}"...')
adata = read_anndata(
    prepare_zarr,
    X='X',
    obs='obs',
    var='var',
)

# compute kNN graph
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.pca(adata)
sc.pp.neighbors(adata, use_rep='X_pca')

logging.info(f'Write "{output_file}"...')
logging.info(adata.__str__())
write_zarr_linked(
    adata,
    in_dir=prepare_zarr,
    out_dir=output_file,
    files_to_keep=['obsm', 'obsp', 'uns']
)
