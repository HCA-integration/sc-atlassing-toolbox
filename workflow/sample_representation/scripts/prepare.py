import numpy as np
from matplotlib import pyplot as plt
import anndata
import pandas as pd
import scanpy as sc
import warnings
warnings.simplefilter("ignore", UserWarning)
import logging
logging.basicConfig(level=logging.INFO)

import patient_representation as pr

from utils.io import read_anndata
from utils.processing import get_pseudobulks


sc.set_figure_params(dpi=100, frameon=False)
input_zarr = snakemake.input.zarr
output_zarr = snakemake.output.zarr
sample_key = snakemake.params.get('sample_key')
cell_type_key = snakemake.params.get('cell_type_key')
aggregate = snakemake.params.get('aggregate')
layer = snakemake.params.get('norm_counts')
min_cells_per_sample = snakemake.params.get('min_cells_per_sample')
min_cells_per_cell_type = snakemake.params.get('min_cells_per_cell_type')

logging.info(f'Read "{input_zarr}"...')
n_obs = read_anndata(input_zarr, obs='obs').n_obs
dask = n_obs > 2e6
adata = read_anndata(
    input_zarr,
    X=layer,
    obs='obs',
    var='var',
    backed=dask,
    dask=dask,
    stride=int(n_obs / 5),
)

# # filter small samples and cell types
# adata = pr.pp.filter_small_samples(
#     adata,
#     sample_key,
#     sample_size_threshold=min_cells_per_sample,
# )
# if cell_type_key is not None:
#     adata = pr.pp.filter_small_cell_types(
#         adata,
#         sample_key,
#         cell_type_key,
#         cluster_size_threshold=min_cells_per_cell_type,
#     )

logging.info(f'Pseudobulk by "{sample_key}"...')
adata_bulk = get_pseudobulks(adata, group_key=sample_key, agg=aggregate)

logging.info(f'Write "{output_zarr}"...')
logging.info(adata_bulk.__str__())
adata_bulk.write_zarr(output_zarr)

# # process pseudobulk adata
# sc.pp.log1p(adata_bulk)
# sc.pp.highly_variable_genes(adata_bulk, n_top_genes=2000)
# sc.pp.pca(adata_bulk)

# # scree plot TODO: move to preprocessing plots
# plt.plot(adata_bulk.uns['pca']['variance_ratio'], marker='.')
# plt.title('Scree plot: PC variance contribution')
# plt.savefig(output_pca_scree)

# # Check if it's enough for at least the first two componenets
# if adata.obs[bulk_by].nunique() < 3:
#     plt.savefig(output_pca_1_2)
#     plt.savefig(output_pca_2_3)
#     exit(0)

# # remove empty columns
# colors = [c for c in colors if not adata_bulk.obs[c].isna().all()]
# # remove colors if zero or too many entries
# colors = [c for c in colors if 0 < adata_bulk.obs[c].nunique() <= 64]
# # plot without colors, if no colors are available
# if len(colors) == 0:
#     colors = None

# n_rows = len(colors)
# height = np.max([7, n_rows * 0.2 * 7])
# plt.rcParams['figure.figsize'] = (10, height)

# sc.pl.pca(
#     adata_bulk[adata_bulk.obs.sample(adata_bulk.n_obs).index],
#     components='1,2',
#     color=colors,
#     title=[f'PCA 1-2 {dataset}: Pseudobulk={bulk_by} color={c}' for c in colors],
#     show=False,
#     ncols=1
# )
# plt.tight_layout()
# plt.savefig(output_pca_1_2, bbox_inches='tight')

# # missing third component
# if len(adata_bulk.uns['pca']['variance_ratio']) < 3:
#     plt.savefig(output_pca_2_3)
#     exit(0)

# sc.pl.pca(
#     adata_bulk[adata_bulk.obs.sample(adata_bulk.n_obs).index],
#     components='2,3',
#     color=colors,
#     title=[f'PCA 2-3{dataset}: Pseudobulk={bulk_by} color={c}' for c in colors],
#     show=False,
#     ncols=1
# )
# plt.tight_layout()
# plt.savefig(output_pca_2_3, bbox_inches='tight')
