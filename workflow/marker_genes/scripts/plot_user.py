from matplotlib import pyplot as plt
import scanpy as sc
import anndata
from pprint import pformat
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata
from utils.misc import dask_compute, ensure_dense


sc.set_figure_params(frameon=False)
plt.rcParams['figure.figsize'] = 30, 25
input_file = snakemake.input[0]
output_png = snakemake.output.dotplot
markers = snakemake.params.markers

wildcards = snakemake.wildcards
group_key = wildcards.group
file_id = wildcards.file_id
title = f'Marker genes for file_id={file_id} group={group_key}'

logging.info(f'Reading {input_file}...')
adata = read_anndata(
    input_file,
    X='X',
    obs='obs',
    var='var',
    uns='uns',
    dask=True,
    backed=True,
)

if 'feature_name' in adata.var.columns:
    adata.var_names = adata.var['feature_name']

# match marker genes and var_names
logging.info(pformat({k: len(v) for k, v in markers.items()}))
markers = {
    k: adata.var[adata.var_names.isin(v)].index.to_list()
    for k, v in markers.items()
}
logging.info(pformat({k: len(v) for k, v in markers.items()}))

# subset to genes for plotting
genes_to_plot = list(set([element for sublist in markers.values() for element in sublist]))
logging.info(f'Load {len(genes_to_plot)} genes to memory...')
adata = adata[:, adata.var_names.isin(genes_to_plot)].copy()
adata = dask_compute(adata)
adata = ensure_dense(adata)
logging.info(adata.__str__())

# if no cells or genes to plot, save empty plots
if min(adata.shape) == 0:
    plt.savefig(output_png)
    exit()

sc.pl.dotplot(
    adata,
    markers,
    groupby=group_key,
    use_raw=False,
    standard_scale='var',
    # title=title,
    show=False,
    swap_axes=False,
)

plt.savefig(output_png, dpi=100, bbox_inches='tight')
