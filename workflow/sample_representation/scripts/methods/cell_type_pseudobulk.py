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

logging.info(f'Read "{input_zarr}"...')
adata = read_anndata(
    input_zarr,
    obs='obs',
    X=f'obsm/{use_rep}',
)
adata.obsm[use_rep] = adata.X

logging.info(f'Calculating Cell Type Pseudobulk representation for "{cell_type_key}", using cell features from "{use_rep}"')
representation_method = pr.tl.CellTypePseudobulk(
    sample_key=sample_key,
    cells_type_key=cell_type_key,
    layer=use_rep
)
representation_method.prepare_anndata(adata)  # Assuming that small cell types and samples are already filtered out

# compute distances
distances = representation_method.calculate_distance_matrix(
    force=True,
    aggregate="mean",
    dist="euclidean"
)

# create new AnnData object for patient representations
adata = sc.AnnData(
    obs=pd.DataFrame(index=representation_method.samples),
    obsm={'distances': distances} | {
        f'{cell_type}_pseudobulk': representation_method.patient_representations[i]
        for i, cell_type in enumerate(representation_method.cell_types)
    }
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
