from pathlib import Path
import pandas as pd

from utils.io import read_anndata, link_zarr

input_zarr = snakemake.input.zarr
output_dir = snakemake.output[0]
batch_key = snakemake.params.get('batch_key')

adata = read_anndata(input_zarr, obs='obs')

output_dir = Path(output_dir)
output_dir.mkdir(exist_ok=True, parents=True)

if batch_key not in adata.obs.columns:
    open(output_dir / 'no_batch.txt', 'w').close()
else:
    value_counts = adata.obs[batch_key].value_counts()
    for batch in adata.obs[batch_key].unique():
        with open(output_dir / f'{batch}.txt', 'w') as f:
            f.write(str(value_counts.loc[batch]))
