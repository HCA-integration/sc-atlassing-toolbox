import logging
from pathlib import Path
import gc

import pandas as pd
import scanpy as sc
import anndata

from utils import SCHEMAS


logging.basicConfig(level=logging.INFO)

def read_adata(file):
    ad = anndata.read_zarr(file)
    ad.var = ad.var[SCHEMAS['CELLxGENE_VARS']]
    # remove data
    del ad.uns
    del ad.layers
    del ad.raw
    del ad.obsm
    return ad


dataset = snakemake.params.dataset
files = snakemake.input
out_file = snakemake.output.zarr
merge_strategy = snakemake.params.merge_strategy

if len(files) == 1:
    in_file = Path(files[0])
    out_dir = Path(out_file)
    if not out_dir.exists():
        out_dir.mkdir()
    for f in in_file.iterdir():
        if f.name == '.snakemake_timestamp':
            continue  # skip snakemake timestamp
        new_file = out_dir / f.name
        new_file.symlink_to(f.resolve())
else:
    logging.info(f'Read first file {files[0]}...')
    adata = read_adata(files[0])
    print(adata)

    for file in files[1:]:
        logging.info(f'Read {file}...')
        _adata = read_adata(file)
        print(_adata)

        if _adata.n_obs == 0:
            logging.info('skip concatenation...')
            continue

        logging.info('Concatenate...')
        # merge genes
        var_map = pd.merge(
            adata.var,
            _adata.var,
            how=merge_strategy,
            on=['feature_id'] + SCHEMAS['CELLxGENE_VARS']
        )
        var_map = var_map[~var_map.index.duplicated()]

        # merge adata
        adata = sc.concat([adata, _adata], join='outer')
        adata = adata[:, var_map.index]
        adata.var = var_map.loc[adata.var_names]

        del _adata
        gc.collect()

    organ = adata.obs['organ'].unique()
    assert not len(organ) > 1
    adata.uns['dataset'] = dataset
    adata.uns['organ'] = organ
    adata.uns['meta'] = {'dataset': dataset, 'organ': organ}
    print(adata)
    print(adata.var)

    assert 'feature_name' in adata.var

    logging.info('Write...')
    adata.write_zarr(out_file)
