import glob
import warnings
from pprint import pformat
import numpy as np
from anndata.experimental import read_elem
from anndata import AnnData
import zarr
import logging
logging.basicConfig(level=logging.INFO)


# direct integration outputs
adata_files = glob.glob('test/out/integration/dataset~*/file_id~*/batch~*/method~scvi*/adata.zarr')
model_files = glob.glob('test/out/integration/dataset~*/file_id~*/batch~*/method~scvi*/model')
if len(adata_files) == 0:
    warnings.warn('No integration outputs found')
    exit()

import scvi
import torch

for adata_file, model_file in zip(adata_files, model_files):
    logging.info(f'Checking {adata_file}...')
    with zarr.open(adata_file) as z:
        adata = AnnData(
            X=read_elem(z["X"]),
            obs=read_elem(z["obs"]),
            var=read_elem(z["var"]),
            uns=read_elem(z["uns"])
        )

    try:
        assert 'model_history' in adata.uns
        model = scvi.model.SCVI.load(model_file, adata)
        logging.info(pformat(model))
        model = pytorch.load(f'{model_file}/model.pt')
        logging.info(pformat(model))
    except AssertionError as e:
        logging.info('AssertionError for:', adata_file)
        logging.info(pformat(adata.uns))
        raise e
