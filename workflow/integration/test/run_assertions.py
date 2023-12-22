import glob
import warnings
from pprint import pprint
import numpy as np
from anndata.experimental import read_elem
import zarr
import logging
logging.basicConfig(level=logging.INFO)


# direct integration outputs
outputs = glob.glob('test/out/integration/dataset~*/file_id~*/batch~*/method~*/adata.zarr')
if len(outputs) == 0:
    warnings.warn('No integration outputs found')

for file in outputs:
    logging.info(f'Checking {file}...')
    z = zarr.open(file)
    uns = read_elem(z["uns"])

    try:
        # Check Metadata
        assert 'dataset' in uns
        assert 'methods' in uns
        assert 'integration' in uns
        for key in ['method', 'label_key', 'batch_key', 'output_type']:
            assert key in uns['integration']

        # Unintegrated data
        # adata_raw = raw.to_adata()
        # assert 'X_pca' in adata_raw.obsm
        # assert 'connectivities' in adata_raw.obsp
        # assert 'distances' in adata_raw.obsp

        # Output type specific outputs
        output_types = uns['integration']['output_type']
        output_types = [output_types] if isinstance(output_types, str) else output_types

        if uns['integration']['method'] == 'scanorama':
            assert len(output_types) == 2

        for ot in output_types:
            if ot not in ['knn', 'embed', 'full']:
                raise ValueError(f'Invalid output type {ot}')

        if 'knn' in output_types:
            obsp = read_elem(z["obsp"])
            assert 'connectivities' in obsp
            assert 'distances' in obsp
        if 'embed' in output_types:
            obsm = read_elem(z["obsm"])
            assert 'X_emb' in obsm
        if 'full' in output_types:
            X = read_elem(z['X'])
            assert X is not None
            # assert 'corrected_counts' in layers
            # check that counts are different from input
            # if uns['integration']['method'] != 'unintegrated':
                # assert adata.X.data.shape != adata_raw.X.data.shape or np.any(X.data != adata_raw.X.data)

    except AssertionError as e:
        print('AssertionError for:', file)
        pprint(uns)
        raise e
