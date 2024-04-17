import warnings
warnings.filterwarnings("ignore")
import logging
logging.basicConfig(level=logging.INFO)
import numpy as np
from utils.processing import sc

from utils.io import read_anndata

input_file = snakemake.input[0]
output_file = snakemake.output[0]
resolution = float(snakemake.wildcards.resolution)
neighbors_key = snakemake.params.get('neighbors_key', 'neighbors')
cluster_alg = snakemake.params.get('algorithm', 'louvain')

logging.info(f'Read anndata file {input_file}...')
adata = read_anndata(input_file, obs='obs', uns='uns', obsp='obsp')

if neighbors_key not in adata.uns:
    assert 'connectivities' in adata.obsp
    assert 'distances' in adata.obsp
    adata.uns[neighbors_key] = {
        'connectivities_key': 'connectivities',
        'distances_key': 'distances',
    }

neighbors = adata.uns[neighbors_key]
adata.uns['neighbors'] = neighbors
adata.obsp['connectivities'] = adata.obsp[neighbors['connectivities_key']]
adata.obsp['distances'] = adata.obsp[neighbors['distances_key']]

logging.info(f'{cluster_alg} clustering with resolution {resolution}...')
cluster_key = f'{cluster_alg}_{resolution}'

cluster_alg_map = {
    'louvain': sc.tl.louvain,
    'leiden': sc.tl.leiden,
}

cluster_alg_map[cluster_alg](
    adata,
    resolution=resolution,
    key_added=cluster_key,
)

logging.info('Write file...')
adata.obs[cluster_key].to_csv(output_file, sep='\t', index=True)