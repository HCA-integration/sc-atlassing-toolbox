import logging
logging.basicConfig(level=logging.INFO)
from anndata.experimental import read_elem
import zarr
from pprint import pformat
from pathlib import Path


input_file = snakemake.input.anndata
output_dir = snakemake.output[0]
sample_key = snakemake.params.get('sample_key')
Path(output_dir).mkdir(parents=True)

covariates = snakemake.params.get('covariates', [])
perm_covariates = snakemake.params.get('permute_covariates')
if perm_covariates is None:
    perm_covariates = covariates
n_perms = snakemake.params.get('n_permute')

logging.info(f'Read {input_file}...')
if input_file.endswith('.h5ad'):
    obs = read_anndata(input_file).obs
else:
    z = zarr.open(input_file)
    obs = read_elem(z['obs'])


def covariate_valid(obs, covariate):
    return (covariate in obs.columns) and (obs[covariate][obs[covariate].notna()].nunique() >= 2)


covariates = [
    covariate for covariate in covariates
    if covariate_valid(obs, covariate)
]
perm_covariates = [
    covariate for covariate in perm_covariates
    if covariate_valid(obs, covariate) and covariate != sample_key
]

logging.info(f'covariates:\n {pformat(covariates)}')
logging.info(f'perm_covariates:\n {pformat(perm_covariates)}')
logging.info(f'n_permutations:\n {pformat(n_perms)}')


logging.info(f'Write covariate setups to {output_dir}...')
for covariate in covariates + perm_covariates:
    covariate_file = f'{output_dir}/{covariate}.yaml'
    logging.info(covariate_file)
    with open(covariate_file, 'w') as f:
        if covariate in perm_covariates:
            f.write(f'n_permute: {n_perms}')
        else:
            f.write('n_permute: 0')
