import pandas as pd
from scipy import sparse
import numpy as np
import scib
from anndata import AnnData
import yaml
from concurrent.futures import ProcessPoolExecutor, as_completed
try:
    from sklearnex import patch_sklearn
    patch_sklearn()
except ImportError:
    print('no hardware acceleration for sklearn', flush=True)

from utils.io import read_anndata


input_file = snakemake.input.anndata
setup_file = snakemake.input.setup
output_file = snakemake.output.tsv
covariate = snakemake.wildcards.covariate
sample_key = snakemake.params.get('sample_key')
n_threads = np.max([snakemake.threads, 1])

print('Read anndata file...', flush=True)
adata = read_anndata(input_file, obsm='obsm', obs='obs', uns='uns')
adata = adata[adata.obs[covariate].notna()].copy()
n_covariate = adata.obs[covariate].nunique()

# set default sample key
if sample_key is None or sample_key == 'None':
    print('Using index as sample key...', flush=True)
    sample_key = 'index'
    adata.obs[sample_key] = adata.obs.index

# make sure the PCA embedding is an array
if not isinstance(adata.obsm['X_pca'], np.ndarray):
    adata.obsm['X_pca'] = adata.obsm['X_pca'].toarray()

# read setup file
with open(setup_file, 'r') as f:
    setup = yaml.safe_load(f)
# n_permute = setup['n_permute']
# n_permute = min(snakemake.params.get('n_permute', 0), n_permute)
n_permute = snakemake.params.get('n_permute', 0)
print(f'n_permute: {n_permute}', flush=True)

# Permute covariates
perm_covariates = []

if covariate == sample_key:
    print('Sample key is the same as covariate, skipping permutation...', flush=True)
    n_permute = 0

series_dict = {}
for i in range(n_permute):
    # aggregate and permutate covariate
    cov_per_sample = adata.obs.groupby(sample_key, observed=True).agg({covariate: 'first'})
    cov_map = dict(zip(cov_per_sample.index, cov_per_sample[covariate].sample(frac=1)))
    
    # apply permutation
    covariate_perm = f'{covariate}-{i}'
    series_dict[covariate_perm] = adata.obs[sample_key].map(cov_map)
    perm_covariates.append(covariate_perm)

# concat new columns to adata to avoid fragmentation
adata.obs = pd.concat([adata.obs, pd.DataFrame(series_dict)], axis=1).copy()

# Determine all covariates
covariates = [covariate]+perm_covariates

with ProcessPoolExecutor(max_workers=n_threads) as pool:
    def compute_pcr(adata, covariate, is_permuted, n_threads=1):
        pcr = scib.me.pcr(
            adata,
            covariate=covariate,
            recompute_pca=False,
            verbose=False,
            linreg_method='numpy',
            n_threads=n_threads,
        )
        return (covariate, is_permuted, pcr)

    futures = [
        pool.submit(
            compute_pcr,
            adata,
            covariate=c,
            is_permuted=c in perm_covariates
        )
        for c in covariates
    ]
    pcr_scores = []
    completed_futures = 0
    total_futures = len(futures)
    
    for future in as_completed(futures):
        completed_futures += 1
        pcr_scores.append(future.result())
        print(f'{completed_futures}/{total_futures} completed for {covariate}', flush=True)

# Set permuted score when covariate is the same as the group variable
if covariate == sample_key:
    # permutations wouln't change values in this case, impute same value as covariate score
    perm_score = (
        f'{covariate}-0',
        True,
        pcr_scores[0][2]
    )
    pcr_scores.append(perm_score)

df = pd.DataFrame.from_records(
    pcr_scores,
    columns=['covariate', 'permuted', 'pcr'],
)

# calculate summary stats
df['covariate'] = df['covariate'].str.split('-', expand=True)[0].astype('category')
df['n_covariates'] = n_covariate
df['perm_mean'] = df.loc[df['permuted'], 'pcr'].mean()
df['perm_std'] = df.loc[df['permuted'], 'pcr'].std()
df['z_score'] = (df['pcr'] - df['perm_mean']) / df['perm_std']
df['signif'] = df['z_score'] > 1.5
df['p-val'] = df.loc[df['permuted'], 'signif'].sum() / n_permute

print(df, flush=True)
df.to_csv(output_file, sep='\t', index=False)