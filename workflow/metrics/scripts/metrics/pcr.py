import numpy as np
import pandas as pd
from scipy.sparse import issparse


def pcr(adata, output_type, meta):
    import scib

    if output_type == 'knn':
        return np.nan

    adata_raw = adata.raw.to_adata()

    return scib.me.pcr_comparison(
        adata_pre=adata_raw,
        adata_post=adata,
        covariate=meta['batch'],
        embed='X_emb' if output_type == 'embed' else 'X_pca',
    )


def pcr_y(adata, output_type, meta):
    import scib_metrics

    if output_type == 'knn':
        return np.nan

    adata_raw = adata.raw.to_adata()
    X_pre = adata_raw.X

    if output_type == 'embed':
        X_post = adata.obsm['X_emb']
    else:
        X_post = adata.X

    X_pre, X_post = [X if isinstance(X, np.ndarray) else X.todense() for X in [X_pre, X_post]]

    return scib_metrics.pcr_comparison(
        X_pre=X_pre,
        X_post=X_post,
        covariate=adata.obs[meta['batch']],
        categorical=True
    )


def cell_cycle(adata, output_type, meta):
    import scib

    if output_type == 'knn':
        return np.nan

    adata_raw = adata.raw.to_adata()
    print(adata_raw)

    if 'feature_name' in adata_raw.var.columns:
        adata_raw.var_names = adata_raw.var['feature_name']
    upper_case_genes = sum(adata.var_names.str.isupper())
    organism = 'mouse' if upper_case_genes <= 0.1 * adata.n_vars else 'human'

    # compute cell cycle score per batch
    batch_key = meta['batch']
    for batch in adata.obs[batch_key].unique():
        scib.pp.score_cell_cycle(
            adata_raw[adata_raw.obs[batch_key] == batch],
            organism=organism
        )

    return scib.me.cell_cycle(
        adata_pre=adata_raw,
        adata_post=adata,
        batch_key=meta['batch'],
        embed='X_emb' if output_type == 'embed' else 'X_pca',
        recompute_cc=False,
        organism=organism,
        verbose=False,
    )
