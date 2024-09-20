import anndata as ad
import scanpy as sc

from utils.accessors import subset_hvg


def morans_i(
    adata,
    output_type,
    batch_key,
    label_key,
    adata_raw,
    var_key='metrics_features',
    n_threads=1,
    # covariate_key='disease',
    **kwargs
):
    def score_per_batch(adata, batch_key, genes, score_name='score', n_threads=1):
        from tqdm import tqdm
        from concurrent.futures import ThreadPoolExecutor, as_completed

        def subset_and_score(adata, batch_key, batch, genes, score_name):
            ad_sub = adata[adata.obs[batch_key] == batch].copy()
            sc.tl.score_genes(
                ad_sub,
                gene_list=genes,
                score_name=score_name,
                copy=False,
            )
            scores = ad_sub.obs[score_name].copy()
            del ad_sub
            return scores
        
        batches = adata.obs[batch_key].unique()
        
        if n_threads == 1:
            for batch in tqdm(batches):
                scores_series = subset_and_score(
                    adata=adata,
                    batch_key=batch_key,
                    batch=batch,
                    genes=genes,
                    score_name=score_name,
                )
                adata.obs.loc[scores_series.index, score_name] = scores_series
        else:
            scores_list = []
            with ThreadPoolExecutor(max_workers=n_threads) as executor, \
                 tqdm(total=len(batches)) as pbar:
                futures = [
                    executor.submit(
                        subset_and_score,
                        adata=adata,
                        batch_key=batch_key,
                        batch=batch,
                        genes=genes,
                        score_name=score_name,
                    )
                    for batch in batches
                ]
                for future in as_completed(futures):
                    scores_series = future.result()
                    scores_list.append(scores_series)
                    pbar.update(1)
            for scores_series in scores_list:
                adata.obs.loc[scores_series.index, score_name] = scores_series
        assert score_name in adata.obs.columns, f'No score found in adata.obs.columns {adata.obs.columns[:10]}'
        return adata


    ifn_yoshida = [
        "IRF7",
        "XAF1",
        "UBE2L6",
        "TRIM22",
        "STAT1",
        "SP110",
        "SAMD9L",
        "SAMD9",
        "PLSCR1",
        "PARP9",
        "OAS2",
        "OAS1",
        "MX2",
        "MX1",
        "LY6E",
        "ISG15",
        "IFIT3",
        "IFI6",
        "IFI44L",
        "IFI35",
        "HERC5",
        "EPSTI1",
        "EIF2AK2",
        "CMPK2",
        "BST2"
    ]
    
    adata = ad.AnnData(
        X=adata_raw.X,
        var=adata_raw.var,
        obs=adata.obs,
        obsp=adata.obsp,
        uns=adata.uns,
    )
    
    if 'feature_name' in adata.var.columns:
        adata.var_names = adata.var['feature_name']

    ifn_yoshida = [g for g in ifn_yoshida if g in adata.var_names]
    assert len(ifn_yoshida) > 0, f'No genes found in adata.var_names {adata.var_names[:10]}'
    
    # subset to HVGs and gene score genes
    adata.var[var_key] = adata.var[var_key] | adata.var_names.isin(ifn_yoshida)
    adata, _ = subset_hvg(
        adata,
        var_column=var_key,
        min_cells=1,
        compute_dask=True,
    )
    
    print(f'Score genes per batch with {n_threads} threads...', flush=True)
    score_per_batch(
        adata,
        batch_key=batch_key,
        genes=ifn_yoshida,
        score_name='score',
        n_threads=n_threads,
        # n_threads=1,
    )
    
    return sc.metrics.morans_i(adata, vals=adata.obs['score'])
