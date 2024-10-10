import numpy as np

def by_sample(
    adata,
    n_cell_max,
    sample_key,
    random_state=42,
    min_cells_per_sample=100,
    **kwargs
):
    """
    Randomly subset complete samples until the maximum number of cells is reached
    """
    samples = []
    n_cells = 0

    shuffled_samples = adata.obs[sample_key].value_counts().sample(frac=1, random_state=random_state)
    
    for sample, count in shuffled_samples.items():
        if len(samples) > 0 and n_cells > n_cell_max:
            break
        if count < min_cells_per_sample:
            # skip samples that are too small
            continue
        n_cells += count
        samples.append(sample)

    return adata.obs[sample_key].isin(samples)


def within_sample(
    adata,
    n_cell_max,
    sample_key,
    n_cells_per_sample=None,
    random_state=42,
    **kwargs
):
    """
    Subset to n_cell_max / n_samples random cells per sample
    """
    if n_cell_max is None:
        n_cell_max = adata.n_obs
    
    if n_cells_per_sample is None:
        n_samples = adata.obs[sample_key].nunique()
        n_cells_per_sample = int(n_cell_max / n_samples)
    
    shuffled_samples = adata.obs[sample_key].sample(frac=1, random_state=random_state)
    subset_indices = []
    n_cells_flagged = 0
    for sample in shuffled_samples:
        # subset to sample
        sample_mask = adata.obs[sample_key] == sample
        # turn mask into subset index
        sample_indices = adata.obs.loc[sample_mask].index
        
        # subset to n_cells_per_sample random cells
        n_cells = min(len(sample_indices), n_cells_per_sample)
        sampled_indices = sample_indices.to_series().sample(n_cells, random_state=random_state)

        subset_indices.append(sampled_indices)
        n_cells_flagged += len(sampled_indices)
        
        # exit loop if n_cell_max is reached
        if n_cells_flagged >= n_cell_max:
            break
    
    adata.obs['subset'] = False
    adata.obs.loc[np.concatenate(subset_indices), 'subset'] = True
    return adata.obs['subset']


def scarf_TopACeDo(adata, n_cell_max, sample_key):
    pass


SUBSET_MAP = {
    'by_sample': by_sample,
    'within_sample': within_sample,
    'scarf_TopACeDo': scarf_TopACeDo,
}