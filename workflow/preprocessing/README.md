# Preprocessing

This module runs a configurable single-cell preprocessing workflow and assembles selected results into one output zarr file.

Implemented steps:

1. Normalize (`normalize.py`)
2. Mark non-zero genes for HVG-safe filtering (`filter_genes.py`)
3. Highly variable genes (`highly_variable_genes.py`)
4. Optional extra HVGs (`extra_hvgs.py`)
5. PCA (`pca.py`)
6. Neighbors graph (`neighbors.py`)
7. UMAP (`umap.py`)
8. PCA/UMAP plots (`plot.py`)
9. Assembly (`assemble.py`)

Rules are declared in `rules/rules.smk`, parameterized in `rules/assemble.smk`, and plotting rules are in `rules/plots.smk`.

## Testing

You can run the module on a small test dataset using:

```bash
snakemake --configfile test/config_no_gpu.yaml --use-conda --cores 8

# or, for the GPU-enabled variant
snakemake --configfile test/config.yaml --use-conda --cores 8
```

## Quickstart configuration

Global settings/defaults plus the `no_hvg` scenario from [`test/config.yaml`](test/config.yaml):

```{eval-rst}
.. literalinclude:: ../../../workflow/preprocessing/test/config.yaml
   :language: yaml
   :lines: 1-19,96-108
```

## Further reading

For the full parameter reference — global configuration knobs, per-step behavior and arguments (Normalize, Filter genes, HVG, Extra HVGs, PCA, Neighbors, UMAP, Plots), assembly semantics, and other notes — see [`docs/modules/preprocessing/preprocessing_reference.md`](preprocessing_reference.md).
