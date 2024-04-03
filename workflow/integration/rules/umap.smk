######### UMAP and embedding plots #########

use rule umap from preprocessing as integration_compute_umap with:
    input:
        anndata=rules.integration_postprocess.output.zarr,
        rep=rules.integration_postprocess.output.zarr,
    output:
        zarr=directory(out_dir / 'umap' / f'{paramspace.wildcard_pattern}.zarr'),
        done=touch(out_dir / 'umap' / f'{paramspace.wildcard_pattern}.done'),
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),


def get_colors(wildcards):
    dataset = wildcards.dataset
    labels = mcfg.get_from_parameters(wildcards, 'label')
    labels = labels if isinstance(labels, list) else [labels]
    batch = mcfg.get_from_parameters(wildcards, 'batch')
    batch = batch if isinstance(batch, list) else [batch]
    umap_colors = mcfg.get_for_dataset(dataset, query=[mcfg.module_name, 'umap_colors'], default=[])
    return [*labels, *batch, *umap_colors]


use rule plots from preprocessing as integration_plot_umap with:
    input:
        anndata=rules.integration_compute_umap.output.zarr,
    output:
        plots=directory(image_dir / 'umap' / paramspace.wildcard_pattern.replace('--output_type', '/output_type')),
    params:
        color=get_colors,
        basis='X_umap',
        ncols=1,
        outlier_factor=10,
        size=2,
    resources:
        partition=mcfg.get_resource(profile='cpu',resource_key='partition'),
        qos=mcfg.get_resource(profile='cpu',resource_key='qos'),
        mem_mb=mcfg.get_resource(profile='cpu',resource_key='mem_mb'),
        gpu=mcfg.get_resource(profile='cpu',resource_key='gpu'),
