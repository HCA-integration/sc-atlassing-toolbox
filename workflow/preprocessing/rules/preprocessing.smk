"""
Preprocessing steps
Only the unique outputs per step are saved for storage efficiency. For assembled zarr files, see `assemble.smk`.
"""

rule normalize:
    input:
        lambda w: get_for_dataset(config, w.dataset, ['input', module_name])
    output:
        zarr=directory(out_dir / '{dataset}' / 'normalized.zarr'),
    params:
        raw_counts=lambda w: get_for_dataset(config, w.dataset, [module_name, 'raw_counts']),
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    conda:
        '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/normalize.py'


rule highly_variable_genes:
    input:
        zarr=rules.normalize.output.zarr
    output:
        zarr=directory(out_dir / '{dataset}' / 'highly_variable_genes.zarr')
    params:
        args=lambda w: get_for_dataset(config, w.dataset, [module_name, 'highly_variable_genes']),
        batch=lambda w: get_for_dataset(config, w.dataset, [module_name, 'batch']),
        lineage=lambda w: get_for_dataset(config, w.dataset, [module_name, 'lineage']),
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    conda:
        '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/highly_variable_genes.py'


rule pca:
    input:
        zarr=rules.highly_variable_genes.output.zarr,
        counts=rules.normalize.output.zarr,
    output:
        zarr=directory(out_dir / '{dataset}' / 'pca.zarr')
    params:
        scale=lambda w: get_for_dataset(config, w.dataset, [module_name, 'scale'])
    resources:
        mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb'),
        disk_mb=get_resource(config,profile='cpu_merged',resource_key='disk_mb'),
    conda: '../envs/scanpy.yaml'
    # shadow: 'minimal'
    script:
        '../scripts/pca.py'


rule neighbors:
    input:
        zarr=rules.pca.output.zarr
    output:
        zarr=directory(out_dir / '{dataset}' / 'neighbors.zarr')
    params:
        args=lambda w: get_for_dataset(config, w.dataset, [module_name, 'neighbors'], default={}),
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
    conda:
        ifelse(
            'use_gpu' not in config.keys() or not config['use_gpu'],
            _if='../envs/scanpy.yaml', _else='../envs/scanpy_rapids.yaml'
        )
    # shadow: 'minimal'
    script:
        '../scripts/neighbors.py'


rule umap:
    input:
        zarr=rules.neighbors.output.zarr,
        rep=rules.pca.output.zarr,
    output:
        zarr=directory(out_dir / '{dataset}' / 'umap.zarr')
    params:
        neighbors_key='neighbors',
    resources:
        partition=get_resource(config,profile='gpu',resource_key='partition'),
        qos=get_resource(config,profile='gpu',resource_key='qos'),
        gpu=get_resource(config,profile='gpu',resource_key='gpu'),
        mem_mb=get_resource(config,profile='gpu',resource_key='mem_mb'),
    conda:
        ifelse(
            'use_gpu' not in config.keys() or not config['use_gpu'],
            _if='../envs/scanpy.yaml', _else='../envs/scanpy_rapids.yaml'
        )
    # shadow: 'minimal'
    script:
        '../scripts/umap.py'
