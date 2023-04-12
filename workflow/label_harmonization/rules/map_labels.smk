rule celltypist:
    """
    Map author labels of datasets
    """
    input:
        anndata=get_input,
    output:
        h5ad=out_dir / 'celltypist' / '{dataset}' / 'adata.h5ad',
        reannotation=out_dir / 'celltypist' / '{dataset}' / 'reannotation.tsv',
        relation=out_dir / 'celltypist' / '{dataset}' / 'relation.tsv',
        model=out_dir / 'celltypist' / '{dataset}' / 'model.pkl',
    params:
        author_label_key=lambda w: get_input(w,query=[module_name,'author_label_key']),
        dataset_key=lambda w: get_input(w,query=[module_name,'dataset_key']),
        params=lambda w: get_input(w,query=[module_name,'celltypist']),
    conda:
        '../envs/celltypist.yaml'
    resources:
        partition=lambda w: get_resource(config,resource_key='partition'),
        qos=lambda w: get_resource(config,resource_key='qos'),
        mem_mb=lambda w: get_resource(config,resource_key='mem_mb'),
    script:
        '../scripts/celltypist.py'


rule celltypist_plots:
    """
    Plots for celltypist output
    """
    input:
        model=rules.celltypist.output.model,
    output:
        treeplot=out_dir / 'celltypist' / '{dataset}' / 'treeplot.pdf',
        heatmap=out_dir / 'celltypist' / '{dataset}' / 'heatmap.pdf',
        # sankeyplot=out_dir / 'celltypist' / '{dataset}_sankeyplot.pdf',
    conda:
        '../envs/celltypist.yaml'
    script:
        '../scripts/celltypist_plots.py'


celltypist_columns = ['dataset', 'dataset_key', 'author_label_key', 'celltypist']
try:
    celltypist_datasets = parameters[celltypist_columns].dropna()['dataset'].unique()
except KeyError:
    celltypist_datasets = []

rule celltypist_all:
    input:
        expand(rules.celltypist_plots.output, dataset=celltypist_datasets)
