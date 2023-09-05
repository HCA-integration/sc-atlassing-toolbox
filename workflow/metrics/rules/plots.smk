module plots:
    snakefile: "../../common/rules/plots.smk"
    config: config


# barplots

use rule barplot from plots as metrics_barplot with:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=image_dir / 'all' /'{metric}-barplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        facet_row='dataset',
        #facet_col='output_type',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,
    group:
        'metrics_plots'


use rule barplot from plots as metrics_barplot_per_dataset with:
    input:
        tsv=rules.merge_per_dataset.output.tsv
    output:
        png=image_dir / 'per_dataset' / '{dataset}' / '{metric}-barplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        #facet_col='output_type',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,
    group:
        'metrics_plots'


use rule barplot from plots as metrics_barplot_per_method with:
    input:
        tsv=rules.merge_per_method.output.tsv
    output:
        png=image_dir / 'per_method' / '{method}' / '{metric}-barplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='dataset',
        #facet_col='output_type',
        title='Metrics computation',
        description=wildcards_to_str,
        dodge=True,
    group:
        'metrics_plots'


# swarm plots

use rule swarmplot from plots as metrics_swarmplot with:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=image_dir / 'all' / '{metric}-swarmplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        facet_col='dataset',
        facet_row='output_type',
        title='Metrics',
        description=wildcards_to_str,
    group:
        'metrics_plots'


use rule swarmplot from plots as metrics_swarmplot_per_dataset with:
    input:
        tsv=rules.merge_per_dataset.output.tsv
    output:
        png=image_dir / 'per_dataset' / '{dataset}' / '{metric}-swarmplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='method',
        facet_row='output_type',
        title='Metrics',
        description=wildcards_to_str,
        dodge=True,
    group:
        'metrics_plots'

use rule swarmplot from plots as metrics_swarmplot_per_method with:
    input:
        tsv=rules.merge_per_method.output.tsv
    output:
        png=image_dir / 'per_method' / '{method}' / '{metric}-swarmplot.png'
    params:
        metric=lambda wildcards: wildcards.metric,
        category='metric',
        hue='dataset',
        facet_row='output_type',
        title='Metrics',
        description=wildcards_to_str,
        dodge=True,
    group:
        'metrics_plots'

rule compare_metrics:
    input:
        tsv=rules.merge.output.tsv
    output:
        time=image_dir / 'comparison_time.png',
        score=image_dir / 'comparison_score.png',
    conda:
        get_env(config, 'plots')
    group:
        'metrics_plots'
    script:
        '../scripts/plots/comparison.py'


# Funkyheatmap

rule funkyheatmap:
    input:
        tsv=rules.merge.output.tsv
    output:
        pdf=image_dir / 'all' / 'funky_heatmap.pdf',
        tsv=image_dir / 'all' / 'funky_heatmap.tsv'
    params:
        id_vars=['dataset', 'method', 'output_type', 'batch', 'label', 'lineage_specific', 'lineage_key', 'lineage'], # TODO: 'hyperparams'
        variable_var='metric',
        value_var='score',
        weight_batch=0.4,
        n_top=50,
    conda:
        get_env(config, 'funkyheatmap')
    singularity:
        'docker://ghcr.io/dynverse/funky_heatmap:latest'
    group:
        'metrics_plots'
    script:
        '../scripts/plots/funkyheatmap.R'


use rule funkyheatmap as funkyheatmap_per_dataset with:
    input:
        tsv=rules.merge_per_dataset.output.tsv
    output:
        pdf=image_dir / 'per_dataset' / '{dataset}' / 'funky_heatmap.pdf',
        tsv=image_dir / 'per_dataset' / '{dataset}' / 'funky_heatmap.tsv',
    params:
        id_vars=['method', 'output_type', 'batch', 'label', 'lineage_specific', 'lineage_key', 'lineage'],
        variable_var='metric',
        value_var='score',
        weight_batch=0.4,
        n_top=50,


use rule funkyheatmap as funkyheatmap_per_method with:
    input:
        tsv=rules.merge_per_method.output.tsv
    output:
        pdf=image_dir / 'per_method' / '{method}' / 'funky_heatmap.pdf',
        tsv=image_dir / 'per_method' / '{method}' / 'funky_heatmap.tsv',
    params:
        id_vars=['method', 'output_type', 'batch', 'label', 'lineage_specific', 'lineage_key', 'lineage'],
        variable_var='metric',
        value_var='score',
        weight_batch=0.4,
        n_top=50,


rule funkyheatmap_standalone:
    input:
        tsv=rules.merge.output.tsv
    output:
        png=image_dir / 'funky_heatmap.png'
    group:
        'metrics_plots'
    shell:
        """
        wget https://github.com/dynverse/funkyheatmap/releases/latest/download/executable.zip
        unzip -o executable.zip -d funky_heatmap
        ./funky_heatmap/funky_heatmap --data {input.tsv} --output {output.png}
        """


rule plots_all:
    input:
        # funky heatmap
        rules.funkyheatmap.output,
        expand(rules.funkyheatmap_per_dataset.output,**get_wildcards(parameters,'dataset')),
        expand(rules.funkyheatmap_per_method.output,**get_wildcards(parameters,'method')),
        # barplot
        expand(rules.metrics_barplot.output,metric=['s', 'max_uss', 'score']),
        expand(rules.metrics_barplot_per_dataset.output,metric=['s', 'max_uss', 'score'],**get_wildcards(parameters,'dataset')),
        expand(rules.metrics_barplot_per_method.output,metric=['s', 'max_uss', 'score'],**get_wildcards(parameters,'method')),
        # swarmplot
        expand(rules.metrics_swarmplot.output,metric='score'),
        expand(rules.metrics_swarmplot_per_dataset.output,metric='score',**get_wildcards(parameters,'dataset')),
        expand(rules.metrics_swarmplot_per_method.output,metric='score',**get_wildcards(parameters,'method')),
        # implementation comparison
        # rules.compare_metrics.output,
