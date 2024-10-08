rule run_method:
    message:
       """
       Sample representation: Run method={wildcards.method} on dataset={wildcards.dataset} file_id={wildcards.file_id}
       input: {input}
       output: {output}
       wildcards: {wildcards}
       resources: gpu={resources.gpu} mem_mb={resources.mem_mb} partition={resources.partition} qos={resources.qos}
       """
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / f'{paramspace.wildcard_pattern}.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        cell_type_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cell_type_key'),
        use_rep=lambda wildcards: mcfg.get_from_parameters(wildcards, 'use_rep'),
        var_mask=lambda wildcards: mcfg.get_from_parameters(wildcards, 'var_mask'),
        hyperparams=lambda wildcards: mcfg.get_from_parameters(wildcards, 'hyperparams_dict'),
        env=lambda wildcards: mcfg.get_from_parameters(wildcards, 'env'),
        script_suffix=lambda wildcards: mcfg.get_from_parameters(wildcards, 'script_suffix'),
        cran_url=config.get('cran_url', 'https://cloud.r-project.org'),
    conda:
        lambda wildcards, params: get_env(config, params.env)
    threads:
        lambda wildcards: max(1, mcfg.get_from_parameters(wildcards, 'threads', default=1)),
    resources:
        partition=lambda w, attempt: mcfg.get_resource(resource_key='partition', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=2),
        qos=lambda w, attempt: mcfg.get_resource(resource_key='qos', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=2),
        mem_mb=lambda w, attempt: mcfg.get_resource(resource_key='mem_mb', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=2, factor=1),
        gpu=lambda w, attempt: mcfg.get_resource(resource_key='gpu', profile=mcfg.get_profile(w), attempt=attempt, attempt_to_cpu=2),
        time="2-00:00:00",
    script:
        '../scripts/methods/{wildcards.method}.{params.script_suffix}'


rule pseudobulk:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'pseudobulk.zarr'),
    params:
        bulk_by=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key')
    conda:
        get_env(config, 'scanpy')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/methods/pseudobulk.py'

# pca_1_2=mcfg.image_dir / paramspace.wildcard_pattern / 'pc_1_2.png',
# pca_2_3=mcfg.image_dir / paramspace.wildcard_pattern / 'pca_2_3.png',
# pca_scree=mcfg.image_dir / paramspace.wildcard_pattern / 'pca_screeplot.png',


# use rule pseudobulk as pseudobulk_organ with:
#     input:
#         zarr=rules.load_data_merge_organ.output.zarr
#     params:
#         dataset=lambda wildcards: wildcards.organ,
#         bulk_by='sample',
#         color=[
#             'dataset',
#             'study',
#             'reference',
#             'sex',
#             'disease',
#             'assay',
#             'modalities',
#             'pipeline_version',
#             'institution',
#             'self_reported_ethnicity',
#             'development_stage'
#         ],
#     output:
#         pca_1_2=images_dir / 'organ' / '{organ}_1_2.png',
#         pca_2_3=images_dir / 'organ' / '{organ}_2_3.png',
#         pca_scree=images_dir / 'organ' / '{organ}_scree.png',
#     resources:
#         mem_mb=get_resource(config,profile='cpu_merged',resource_key='mem_mb')

rule composition:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'composition.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        cell_type_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cell_type_key'),
        use_rep=lambda wildcards: mcfg.get_from_parameters(wildcards, 'use_rep')
    conda:
        get_env(config, 'sample_representation')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/methods/composition.py'


rule cell_type_pseudobulk:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'cell_type_pseudobulk.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        cell_type_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cell_type_key'),
        use_rep=lambda wildcards: mcfg.get_from_parameters(wildcards, 'use_rep')
    conda:
        get_env(config, 'sample_representation')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/methods/cell_type_pseudobulk.py'


rule pilot:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'pilot.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        cell_type_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cell_type_key'),
        use_rep=lambda wildcards: mcfg.get_from_parameters(wildcards, 'use_rep')    
    conda:
        get_env(config, 'sample_representation')
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/methods/pilot.py'

rule scpoli:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'scpoli.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        cell_type_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cell_type_key'),
        use_rep=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts'),
        n_epochs=lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_epochs'),
    conda:
        get_env(config, 'sample_representation')
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt,attempt_to_cpu=2),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt,attempt_to_cpu=2),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt,attempt_to_cpu=2),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt,attempt_to_cpu=2),
    script:
        '../scripts/methods/scpoli.py'

rule mrvi:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'mrvi.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        cell_type_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cell_type_key'),
        use_rep=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts'),
        n_epochs=lambda wildcards: mcfg.get_from_parameters(wildcards, 'n_epochs'),
    conda:
        get_env(config, 'scvi-tools')
    resources:
        partition=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='partition',attempt=attempt,attempt_to_cpu=2),
        qos=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='qos',attempt=attempt,attempt_to_cpu=2),
        gpu=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='gpu',attempt=attempt,attempt_to_cpu=2),
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='gpu',resource_key='mem_mb',attempt=attempt,attempt_to_cpu=2),
    script:
        '../scripts/methods/mrvi.py'


rule gloscope:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'gloscope.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        use_rep=lambda wildcards: mcfg.get_from_parameters(wildcards, 'use_rep'),
        k=25,
        seed=0,
    conda:
        get_env(config, 'gloscope')
    threads: 10
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/methods/gloscope.R'


rule scitd:
    input:
        zarr=lambda wildcards: mcfg.get_input_file(**wildcards),
        prepare=rules.prepare.output.zarr,
    output:
        zarr=directory(mcfg.out_dir / paramspace.wildcard_pattern / 'scitd.zarr'),
    params:
        sample_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'sample_key'),
        cell_type_key=lambda wildcards: mcfg.get_from_parameters(wildcards, 'cell_type_key'),
        use_rep=lambda wildcards: mcfg.get_from_parameters(wildcards, 'raw_counts'),
        seed=0,
        cran_url=config.get('cran_url', 'https://cloud.r-project.org'),
    conda:
        get_env(config, 'scitd')
    threads: 10
    resources:
        mem_mb=lambda w, attempt: mcfg.get_resource(profile='cpu',resource_key='mem_mb', attempt=attempt)
    script:
        '../scripts/methods/scitd.R'
