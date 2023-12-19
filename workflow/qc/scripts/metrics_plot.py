from pathlib import Path
import pandas as pd
from pandas.api.types import is_numeric_dtype
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.axes import Axes
import seaborn as sns
import scanpy as sc
import anndata
import logging
logging.basicConfig(level=logging.INFO)

from utils.io import read_anndata


def plot_qc_joint(
    df,
    x,
    y,
    log_x=1,
    log_y=1,
    hue=None,
    main_plot_function=None,
    marginal_hue=None,
    marginal_legend=False,
    x_threshold=None,
    y_threshold=None,
    title='',
    return_df=False,
    **kwargs,
):
    """
    Plot scatter plot with marginal histograms from df columns.

    :param df: observation dataframe
    :param x: df column for x axis
    :param y: df column for y axis
    :param log: log base for transforming values. Default 1, no transformation
    :param hue: df column with annotations for color coding scatter plot points
    :param marginal_hue: df column with annotations for color coding marginal plot distributions
    :param x_threshold: tuple of upper and lower filter thresholds for x axis
    :param y_threshold: tuple of upper and lower filter thresholds for y axis
    :param title: Title text for plot
    :return:
        seaborn plot (and df dataframe with updated values, if `return_df=True`)
    """
    if main_plot_function is None:
        main_plot_function = sns.scatterplot
    if not x_threshold:
        x_threshold=(0, np.inf)
    if not y_threshold:
        y_threshold=(0, np.inf)

    def log1p_base(_x, base):
        return np.log1p(_x) / np.log(base)

    if log_x > 1:
        x_log = f'log{log_x} {x}'
        df[x_log] = log1p_base(df[x], log_x)
        x_threshold = log1p_base(x_threshold, log_x)
        x = x_log
    
    if log_y > 1:
        y_log = f'log{log_y} {y}'
        df[y_log] = log1p_base(df[y], log_y)
        y_threshold = log1p_base(y_threshold, log_y)
        y = y_log

    g = sns.JointGrid(
        data=df,
        x=x,
        y=y,
        xlim=(0, df[x].max()),
        ylim=(0, df[y].max()),
    )
    # main plot
    g.plot_joint(
        main_plot_function,
        data=df,
        hue=hue,
        **kwargs,
    )
    
    # marginal hist plot
    if marginal_hue in df.columns:
        marginal_hue = None if df[marginal_hue].nunique() > 100 else marginal_hue
    use_marg_hue = marginal_hue is not None
    g.plot_marginals(
        sns.histplot,
        data=df,
        hue=marginal_hue,
        legend=marginal_legend,
        element='step' if use_marg_hue else 'bars',
        fill=False,
        bins=100
    )

    g.fig.suptitle(title, fontsize=12)

    # x threshold
    for t, t_def in zip(x_threshold, (0, np.inf)):
        if t != t_def:
            g.ax_joint.axvline(x=t, color='red')
            g.ax_marg_x.axvline(x=t, color='red')

    # y threshold
    for t, t_def in zip(y_threshold, (0, np.inf)):
        if t != t_def:
            g.ax_joint.axhline(y=t, color='red')
            g.ax_marg_y.axhline(y=t, color='red')

    if return_df:
        return g, df
    return g


input_zarr = snakemake.input.zarr
output_joint = snakemake.output.joint
output_density = snakemake.output.density
# output_violin = snakemake.output.violin
# output_avg = snakemake.output.average_jitter

output_joint = Path(output_joint)
output_joint.mkdir(parents=True, exist_ok=True)


file_id = snakemake.wildcards.file_id
hues = snakemake.params.hue
sample = snakemake.params.sample
dataset = snakemake.params.dataset

logging.info(f'Read {input_zarr}...')
adata = read_anndata(input_zarr, obs='obs', uns='uns')

# if no cells filtered out, save empty plots
if adata.obs.shape[0] == 0:
    logging.info('Save empty plot...')
    plt.savefig(output_density)
    exit()

# set default thresholds
threshold_keys = ['n_counts', 'n_genes', 'percent_mito']
thresholds = {f'{key}_min': 0 for key in threshold_keys} | {f'{key}_max': np.inf for key in threshold_keys}

# get autoqc thresholds
autoqc_thresholds = adata.uns['scautoqc_ranges']
thresholds |= {
    f'{key}_min': autoqc_thresholds.loc[key, 'low']
    for key in autoqc_thresholds.index
}
thresholds |={
    f'{key}_max': autoqc_thresholds.loc[key, 'high']
    for key in autoqc_thresholds.index
}

# update to user thresholds
user_thresholds = snakemake.params.get('thresholds')
if user_thresholds is None:
    user_thresholds = {}
elif isinstance(user_thresholds, str):
    import ast
    user_thresholds = ast.literal_eval(user_thresholds)
elif isinstance(user_thresholds, dict):
    pass
else:
    ValueError('thresholds must be a dict or string')
thresholds |= user_thresholds.get(file_id, {})

# transform to shape expected by plot_qc_joint
thresholds = {
    key: (thresholds[f'{key}_min'], thresholds[f'{key}_max'])
    for key in threshold_keys
}
logging.info(thresholds)

sns.set_theme(style='white')
sc.set_figure_params(frameon=False, fontsize=10, dpi_save=300)
plt.rcParams['figure.figsize'] = 12, 12

split_datasets = dataset.split('--')
if len(split_datasets) > 1:
    dataset = ' '.join([split_datasets[0], split_datasets[-1]])

if isinstance(hues, str):
    hues = [hues]
hues = [hue for hue in hues if hue in adata.obs.columns]
hues = [hue for hue in hues if adata.obs[hue].nunique() > 1]
if len(hues) == 0:
    hues = [None]

logging.info('Joint QC plots...')

# TODO
# density plot
# log scale n_genes
# color by donor
# color by disease

# n_counts vs n_features
x = 'n_counts'
y = 'n_genes'

for hue in hues+['percent_mito']:
    joint_title = f'Joint QC for\n{dataset}\nmargin hue: {hue}'
    if is_numeric_dtype(adata.obs[hue]):
        palette = 'plasma'
        legend = 'brief'
    else:
        palette = None # if adata.obs[hue].nunique() > 100 else 'plasma'
        legend = adata.obs[hue].nunique() <= 20

    plot_qc_joint(
        adata.obs,
        x=x,
        y=y,
        hue=hue,
        marginal_hue=hue,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        title=joint_title,
        s=5,
        alpha=.6,
        linewidth=0,
        palette=palette,
        legend=legend,
    )
    plt.tight_layout()
    plt.savefig(output_joint / f'counts_vs_genes_hue={hue}.png')

    _, adata.obs = plot_qc_joint(
        adata.obs,
        x=x,
        y=y,
        log_x=10,
        log_y=10,
        hue=hue,
        marginal_hue=hue,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        title=joint_title,
        return_df=True,
        s=5,
        alpha=.6,
        linewidth=0,
        palette=palette,
        legend=legend,
    )
    plt.tight_layout()
    plt.savefig(output_joint / f'log_counts_vs_log_genes_hue={hue}.png')


# n_feature x mito frac
x = 'n_genes'
y = 'percent_mito'

for hue in hues:
    joint_title = f'Joint QC for\n{dataset}\nmargin hue: {hue}'
    if is_numeric_dtype(adata.obs[hue]):
        palette = 'plasma'
        legend = 'brief'
    else:
        palette = None # if adata.obs[hue].nunique() > 100 else 'plasma'
        legend = adata.obs[hue].nunique() <= 20
    plot_qc_joint(
        adata.obs,
        x=x,
        y=y,
        # log_x=10,
        hue=hue,
        marginal_hue=hue,
        x_threshold=thresholds[x],
        y_threshold=thresholds[y],
        title=joint_title,
        s=5,
        alpha=.6,
        linewidth=0,
        palette=palette,
        legend=legend,
    )
    plt.tight_layout()
    plt.savefig(output_joint / f'genes_vs_mito_frac_hue={hue}.png')

plot_qc_joint(
    adata.obs.sample(n=int(min(1e5, adata.n_obs)), random_state=42),
    x=x,
    y=y,
    # main_plot_function=Axes.hexbin,
    main_plot_function=sns.kdeplot,
    # log_x=10,
    x_threshold=thresholds[x],
    y_threshold=thresholds[y],
    marginal_hue=sample,
    title=joint_title,
    fill=True,
    cmap='plasma',
    alpha=.8,
)
plt.tight_layout()
plt.savefig(output_density)


# logging.info('Violin plots...')
# fig, axes = plt.subplots(nrows=2, ncols=1)
# hue = hues[0]
# adata.obs[hue] = adata.obs[hue].astype('category')

# sc.pl.violin(
#     adata,
#     keys='log10 n_counts',
#     groupby=hue,
#     rotation=90,
#     show=False,
#     fill=False,
#     ax=axes[0],
# )

# sc.pl.violin(
#     adata,
#     keys='percent_mito',
#     groupby=hue,
#     rotation=90,
#     show=False,
#     fill=False,
#     ax=axes[1],
# )
# fig.suptitle(dataset)

# # g = sns.violinplot(x="donor", y="n_counts", data=adata.obs)
# # g.set_xticklabels(g.get_xticklabels(), rotation=90)
# # g.get_figure().savefig(output_violin)

# plt.tight_layout()
# plt.savefig(output_violin)

# logging.info('Average scatter plots...')
# fig, axes = plt.subplots(nrows=3, ncols=1)

# df = (
#     adata.obs[['log10 n_genes', 'percent_mito', 'log10 n_counts', hue, sample]]
#     .groupby([hue, sample])
#     .median()
#     .reset_index(hue)
#     .dropna()
# )

# sns.stripplot(
#     data=df,
#     x='log10 n_counts',
#     hue=hue,
#     legend=False,
#     ax=axes[0],
# )
# axes[0].set_title('Mean counts per cell, median over sample')
# axes[0].set_xlim(0, None)

# sns.stripplot(
#     data=df,
#     x='log10 n_genes',
#     hue=hue,
#     legend=False,
#     ax=axes[1],
# )
# axes[1].set_title('Mean genes per cell, median over sample')
# axes[1].set_xlim(0, None)

# sns.stripplot(
#     data=df,
#     x='percent_mito',
#     hue=hue,
#     legend=False,
#     ax=axes[2],
# )
# axes[2].set_title('Percent of mitochondrial genes per cell, median over sample')
# axes[2].set_xlim(0, 100)

# fig.suptitle(dataset)

# plt.tight_layout()
# plt.savefig(output_avg)
