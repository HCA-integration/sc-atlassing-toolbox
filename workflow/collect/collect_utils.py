import pandas as pd

ALL_SLOTS = [
    'X',
    'layers',
    'raw',
    'obs',
    'obsm',
    'obsp',
    'var',
    'varm',
    'varp',
    'uns',
]

def get_same_columns(adatas):
    _ad1 = next(iter(adatas.values()))
    # collect obs columns that exist in all files
    same_obs_columns = set(
        _ad1.obs.columns.tolist()
    ).intersection(
        *[_ad.obs.columns.tolist() for _ad in adatas.values()]
    )
    # determine which columns do not have the same values across all files
    obs_to_remove = []
    for col in same_obs_columns:
        same_across = all(
            _ad1.obs[col].equals(_ad.obs[col])
            for _ad in adatas.values()
        )
        if not same_across:
            obs_to_remove.append(col)
    # remove columns that are not the same across all files
    return [
        col for col in same_obs_columns
        if col not in obs_to_remove
    ]


def merge_df(
    df_current,
    file_id,
    df_previous,
    same_columns,
    sep='_',
):
    if df_previous is None:
        return df_current.rename(
            columns={
                col: f'{col}{sep}{file_id}'
                for col in df_current.columns.tolist()
                if col not in same_columns
            }
        )
    
    assert df_current.index.equals(df_previous.index), \
        f'Index must be the same\n {df_current.index}\n{df_previous.index}'
    unique_columns = [col for col in df_current.columns if col not in same_columns]
    df_current = df_current[unique_columns].rename(
        columns={col: f'{col}{sep}{file_id}' for col in unique_columns}
    )
    df = pd.concat([df_previous, df_current], axis=1)
    return df
