import warnings
import os
from pathlib import Path
import shutil
import anndata as ad
import zarr
import h5py
from scipy.sparse import csr_matrix
from anndata.experimental import read_elem, sparse_dataset

zarr.default_compressor = zarr.Blosc(shuffle=zarr.Blosc.SHUFFLE)


def get_file_reader(file):
    if file.endswith(('.zarr', '.zarr/')):
        func = zarr.open
        file_type = 'zarr'
    elif file.endswith('.h5ad'):
        func = h5py.File
        file_type = 'h5py'
    else:
        raise ValueError(f'Unknown file format: {file}')
    return func, file_type


def check_slot_exists(file, slot):
    func, file_type = get_file_reader(file)
    with func(file, 'r') as f:
        exists = slot in f
    return exists


def to_memory(matrix):
    if isinstance(matrix, (ad.experimental.CSRDataset, ad.experimental.CSCDataset)):
        print('Convert to memory...')
        return matrix.to_memory()
    return matrix


def read_anndata(
    file: str,
    dask: bool = False,
    backed: bool = False,
    **kwargs
) -> ad.AnnData:
    """
    Read anndata file
    :param file: path to anndata file in zarr or h5ad format
    :param kwargs: AnnData parameter to zarr group mapping
    """
    # assert Path(file).exists(), f'File not found: {file}'
    
    if dask:
        read_func = read_dask
    elif backed:
        read_func = read_partial
    else:
        read_func = read_partial
    
    func, file_type = get_file_reader(file)
    store = func(file, 'r')
    
    # set default kwargs
    kwargs = {x: x for x in store} if not kwargs else kwargs
    # set key == value if value is None
    kwargs |= {k: k for k, v in kwargs.items() if v is None}
    
    # return an empty AnnData object if no keys are available
    if len(store.keys()) == 0:
        return ad.AnnData()
    
    # check if keys are available
    for name, slot in kwargs.items():
        if slot not in store:
            warnings.warn(
                f'Cannot find "{slot}" for AnnData parameter `{name}`'
                ' from "{file}", will be skipped'
            )
    adata = read_func(store, backed=backed, **kwargs)
    if not backed and file_type == 'h5py':
        store.close()
    
    return adata


def read_partial(
    group: [h5py.Group, zarr.Group],
    backed: bool = False,
    force_sparse_types: [str, list] = None,
    **kwargs
) -> ad.AnnData:
    """
    Partially read zarr or h5py groups
    :params group: file group
    :params force_sparse_types: encoding types to convert to sparse_dataset via csr_matrix
    :params backed: read sparse matrix as sparse_dataset
    :params **kwargs: dict of slot_name: slot, by default use all available slot for the zarr file
    :return: AnnData object
    """
    if force_sparse_types is None:
        force_sparse_types = []
    elif isinstance(force_sparse_types, str):
        force_sparse_types = [force_sparse_types]
    slots = {}
    if backed:
        print('Read as backed sparse matrix...')
    
    for slot_name, slot in kwargs.items():
        print(f'Read slot "{slot}", store as "{slot_name}"...')
        if slot not in group:
            warnings.warn(f'Slot "{slot}" not found, skip...')
            slots[slot_name] = None
        else:
            elem = group[slot]
            iospec = ad._io.specs.get_spec(elem)
            if iospec.encoding_type in ("csr_matrix", "csc_matrix") and backed:
                slots[slot_name] = sparse_dataset(elem)
            elif iospec.encoding_type in force_sparse_types:
                slots[slot_name] = csr_matrix(read_elem(elem))
                if backed:
                    slots[slot_name] = sparse_dataset(slots[slot_name])
            else:
                slots[slot_name] = read_elem(elem)
    return ad.AnnData(**slots)


def read_dask(
    group: [h5py.Group, zarr.Group],
    backed: bool = False,
    obs_chunk: int = 1000,
    **kwargs
) -> ad.AnnData:
    """
    Modified from https://anndata.readthedocs.io/en/latest/tutorials/notebooks/%7Bread%2Cwrite%7D_dispatched.html
    """
    from anndata.experimental import read_dispatched
    from utils.sparse_dask import sparse_dataset_as_dask
    
    def callback(func, elem_name: str, elem, iospec):
        import re
        import dask.array as da
        import sparse
        
        elem_matches = [
            not (
                bool(re.match(f'/{e}(/.?|$)', elem_name)) 
                or f'/{e}'.startswith(elem_name) 
            )
            for e in kwargs.values()
        ]
        if elem_name != '/' and all(elem_matches):
            print('skip reading', elem_name)
            return None
        else:
            print('read', elem_name)
        
        if elem_name != '/' and all(elem_matches):
            print('skip reading', elem_name)
            return None
        elif iospec.encoding_type in (
            "dataframe",
            "awkward-array",
        ):
            # Preventing recursing inside of these types
            return read_elem(elem)
        elif iospec.encoding_type in ("csr_matrix", "csc_matrix"):
            # return da.from_array(read_elem(elem))
            matrix = sparse_dataset_as_dask(sparse_dataset(elem), obs_chunk)
            return matrix.map_blocks(sparse.COO)
        elif iospec.encoding_type == "array":
            return da.from_zarr(elem)
        return func(elem)

    return read_dispatched(group, callback=callback)


# deprecated
def read_anndata_or_mudata(file):
    if file.endswith('.h5mu'):
        import mudata as mu
        print('Read as mudata...')
        return mu.read(file)
    elif file.endswith('.h5mu.zarr'):
        import mudata as mu
        print('Read as mudata from zarr...')
        return mu.read_zarr(file)
    else:
        print('Read as anndata...')
        return read_anndata(file)


def write_zarr(adata, file):
    def sparse_coo_to_csr(matrix):
        from dask.array import Array as DaskArray
        import sparse
        
        if isinstance(matrix, DaskArray) and isinstance(matrix._meta, sparse.COO):
            matrix = matrix.map_blocks(csr_matrix, dtype='float32')
        return matrix
    
    adata.X = sparse_coo_to_csr(adata.X)
    for layer in adata.layers:
        adata.layers[layer] = sparse_coo_to_csr(adata.layers[layer])
    adata.write_zarr(file) # doesn't seem to work with dask array


def link_zarr(
    in_dir: [str, Path],
    out_dir: [str, Path],
    file_names: list = None,
    overwrite: bool = False,
    relative_path: bool = True,
    slot_map: dict = None,
):
    """
    Link to existing zarr file
    :param in_dir: path to existing zarr file
    :param out_dir: path to output zarr file
    :param file_names: list of files to link, if None, link all files
    :param overwrite: overwrite existing output files
    :param relative_path: use relative path for link
    :param kwargs: custom mapping of output slot to input slot,
        will update default mapping of same input and output naming
    """
    def resolved_nested_links(d):
        # determine equivalent classes of slots (top hierarchy)
        eq_classes = {}
        for out_slot in d:
            if '/' not in out_slot:
                continue
            eq = out_slot.rsplit('/', 1)[0]
            if eq not in d:
                continue
            eq_classes.setdefault(eq, []).append(out_slot)

        for out_slot in eq_classes.keys():
            in_slot = d[out_slot]
            for f in (in_dir / in_slot).iterdir():
                new_out_slot = f'{out_slot}/{f.name}'
                if new_out_slot in d or f.name == '.snakemake_timestamp':
                    continue
                d[new_out_slot] = f'{in_slot}/{f.name}'
            del d[out_slot]
        return d
    
    def link_file(in_file, out_file, relative_path=True):
        in_file = in_file.resolve()
        out_dir = out_file.parent.resolve()
        out_dir.mkdir(parents=True, exist_ok=True)
        
        if relative_path:
            in_file = Path(os.path.relpath(in_file, out_dir))
        
        if overwrite and out_file.exists():
            if out_file.is_dir() and not out_file.is_symlink():
                shutil.rmtree(out_file)
            else:
                out_file.unlink()
        
        out_file.symlink_to(in_file)

    in_dir = Path(in_dir)
    out_dir = Path(out_dir)
    
    if not in_dir.exists():
        raise ValueError(f'Input directory {in_dir} does not exist')
    
    if file_names is None:
        file_names = [f.name for f in in_dir.iterdir()]
    file_names = [
        file for file in file_names
        if file not in ('.snakemake_timestamp')
    ]
    
    if slot_map is None:
        slot_map = {}
    slot_map = {file: file for file in file_names} | slot_map
    slot_map |= {k: k for k, v in slot_map.items() if v is None}
    # deal with nested mapping
    slot_map = resolved_nested_links(slot_map)

    for out_slot, in_slot in slot_map.items():
        in_file_name = str(in_dir).split('.zarr/')[-1] + '/' + in_slot
        out_file_name = str(out_dir).split('.zarr/')[-1] + '/' + out_slot
        print(f'Link {out_file_name} -> {in_file_name}')
        link_file(
            in_file=in_dir / in_slot,
            out_file=out_dir / out_slot,
            relative_path=relative_path
        )


def link_zarr_partial(in_dir, out_dir, files_to_keep=None, overwrite=True, relative_path=True):
    """
    Link zarr files excluding defined slots
    """
    if not in_dir.endswith('.zarr'):
        return
    if files_to_keep is None:
        files_to_keep = []
    in_dirs = [f.name for f in Path(in_dir).iterdir()]
    files_to_link = [f for f in in_dirs if f not in files_to_keep]
    link_zarr(
        in_dir=in_dir,
        out_dir=out_dir,
        file_names=files_to_link,
        overwrite=overwrite,
        relative_path=relative_path,
    )


def write_zarr_linked(
    adata: ad.AnnData,
    in_dir: [str, Path],
    out_dir: [str, Path],
    relative_path: bool = True,
    files_to_keep: list = None,
    slot_map: dict = None,
):
    """
    Write adata to linked zarr file
    :param adata: AnnData object
    :param in_dir: path to existing zarr file
    :param out_dir: path to output zarr file
    :param files_to_keep: list of files to keep and not overwrite
    :param relative_path: use relative path for link
    :param slot_map: custom mapping of output slot to input slot, for slots that are not in files_to_keep
    """
    in_dir = Path(in_dir)
    
    if not in_dir.name.endswith('.zarr'):
        adata.write_zarr(out_dir)
        return
    
    if files_to_keep is None:
        files_to_keep = []
    in_dirs = [f.name for f in in_dir.iterdir()]
    files_to_link = [f for f in in_dirs if f not in files_to_keep]
    
    if slot_map is None:
        slot_map = {}
    extra_slots_to_link = list(slot_map.keys())
    
    # keep only slots that are not explicitly in files_to_keep
    slot_map = {
        in_slot: out_slot 
        for in_slot, out_slot in slot_map.items()
        if in_slot not in files_to_keep
    }
    
    # remove slots that will be overwritten anyway
    for slot in files_to_link+extra_slots_to_link:
        if slot in adata.__dict__:
            print(f'remove {slot}...')
            delattr(adata, slot)
    
    # write zarr file
    adata.write_zarr(out_dir)
    
    # link files
    link_zarr(
        in_dir=in_dir,
        out_dir=out_dir,
        file_names=files_to_link,
        overwrite=True,
        relative_path=relative_path,
        slot_map=slot_map,
    )
