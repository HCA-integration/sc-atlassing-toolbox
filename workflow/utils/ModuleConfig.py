import warnings
from pprint import pprint
from typing import Union
from pathlib import Path
import pandas as pd
from snakemake.io import expand, Wildcards
from snakemake.rules import Rule

from .WildcardParameters import WildcardParameters
from .InputFiles import InputFiles
from .config import get_from_config, _get_or_default_from_config
from .misc import create_hash, get_use_gpu


class ModuleConfig:

    # module_name = ''
    # config = {}
    # datasets = {}
    # wildcard_names = []
    # parameters = WildcardParameters()
    # input_files = InputFiles()
    # default_target = None
    # out_dir = Path()
    # image_dir = Path()


    def __init__(
        self,
        module_name: str,
        config: dict,
        parameters: [pd.DataFrame, str] = None,
        default_target: [str, Rule] = None,
        wildcard_names: list = None,
        config_params: list = None,
        rename_config_params: dict = None,
        explode_by: [str, list] = None,
        paramspace_kwargs: dict = None,
    ):
        """
        :param module_name: name of module
        :param config: complete config from Snakemake
        :param parameters: pd.DataFrame or path to TSV with module specific parameters
        :param default_target: default output pattern for module
        :param wildcard_names: list of wildcard names for expanding rules
        :param config_params: list of parameters that a module should consider as wildcards, order and length must match wildcard_names, by default will take wildcard_names
        :param explode_by: column(s) to explode wildcard_names extracted from config by
        :param paramspace_kwargs: arguments passed to WildcardParameters
        """
        self.module_name = module_name
        self.config = config
        self.set_defaults()
        self.set_datasets()

        self.out_dir = Path(self.config['output_dir']) / self.module_name
        self.out_dir.mkdir(parents=True, exist_ok=True)
        self.image_dir = Path(self.config['images']) / self.module_name

        self.input_files = InputFiles(
            module_name=self.module_name,
            dataset_config=self.datasets,
            output_directory=self.out_dir,
        )

        if wildcard_names is None:
            wildcard_names = []

        if isinstance(parameters, str):
            parameters = pd.read_table(parameters)

        self.parameters = WildcardParameters(
            module_name=module_name,
            parameters=parameters,
            input_file_wildcards=self.input_files.get_wildcards(),
            dataset_config=self.datasets,
            default_config=self.config.get('defaults'),
            wildcard_names=wildcard_names,
            config_params=config_params,
            rename_config_params=rename_config_params,
            explode_by=explode_by,
            paramspace_kwargs=paramspace_kwargs,
        )

        self.set_default_target(default_target)


    def set_defaults(self, warn: bool = False):
        self.config['DATASETS'] = self.config.get('DATASETS', {})
        self.config['defaults'] = self.config.get('defaults', {})
        self.config['defaults'][self.module_name] = self.config['defaults'].get(self.module_name, {})
        datasets = list(self.config['DATASETS'].keys())
        self.config['defaults']['datasets'] = self.config['defaults'].get('datasets', datasets)
        self.config['output_map'] = self.config.get('output_map', {})


    def set_defaults_per_dataset(self, dataset: str, warn: bool = False):
        # update entries for each dataset
        entry = _get_or_default_from_config(
            config=self.config['DATASETS'],
            defaults=self.config['defaults'],
            key=dataset,
            value=self.module_name,
            return_missing={},
            warn=warn,
            update=True,
        )
        
        # for TSV input make sure integration methods have the proper types TODO: deprecate
        if self.module_name == 'integration' and isinstance(entry, list):
            entry = {k: self.config['defaults'][self.module_name][k] for k in entry}
        
        # set entry in config
        self.config['DATASETS'][dataset][self.module_name] = entry
        self.datasets[dataset][self.module_name] = entry


    def set_datasets(self):
        """
        Set dataset configs
        """
        self.datasets = {
            dataset: entry for dataset, entry in self.config['DATASETS'].items()
            if self.module_name in entry.get('input', {}).keys()
        }
        
        for dataset in self.datasets:
            self.set_defaults_per_dataset(dataset)


    def set_default_target(self, default_target: [str, Rule] = None):
        if default_target is not None:
            self.default_target = default_target
        elif self.module_name in self.config['output_map']:
            self.default_target = self.config['output_map'][self.module_name]
        else:
            default_target = self.out_dir / self.parameters.get_paramspace().wildcard_pattern / f'{self.module_name}.tsv'
            warnings.warn(f'\nNo default target specified for module "{self.module_name}", using "{default_target}"')
            self.default_target = default_target


    def get_defaults(self, module_name: str = None):
        """
        Get defaults for module
        """
        if module_name is not None:
            module_name = self.module_name
        return self.config['defaults'].get(module_name, {}).copy()


    def get_for_dataset(
        self,
        dataset: str,
        query: list,
        default: Union[str,bool,float,int,dict,list, None] = None,
        # module_name: str = None,
        warn: bool = False,
    ) -> Union[str,bool,float,int,dict,list, None]:
        """Get any key from the config via query

        Args:
            dataset (str): dataset key in config['DATASETS']
            query (list): list of keys to walk down the config
            default (Union[str,bool,float,int,dict,list, None], optional): default value if key not found. Defaults to None.

        Returns:
            Union[str,bool,float,int,dict,list, None]: value of query in config
        """
        return get_from_config(
            config=self.config['DATASETS'][dataset],
            query=query,
            default=default,
            warn=warn
        )


    def get_config(self) -> dict:
        """
        Get complete config with parsed input files
        """
        return self.config


    def get_datasets(self) -> dict:
        """
        Get config for datasets that use the module
        """
        return self.datasets


    def get_wildcards(self, **kwargs) -> [dict, pd.DataFrame]:
        """
        Retrieve wildcard instances as dictionary
        
        :param **kwargs: arguments passed to WildcardParameters.get_wildcards
        :return: dictionary of wildcards that can be applied directly for expanding target files
        """
        return self.parameters.get_wildcards(**kwargs)


    def get_wildcard_names(self):
        return self.parameters.wildcard_names


    def get_output_files(
        self,
        pattern: [str, Rule] = None,
        allow_missing: bool=False,
        as_dict: bool=False,
        **kwargs
    ) -> list:
        """
        Get output file based on wildcards
        :param pattern: output pattern, defaults to self.default_target
        :param **kwargs: arguments passed to WildcardParameters.get_wildcards
        """
        if pattern is None:
            pattern = self.default_target
        wildcards = self.get_wildcards(**kwargs)
        targets = expand(pattern, zip, **wildcards, allow_missing=allow_missing)
        if as_dict:
            
            def get_wildcard_string(wildcard_name, wildcard_value):
                if wildcard_name == 'file_id':
                    if '/' in wildcard_value:
                        wildcard_value = create_hash(wildcard_value)
                    else:
                        return f'{self.module_name}:{wildcard_value}'
                return f'{self.module_name}_{wildcard_name}={wildcard_value}'
            
            def shorten_name(name):
                split_values = name.split(f'--{self.module_name}_', 1)
                if len(split_values) > 1 and len(name) > 300:
                    name = f'{split_values[0]}--{self.module_name}={create_hash(split_values[1])}'
                return name
            
            task_names = [
                '--'.join([get_wildcard_string(k, v) for k, v in zip(wildcards.keys(), w) if k != 'dataset'])
                for w in zip(*wildcards.values())
            ]
            task_names = [shorten_name(name) for name in task_names]
            targets = dict(zip(task_names, targets))
        return targets


    def get_input_file_wildcards(self):
        return self.input_files.get_wildcards()


    def get_input_files(self):
        return self.input_files.get_files()


    def get_input_files_per_dataset(self, dataset):
        return self.input_files.get_files_per_dataset(dataset)


    def get_input_file(self, dataset, file_id, **kwargs):
        return self.input_files.get_file(dataset, file_id)


    def get_parameters(self):
        return self.parameters.get_parameters()


    def get_paramspace(self, **kwargs):
        return self.parameters.get_paramspace(**kwargs)


    def get_from_parameters(self, query_dict: [dict, Wildcards], parameter_key: str, **kwargs):
        return self.parameters.get_from_parameters(dict(query_dict), parameter_key, **kwargs)


    def update_inputs(self, dataset: str, input_files: [str, dict]):
        self.config['DATASETS'][dataset]['input'][self.module_name] = input_files
        self.set_datasets()
        self.input_files = InputFiles(
            module_name=self.module_name,
            dataset_config=self.datasets,
            output_directory=self.out_dir,
        )


    def update_parameters(
        self,
        wildcards_df: pd.DataFrame = None,
        parameters_df: pd.DataFrame = None,
        wildcard_names: list = None,
        **paramspace_kwargs
    ):
        """
        Update parameters and all dependent attributes
        
        :param wildcards_df: dataframe with updated wildcards
        :param parameters_df: dataframe with updated parameters
        :param wildcard_names: list of wildcard names to subset the paramspace by
        :param paramspace_kwargs: additional arguments for snakemake.utils.Paramspace
        """
        self.parameters.update(
            wildcards_df=wildcards_df,
            wildcard_names=wildcard_names,
            **paramspace_kwargs,
        )
        self.set_default_target(self.default_target)


    def get_profile(self, wildcards: [dict, Wildcards]):
        return self.get_from_parameters(wildcards, 'resources')


    def get_resource(
        self,
        resource_key: str,
        profile: str = 'cpu',
        attempt: int = 1,
        factor: float = 0.5
    ) -> [str, int, float]:
        """
        Retrieve resource information from config['resources']
        
        :param profile: resource profile, key under config['resources']
        :param resource_key: resource key, key under config['resources'][profile]
        """
        if 'resources' not in self.config or not profile:
            return ''
        # overwrite profile to cpu if turned off in config
        profile = profile if get_use_gpu(self.config) else 'cpu'
        
        if attempt > 2:
            profile = 'cpu'
        
        resources = self.config['resources']
        try:
            res = resources[profile][resource_key]
        except KeyError:
            print(
                f'WARNING: Invalid profile "{profile}" or resource key "{resource_key}". '
                'Please check that your config contains the correct entries under config["resources"]'
            )
            return ''
        if resource_key == 'mem_mb':
            return int(res + (attempt - 1) * factor * res)
        return res
