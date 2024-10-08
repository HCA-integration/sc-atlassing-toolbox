from .misc import ifelse, get_use_gpu


def get_env(
    config: dict,
    env_name: str,
    mode: bool = None,
    env_dir: str = '../envs',
    gpu_env: bool = None,
    no_gpu: bool = False,
):
    if mode is None:
        mode = config.get('env_mode', 'local')
    if no_gpu or gpu_env is None:
        gpu_env = env_name
    env_name = ifelse(get_use_gpu(config), _if=gpu_env, _else=env_name)

    if mode == 'from_yaml':
        return f'{env_dir}/{env_name}.yaml'
    elif mode == 'local':
        return env_name
    else:
        raise ValueError(f'Unknown environment mode: {mode}')

