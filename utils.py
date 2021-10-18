import functools
import glob
from pathlib import Path

from astropy.time import Time
import numpy as np

from config import fits_dir


def remove_duplicates(map_list):
    fnames = [Path(f).name for f in map_list]
    # Because fnames is sorted, only the larger version number will be saved
    fname_versions = {name[:-7]: int(name[-7:-5]) for name in sorted(fnames)}
    names = [fits_dir / (name + str(version).zfill(2) + '.fits') for
             name, version in fname_versions.items()]
    return names


@functools.lru_cache()
def get_all_files(product):
    files = glob.glob(str(fits_dir / f'solo_L2_eui-{product.name}-image*.fits'))
    return remove_duplicates(files)


def get_eui_date(fname, product):
    return Time.strptime(fname[:-12], f'solo_L2_eui-{product.name}-image_%Y%m%dT%H%M%S')


@functools.lru_cache()
def all_eui_dates(product):
    """
    Get the dates of all available files for a specific product.

    Parameters
    ----------
    product : Products

    Returns
    -------
    astropy.time.Time
    """
    fnames = [f.name for f in get_all_files(product)]
    dates = [get_eui_date(f, product) for f in fnames]
    dates = Time(dates)
    return dates


def get_closest_map(product, t, max_delta):
    dates = all_eui_dates(product)
    idx = np.argmin(np.abs(dates - t))
    if np.abs(dates[idx] - t) > max_delta:
        return
    return get_all_files(product)[idx]
