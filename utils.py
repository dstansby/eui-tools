import glob
from pathlib import Path

from astropy.time import Time

from config import fits_dir


def remove_duplicates(map_list):
    fnames = [Path(f).name for f in map_list]
    # Because fnames is sorted, only the larger version number will be saved
    fname_versions = {name[:-7]: int(name[-7:-5]) for name in sorted(fnames)}
    names = [fits_dir / (name + str(version).zfill(2) + '.fits') for
             name, version in fname_versions.items()]
    return names


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
    files = glob.glob(str(fits_dir / f'solo_L2_eui-{product.name}-image*.fits'))
    files = remove_duplicates(files)
    fnames = [f.name for f in files]
    dates = [Time.strptime(f[:-12], f'solo_L2_eui-{product.name}-image_%Y%m%dT%H%M%S') for f in fnames]
    dates = Time(dates)
    return dates
