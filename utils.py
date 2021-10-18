import functools
import glob
from pathlib import Path

from aiapy.calibrate import (update_pointing, fix_observer_location,
                             correct_degradation, normalize_exposure)
from aiapy.calibrate.util import get_correction_table
from astropy.time import Time
import astropy.units as u
import numpy as np
from sunpy.net import Fido, attrs as a
from sunpy.map import Map
import sunpy.sun.constants

from config import fits_dir
from products import *

correction_table = get_correction_table()


def remove_duplicates(map_list):
    fnames = [Path(f).name for f in map_list]
    # Because fnames is sorted, only the larger version number will be saved
    fname_versions = {name[:-7]: int(name[-7:-5]) for name in sorted(fnames)}
    names = [fits_dir / (name + str(version).zfill(2) + '.fits') for
             name, version in fname_versions.items()]
    return names


@functools.lru_cache(maxsize=4)
def get_all_files(product):
    files = glob.glob(str(fits_dir / f'solo_L2_eui-{product.name}-image*.fits'))
    return remove_duplicates(files)


def get_eui_date(fname, product):
    return Time.strptime(fname[:-12], f'solo_L2_eui-{product.name}-image_%Y%m%dT%H%M%S')


@functools.lru_cache(maxsize=4)
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
    """
    Get closest EUI map to a given time.

    Parameters
    ----------
    product : Product
    t : astropy.time.Time
    max_delta : astropy.time.TimeDelta
        Max timedelta threshold.

    Returns
    -------
    sunpy.map.Map
    """
    dates = all_eui_dates(product)
    idx = np.argmin(np.abs(dates - t))
    if np.abs(dates[idx] - t) > max_delta:
        return
    fname = get_all_files(product)[idx]
    return Map(fname)


def get_aia_map(product, t, max_delta):
    """
    Get closest AIA map to a given time that matches an EUI product.

    Parameters
    ----------
    product : Product
    t : astropy.time.Time
    max_delta : astropy.time.TimeDelta
        Max timedelta threshold.

    Returns
    -------
    sunpy.map.Map
    """
    t = a.Time(t - max_delta, t + max_delta, t)
    if product is FSI174:
        wlen = a.Wavelength(171 * u.Angstrom)
    result = Fido.search(t & wlen & a.Instrument('AIA'))
    if len(result[0]) == 0:
        return

    fname = Fido.fetch(result)[0]
    return prep_aia(Map(fname))


def prep_aia(m):
    """
    Prep an AIA map.

    This runs (in order):
        - `aiapy.calibrate.update_pointing`
        - `aiapy.calibrate.fix_observer_location`
        - `aiapy.caibrate.correct_degradation`
        - `aiapy.calibrate.normalize_exposure`
    """
    print('Prepping map')
    if m.exposure_time <= 0 * u.s:
        raise RuntimeError('Exposure time <= 0')
    m = update_pointing(m)
    m = fix_observer_location(m)
    m = correct_degradation(m, correction_table=correction_table)
    m = normalize_exposure(m)
    return m
