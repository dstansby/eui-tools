import glob
from pathlib import Path
import os

from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.colors as mcolor
import numpy as np
import sunpy.map

from products import *

fits_dir = Path('/Volumes/Work/Data/solo/eui_fits')


def remove_duplicates(map_list):
    fnames = [Path(f).name for f in map_list]
    # Because fnames is sorted, only the larger version number will be saved
    fname_versions = {name[:-7]: int(name[-7:-5]) for name in sorted(fnames)}
    names = [fits_dir / (name + str(version).zfill(2) + '.fits') for name, version in fname_versions.items()]
    return names


def get_map_paths(product):
    """
    Get all maps of a given data product, remove version duplicates.

    Parameters
    ----------
    product: product

    Returns
    -------
    list[str]
        List of all detected filepaths.
    """
    check_product(product)
    euimaps = sorted(glob.glob(f'{fits_dir}/solo_L2_eui-{product.name}-image_*.fits'))
    euimaps = remove_duplicates(euimaps)
    return euimaps


def all_fits_to_png(product):
    euimaps = get_map_paths(product)

    cmap = plt.get_cmap(product.cmap)
    png_dir = f'/Volumes/Work/Data/solo/eui_png/{product.name}'
    broken_files = []

    for f in euimaps:
        date = Time.strptime(f.name[-27:-12], '%Y%m%dT%H%M%S')
        datestr = date.strftime('%Y_%m_%d-%H%M%S')

        png_fname = f'{png_dir}/EUI_{datestr}.png'
        if Path(png_fname).exists():
            print(f'✅ {date.isot}')
            continue

        try:
            m = sunpy.map.Map(f)
        except Exception:
            print(f'❌ {f}')
            broken_files.append(f)
            continue
        m.meta['crota2'] = m.meta.pop('crota')

        print(m.date.isot)
        m.data[m.data < 1] = 1
        fov = m.rsun_obs * 1.3
        bottom_left = m.world_to_pixel(
            SkyCoord(-fov, -fov, frame=m.coordinate_frame))
        top_right = m.world_to_pixel(
            SkyCoord(fov, fov, frame=m.coordinate_frame))

        fig = plt.figure()

        vmin, vmax = product.vmin, product.vmax
        if vmin is None:
            vmin = max(np.nanpercentile(m.data, 1), 10)
        if vmax is None:
            vmax = max(np.nanpercentile(m.data, 100), 100)
            vmax = 10**(np.ceil(np.log10(vmax)))

        norm = mcolor.LogNorm(vmin=vmin, vmax=vmax)
        m.plot(cmap=cmap, norm=norm)
        m.draw_limb()

        plt.colorbar(label=m.unit.to_string(), extend='both')
        plt.gca().set_facecolor('black')
        plt.xlim(bottom_left.x.value, top_right.x.value)
        plt.ylim(bottom_left.y.value, top_right.y.value)
        fig.savefig(png_fname)
        plt.close('all')

    if len(broken_files):
        print(f'Found {len(broken_files)} broken files')
        if input('Remove files? [y/N]: ') == 'y':
            for f in broken_files:
                os.remove(f)
