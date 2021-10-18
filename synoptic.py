import astropy.units as u
import matplotlib.colors as mcolor
import matplotlib.pyplot as plt
import numpy as np
import sunpy.sun.constants

from astropy.coordinates import SkyCoord
from astropy.time import Time, TimeDelta
from astropy.wcs import WCS
from reproject import reproject_interp
from sunpy.map import Map, make_fitswcs_header

from products import *
from utils import get_closest_map, get_aia_map

SHAPE_OUT = [360, 720]


def generate_synoptic_images(product):
    assert product in [FSI174, FSI304], 'Must be FSI'

    t = Time('2020-11-19 12:00:00')
    endtime = Time('2020-12-01 00:00:00')
    delta = TimeDelta(12 * u.h)

    while t < endtime:
        eui_map = get_closest_map(product, t, delta / 2)
        if eui_map is None:
            print(f'💔  No EUI data available near {t.iso}')
            t += delta
            continue

        aia_map = get_aia_map(product, t, delta)
        if aia_map is None:
            print(f'💔💔💔  No AIA data available near {t.iso}')
            t += delta
            continue

        eui_reproj = reproject_carrington(eui_map, t)
        eui_reproj.plot_settings['norm'] = mcolor.LogNorm(vmin=1e1, vmax=1e4)
        eui_reproj.plot_settings['cmap'] = 'sdoaia171'

        aia_reproj = reproject_carrington(aia_map, t)

        total = np.nanmean(np.stack((eui_reproj.data,
                                     aia_reproj.data * 2),
                                    axis=-1),
                           axis=-1)
        total_map = Map(total, aia_reproj.wcs)
        total_map.plot_settings['norm'] = mcolor.LogNorm(vmin=1e1, vmax=1e4)
        total_map.plot_settings['cmap'] = 'sdoaia171'

        fig = plt.figure()
        ax = fig.add_subplot(111, projection=total_map)
        total_map.plot(axes=ax)

        eui_map['date-obs'] = t.isot
        aia_map['date-obs'] = t.isot
        limb_eui, _ = eui_map.draw_limb(axes=ax, color='tab:blue')
        limb_aia, _ = aia_map.draw_limb(axes=ax, color='tab:green')
        ax.legend([limb_eui, limb_aia], ['EUI limb', 'AIA limb'])
        fig.savefig(t.strftime('synoptic_%Y%m%d_%H%M%S.png'), dpi=300)

        t += delta


def synop_header(shape_out, dtime):
    """
    Generate a synoptic map header in a Carrington frame of reference.

    Parameters
    ----------
    shape_out : [int, int]
    dtime : astropy.time.Time

    Returns
    -------
    dict
    """
    frame_out = SkyCoord(180, 0, unit=u.deg,
                         frame="heliographic_carrington",
                         obstime=dtime,
                         observer='earth')
    header = make_fitswcs_header(shape_out, frame_out,
                                 scale=[180 / shape_out[0],
                                        360 / shape_out[1]] * u.deg / u.pix,
                                 projection_code="CAR")
    return header


def reproject_carrington(m, date):
    """
    Reproject a map into Heliographic Carrington coordinates.

    Parameters
    ----------
    m : sunpy.map.Map

    Returns
    -------
    m : sunpy.map.Map
    """
    m.meta['rsun_ref'] = sunpy.sun.constants.radius.to_value(u.m)
    m.meta['date-obs'] = date.isot
    header = synop_header(SHAPE_OUT, m.date)
    wcs = WCS(header)
    # with np.errstate(invalid='ignore'):
    array, footprint = reproject_interp(m, wcs, shape_out=SHAPE_OUT)
    return Map((array, header))
