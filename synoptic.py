from astropy.time import Time, TimeDelta
import astropy.units as u

from products import *
from utils import get_closest_map


def generate_synoptic_images(product):
    assert product in [FSI174, FSI304], 'Must be FSI'

    t = Time('2020-11-19 12:00:00')
    endtime = Time('2020-12-01 00:00:00')
    delta = TimeDelta(12 * u.h)

    while t < endtime:
        m = get_closest_map(product, t, delta / 2)
        if m is None:
            print(f'💔  No EUI data available near {t.iso}')
            continue
        print(m)

        t += delta
