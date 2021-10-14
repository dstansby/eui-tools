import sunpy.net.attrs as a
from sunpy.net import Fido
from sunpy_soar.attrs import Identifier


def get_eui(product, start_time, end_time, dl_path):
    """
    Dowload EUI data.

    Parameters
    ----------
    product: prodcuts.Product
        Data product. See ``products.py`` for the four availalble products.
    start_time:
        Search start time.
    end_time:
        Search end time.
    dl_path:
        Directory to download files to.

    Returns
    -------
    files: list[str]
        List of downloaded file paths.
    """
    id = a.Instrument('EUI')
    time = a.Time(start_time, end_time)
    level = a.Level(2)
    identifier = Identifier(f'EUI-{product.name.upper()}.-IMAGE')

    res = Fido.search(id, time, level, identifier)
    print(res)
    files = Fido.fetch(res, path=dl_path)
    return files
