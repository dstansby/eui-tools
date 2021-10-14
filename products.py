from dataclasses import dataclass


@dataclass
class Product:
    name: str
    cmap: str


FSI174 = Product('fsi174', 'solar orbiterfsi174')
FSI304 = Product('fsi304', 'solar orbiterfsi304')
HRIEUV174 = Product('hrieuv174', 'solar orbiterhri_lya1216')
HRILYA1216 = Product('hrilya1216', 'solar orbiterhri_euv174')

all_products = [FSI174, FSI304, HRIEUV174, HRILYA1216]


def check_product(product):
    """
    Parameters
    ----------
    product: Product
    """
    if product not in all_products:
        raise ValueError(f'{product} must be one of {all_products}')
