from dataclasses import dataclass


@dataclass
class Product:
    name: str
    cmap: str


FSI174 = Product('FSI174', 'solar orbiterfsi174')
