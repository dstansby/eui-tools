from pathlib import Path
import subprocess
from products import *

for product in all_products:
    print(product.name)
    eui_dir = Path('/Volumes/Work/Data/solo')
    png_glob = eui_dir / 'eui_png' / product.name / '*.png'
    output = eui_dir / f'eui_mp4/eui_{product.name}.mp4'
    args = f'ffmpeg -framerate 24 -pattern_type glob -i {png_glob} -pix_fmt yuv420p -vf scale=1000:800 {output}'.split(' ')
    subprocess.run(args)
