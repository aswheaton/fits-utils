"""Script to combine stacked frames into color images.

"""
from fits_utils import *

def main():
    unaligned_images = load_fits(path="sta/", target="m52", band="")
    # unaligned_images += load_fits(path="sta/", target="m52", band="g")
    # unaligned_images += load_fits(path="sta/", target="m52", band="u")
    aligned_images = align(unaligned_images, centroid=hybrid_centroid, mask=True)
    rgb(aligned_images[0],aligned_images[1],aligned_images[2])

if __name__ == '__main__':
    main()
