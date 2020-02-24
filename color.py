"""Script to combine stacked frames into color images.

"""
from fits_utils import *

def main():
    unaligned_images = load_fits(path="sta/", target="m52", band="")
    aligned_images = align(unaligned_images, centroid=hybrid_centroid, filter="combined")
    rgb(smooth(aligned_images[0]["data"],sigma=2.5),
        smooth(aligned_images[1]["data"],sigma=2.5),
        smooth(aligned_images[2]["data"],sigma=2.5)
        )

if __name__ == '__main__':
    main()
