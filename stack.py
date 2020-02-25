"""
Script to align and stack images.
"""
from fits_utils import *

def main():

    # for target in ["bd25", "bd62"]:
    #     for band in ["r", "g"]:
    #         unaligned_images = load_fits(path="sci/", target=target, band=band)
    #         aligned_images = align(unaligned_images, centroid=hybrid_centroid, filter="combined")
    #         stacked_image = stack(aligned_images, correct_exposure=True)
    #         write_out_fits_2(stacked_image, "sta/{}_{}_stacked.fits".format(target, band))
    # for target in ["bd25", "bd62"]:
    #     for band in ["u"]:
    #         unaligned_images = load_fits(path="sci/", target=target, band=band)
    #         aligned_images = align(unaligned_images, centroid=hybrid_centroid, filter="combined")
    #         stacked_image = stack(aligned_images, correct_exposure=True)
    #         write_out_fits_2(stacked_image, "sta/{}_{}_stacked.fits".format(target, band))
    for target in ["m52"]:
        for band in ["r", "g"]:
            unaligned_images = load_fits(path="sci/", target=target, band=band)
            aligned_images = align(unaligned_images, centroid=hybrid_centroid, filter="none")
            stacked_image = stack(aligned_images, correct_exposure=False)
            write_out_fits_2(stacked_image, "sta/{}_{}_stacked.fits".format(target, band))
    for target in ["m52"]:
        for band in ["u"]:
            unaligned_images = load_fits(path="sci/", target=target, band=band)
            aligned_images = align(unaligned_images, centroid=hybrid_centroid, filter="combined")
            stacked_image = stack(aligned_images, correct_exposure=False)
            write_out_fits_2(stacked_image, "sta/{}_{}_stacked.fits".format(target, band))

if __name__ == '__main__':
    main()
