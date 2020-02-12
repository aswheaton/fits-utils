"""Script to align and stack images.

"""
import fits_utils as fu

def main():
    for target in ["m52"]:
        for band in ["r","g"]:
            unaligned_images = fu.load_fits(path="sci/", target=target, band=band)
            aligned_images = fu.align(unaligned_images, centroid=hybrid_centroid, filter="none")
            stacked_image = fu.stack(aligned_images, correct_exposure=True)
            fu.write_out_fits(stacked_image, "sta/{}_{}_stacked.fits".format(target, band))
    for target in ["m52"]:
        for band in ["u"]:
            unaligned_images = fu.load_fits(path="sci/", target=target, band=band)
            aligned_images = fu.align(unaligned_images, centroid=hybrid_centroid, filter="combined")
            stacked_image = fu.stack(aligned_images, correct_exposure=True)
            fu.write_out_fits(stacked_image, "sta/{}_{}_stacked.fits".format(target, band))

if __name__ == '__main__':
    main()
