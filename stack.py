"""Script to align and stack images.

"""
import fits_utils as fu

def main():
    for target in ["m52"]:
        for band in ["r","g"]:
<<<<<<< HEAD
            unaligned_images = load_fits(path="sci/", target=target, band=band)
            aligned_images = align(unaligned_images, centroid=hybrid_centroid, filter="none")
            stacked_image = stack(aligned_images, correct_exposure=True)
            write_out_fits(stacked_image, "sta/{}_{}_stacked.fits".format(target, band))
    for target in ["m52"]:
        for band in ["u"]:
            unaligned_images = load_fits(path="sci/", target=target, band=band)
            aligned_images = align(unaligned_images, centroid=hybrid_centroid, filter="combined")
            stacked_image = stack(aligned_images, correct_exposure=True)
            write_out_fits(stacked_image, "sta/{}_{}_stacked.fits".format(target, band))
=======
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
>>>>>>> 1664a866f074d077c935f982d75ba1bf7bb8f70f

if __name__ == '__main__':
    main()
