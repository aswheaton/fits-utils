#! /usr/bin/env python

"""
reduction.py is used to reduce astro images blah blah blah
"""

# import numpy as np
# from astropy.io import fits
# from os import listdir
# from os.path import isfile, join

import sys
import math as m

from astropy.io import fits
from astropy.nddata import CCDData
#from pyraf import iraf

# def get_fits(filename):
#     """
#     function for importing fits files to fits objects
#     creates a list of fits objects
#     """
#     hdul = fits.open(filename)
#     hdul.info()
    # file_list = listdir(dir)
    # fits_list = []
    # for item in file_list:
    #     if item.endswith(".fits") == false:
    #         file_list.remove(item)
    #     fits_image_filename = fits_util.get_testdata_filepath(item)
    #     fits_list.append(fits.open(fits_image_filename))

def get_image(filename):
    """
    Creates and returns a CCDData object from a .fits filename.
    """
    image = CCDData.read(filename)
    return(image)

# def subtract_dark(master_dark_frame, raw_frame):
#     """
#     subtracts a master dark from a fits file depending on the integration time
#     """
#     return(raw_frame - master_dark_frame)

# def normalise_flat(flat_frame):
#     """
#     divides a flat object by the mode of the data to normalise it
#     """
#     # Obtain mode of the flat frame


#     return(norm_flat_frame)

# def divide_flat(master_flat_frame, raw_frame):
#     """
#     divides a fits object by a flat depending on the integraton time and band
#     """
#     return(raw_frame)

# def align_images():
#     """
#     aligns multiple images of the same object in the same band so they can be
#     combined correctly
#     """
def median(list):
    """
    Returns the median value of a given list. Used for averaging frames without
    the influence of outlier count values.
    """
    list.sort()
    if len(list)%2 == 0:
        index = int(m.floor(len(list)/2)) - 1
        median = int(m.floor((list[index] + list[index+1])/2))
    if len(list)%2 == 1:
        index = int(m.ceil(len(list)/2))
        median = list[index]
    return(median)

def average_frame(filelist, **kwargs):
    """
    Recieves a list of .fits images and returns either their mean or median as a
    CCDData object, depending upon the value of keyword argument "average=".
    Uses either the astropy or the pyraf module depending upon the value of the
    keyword argument "mode=".

    Utilises the CCDData class; documentation can be found here:
    https://docs.astropy.org/en/stable/io/unified.html
    """
    if kwargs.get('mode') == 'astropy':
        if kwargs.get('average') == 'mean':
            for filename in filelist:
                average += get_image(filename)
            average = average / len(filelist)
        if kwargs.get('average') == 'median':
            for filename in filelist:
                images.append(get_image(filename))
            # Check that all the images are of the same dimensions.
            ra_pix = # ra width of the first CCDData object
            dec_pix = # dec width of the first CCDData object
            for image in images:
                if ra_pix != # ra width of image:
                    # Raise Exception
                if dec_pix != # dec width of image:
                    # Raise Exception
            # Initialise empty CCDData object.
            average = CCDData(np.zeros((ra_pix, dec_pix)))
            for i in range(ra_pix):
                for j in range(dec_pix):
                    counts = []
                    for image in images:
                        counts.append(image[i][j])
                    average[i][j] = median(counts)
        return average
    if kwargs.get('mode') == 'pyraf':
        # lol, don't do this.

        return average

def write_out_fits(image, filename):
    image.write(filename)

def main():
    # Example dark frame averaging.
    dark_list = ['dat/dark_5sec_001.fits',
                'dat/dark_5sec_002.fits',
                'dat/dark_5sec_003.fits',
                'dat/dark_5sec_004.fits',
                'dat/dark_5sec_005.fits'
                ]
    # Example flat frame averaging
    flat_list = ['dat/Flat_g5sec_001.fits',
                'dat/Flat_g5sec_002.fits',
                'dat/Flat_g5sec_003.fits',
                'dat/Flat_g5sec_004.fits',
                'dat/Flat_g5sec_005.fits',
                'dat/Flat_g5sec_006.fits',
                'dat/Flat_g5sec_007.fits',
                ]
    #Example raw list
    raw_list = ['dat/bd71_g10sec_001.fits',
                'dat/bd71_g10sec_002.fits',
                'dat/bd71_g10sec_003.fits'
                ]

    master_dark_frame = average_frame(dark_list, mode='astropy', average='median')
    master_flat_frame = average_frame(filelist, mode='astropy', average='median')

    # Empty list of sciences
    science_list = []
    for raw_filename in raw_list:
        # Get CCDData from fits file
        raw_image = get_image(raw_filename)
        # Dark subtract
        raw_image.subtract(master_dark_frame)
        # Flat divide
        raw_image.divide(master_flat_frame)
        # Append to science list
        science_list.append(raw_image)

    write_out_fits(master_dark_frame, 'master_dark_frame.fits')
main()
