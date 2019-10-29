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

def get_fits(filename):
    """
    function for importing fits files to fits objects
    creates a list of fits objects
    """
    hdul = fits.open(filename)
    hdul.info()
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

def subtract_dark(master_dark_frame, raw_frame):
    """
    subtracts a master dark from a fits file depending on the integration time
    """
    dark_adjusted_frame = raw_frame - master_dark_frame
    return(adjusted_frame)

# def normalise_flat():
#     """
#     divides a flat object by the mode of the data to normalise it
#     """
#
# def divide_flat():
#     """
#     divides a fits object by a flat depending on the integraton time and band
#     """
#
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
    filelist = ['dat/dark_5sec_001.fits',
                'dat/dark_5sec_002.fits',
                'dat/dark_5sec_003.fits',
                'dat/dark_5sec_004.fits',
                'dat/dark_5sec_005.fits'
                ]
    master_dark_frame = average_frame(filelist, mode='astropy', average='median')
    write_out_fits(master_dark_frame, 'master_dark_frame.fits')
main()
