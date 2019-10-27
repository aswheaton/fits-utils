#! /usr/bin/env python

"""
reduction.py is used to reduce astro images blah blah blah
"""

# import numpy as np
# from astropy.io import fits
# from os import listdir
# from os.path import isfile, join

import sys

from astropy.io import fits
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

def combine_dark(dark_list):
    """
    combines darks by taking the median of the data
    """

def subtract_dark(dark_frame, raw_frame):
    """
    subtracts a master dark from a fits file depending on the integration time
    """

    return(adjusted_frame)

def normalise_flat():
    """
    divides a flat object by the mode of the data to normalise it
    """

def divide_flat():
    """
    divides a fits object by a flat depending on the integraton time and band
    """

def align_images():
    """
    aligns multiple images of the same object in the same band so they can be
    combined correctly
    """

def combine_images():
    """
    combines multiple images of the same object in the same band to be used for
    science
    """

def main():
    get_fits("dat/bd71_g10sec_001.fits")
main()
