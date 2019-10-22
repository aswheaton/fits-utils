"""
reduction.py is used to reduce astro images blah blah blah
"""

import numpy as np
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
from astropy.nddata import Cutout2D
from astropy.nddata import CCDData

def get_fits():
    """
    function for importing fits files to fits objects
    """

def combine_dark():
    """
    combines darks by taking the median of the data
    """

def subtract_dark():
    """
    subtracts a master dark from a fits file depending on the integration time
    """

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
    aligns multiple images of the same object in the same band so they can be combined correctly
    """

def combine_images():
    """
    combines multiple images of the same object in the same band to be used for science
    """

def main():
    """
    Main wrapper function
    """

if __name__ == "__main__": main()
