#! /usr/bin/env python

"""
reduction.py is used to reduce astro images blah blah blah
"""

# import numpy as np
# from astropy.io import fits

import sys
import math as m
import numpy as np

from astropy.io import fits
from astropy.nddata import CCDData
#from pyraf import iraf

from os import walk
import configparser

def get_lists(dir):
    """

    Function to create lists of darks, flats, targets, and standard stars. Each
    list contains a dictionary and is constructed by the add_to_list method.
    The configparser package is used to read a .ini file containing settings
    that may be changed by the user in a text editor.

    The function loops through the specified directory and assigns .fits files
    to the lists if they match data in the config.ini file.

    Args:
        dir (str): The path of the directory to be used.

    Returns:
        dark_list (list): The list of 'dark' .fits files.
        flat_list (list): The list of 'flat' .fits files.
        target_list (list): The list of 'target' .fits files.
        standard_star_list (list): The list of 'standard stars' .fits files.

    """
    dark_list = []
    flat_list = []
    target_list = []
    standard_star_list = []

    def add_to_list(array, filename, **kwargs):
        """

        Function to construct a dictionary for a fits file containing a filename
        key and keys for other optional information: integration time and band.
        The dictionary is appended to a specific list.

        Args:
            array (list): The list to append the dictionary to.
            filename (str): The name of the fits file.
            **kwargs: Arbitary keyword arguments.

        """
        new_dict = {
            'filename': filename
        }
        for key, value in kwargs.items():
             new_dict[key] = value

        array.append(new_dict)

    config = configparser.ConfigParser()
    config.read('config.ini')
    data_settings = config['DATA SETTINGS']

    _, _, filenames = next(walk(dir), (None, None, []))
    for filename in filenames:
        if filename.endswith('.fits'):
            name_str = filename.split('_')
            if name_str[0].lower() == 'dark':
                add_to_list(dark_list, filename, integration_time = name_str[1])
            elif name_str[0].lower() == 'flat':
                add_to_list(flat_list, filename, intgeration_time = name_str[2], band = name_str[1])
            elif name_str[0].lower() == data_settings['standard star a'] or name_str[0].lower() == data_settings['standard star b']:
                add_to_list(standard_star_list, filename, integration_time = name_str[2], band = name_str[1])
            elif name_str[0].lower() == data_settings['target id']:
                add_to_list(target_list, filename, integration_time = name_str[2], band = name_str[1])

    return dark_list, flat_list, target_list, standard_star_list

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
        images = []
        if kwargs.get('average') == 'median':
            for filename in filelist:
                images.append(get_image(filename))
            # Check that all the images are of the same dimensions.
            ra_pix = len(images[0])
            dec_pix = len(images[0][0])
            for image in images:
                if ra_pix != len(image):
                    print('Warning! Image dimensions do not match!\n')
                    break
                if dec_pix != len(image[0]):
                    print('Warning! Image dimensions do not match!\n')
                    break
            # Initialise empty array and populates it with median values.
            average = np.zeros((ra_pix, dec_pix))
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
    hdu = fits.PrimaryHDU(image)
    hdul = fits.HDUList([hdu])
    hdul.writeto(filename)

def main():
    # Create lists of darks/flats/raws for specific integration times and bands.
    raw_dark_list, raw_flat_list, raw_target_list, raw_stand_star_list = get_lists('dat/')

    # Create a new list of normalised flats in the tmp/ directory.
    norm_flat_list = normalise_flats(flat_list)

    # Create master dark and flat frame in the tmp/ directory.
    master_dark_frame = average_frame(dark_list, mode='astropy', average='median')
    master_flat_frame = average_frame(norm_flat_list, mode='astropy', average='median')
    write_out_fits(master_dark_frame, 'master_dark_frame.fits')
    write_out_fits(master_flat_frame, 'master_flat_frame.fits')

    # Create a list of reduced science frames for alignment.
    science_list = []
     for raw_filename in raw_list:
         raw_image = get_image(raw_filename)
         raw_image.subtract(master_dark_frame)
         raw_image.divide(master_flat_frame)
         science_list.append(raw_image)

    # Align images.
    aligned_science_list = align_images(science_list)

    # Stack the frames and write out.
    stacked_science_frame = stack_images(aligned_science_list)
    write_out_fits(stacked_science_frame, filename)

main()
