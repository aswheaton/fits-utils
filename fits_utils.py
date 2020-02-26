#! /usr/bin/env python

"""
fits-utils module

This module contains various functions for the reduction, alignment and stacking
of raw astronomical frames using the FITS file formatting standard. It is an
extension of the astropy module, and utilises many of its basic functions to
perform more complicated operations.

The raw files are stored in the "dat/" folder, which is read-only; science
frames are written to and read from in the "sci/" folder; temporary files
(including master dark and flat fields) are written to and read from the "tmp/"
directory; and stacked images are written to and read from the "sta/" directory.

This module utilises the configparser module in order to generate a "config.ini"
file. This file can be populated by the user with a target ID, standard star ID,
and observing bands. The module can then use these parameters to parse the
"raw/" directory and sort FITS files into lists with specific integration times
and bands. The images can then be reduced and written out to FITS files.

By convention, imported:

import fits-utils as fu
"""

import numpy as np
import matplotlib.pyplot as plt
import configparser

from matplotlib.colors import LogNorm
from scipy.stats import mode
from scipy.ndimage import gaussian_filter
from pathlib import Path
from astropy.io import fits
from os import walk

def gen_config():
    config = configparser.ConfigParser()
    config["TELESCOPE"] = {"Size": "15"}
    config["DATA SETTINGS"] = {"Standard Star A": "bd62",
                               "Standard Star B": "bd25",
                               "Target ID": "m52",
                               "Bands": "g, r, u"}
    with open("config.ini", "w") as configfile:
        config.write(configfile)

def get_lists(dir):
    """
    Sorts fits files into lists.
    Function to create lists of darks, flats, targets, and standard stars. Each
    list contains a dictionary and is constructed by the add_to_list method.
    The configparser package is used to read a .ini file containing settings
    that may be changed by the user in a text editor.

    The function loops through the specified directory and assigns .fits files
    to the lists if they match data in the config.ini file.

    Args:
        dir (str): The path of the directory to be used.

    Returns:
        dark_list (list): The list of "dark" .fits files.
        flat_list (list): The list of "flat" .fits files.
        target_list (list): The list of "target" .fits files.
        standard_star_list (list): The list of "standard stars" .fits files.
    """
    #: list: Initially empty lists for containing dicts
    dark_list = []
    flat_list = []
    target_list = []
    standard_star_list = []
    def add_to_list(array, filename, **kwargs):
        """
        Adds file data to lists of dictionaries.

        Function to construct a dictionary for a fits file containing a filename
        key and keys for other optional information: integration time and band.
        The dictionary is appended to a specific list.

        Args:
            array (list): The list to append the dictionary to.
            filename (str): The name of the fits file.
            **kwargs: Arbitrary keyword arguments.
        """
        #: dict of str: Holds keys and values for filename and optional keys
        new_dict = {
            "filename": filename
        }
        for key, value in kwargs.items():
             new_dict[key] = value
        array.append(new_dict)
    #: ConfigParser: Contains setup data stored in .ini file.
    config = configparser.ConfigParser()
    config.read("config.ini")
    data_settings = config["DATA SETTINGS"]
    _, _, filenames = next(walk(dir), (None, None, []))
    for filename in filenames:
        if filename.endswith(".fits"):
            name_str = filename.split("_")
            if name_str[0].lower() == "dark":
                add_to_list(dark_list, filename, integration_time = name_str[1])
            elif name_str[0].lower() == "flat":
                add_to_list(flat_list, filename, integration_time = name_str[2], band = name_str[1])
            elif name_str[0].lower() == data_settings["standard star a"] or name_str[0].lower() == data_settings["standard star b"]:
                add_to_list(standard_star_list, filename, integration_time = name_str[2], band = name_str[1])
            elif name_str[0].lower() == data_settings["target id"]:
                add_to_list(target_list, filename, integration_time = name_str[2], band = name_str[1])
    return dark_list, flat_list, target_list, standard_star_list

def load_fits(**kwargs):
    """
    Receives a directory path and .fits filename parameters. Parses the
    directory for files matching the naming parameters and loads matched files
    into a list of astropy.fits objects. Returns the list.
    """
    path = kwargs.get("path")
    # year = kwargs.get("year")
    band = kwargs.get("band")
    target_id = kwargs.get("target")
    images = []
    for root, dir, files in walk(path):
        for filename in files:
            if "left" in filename:
                pass
            elif "at" in filename:
                pass
            elif "end" in filename:
                pass
            elif target_id in filename and band in filename:
                int_time = int(filename.split("_")[2][:-1])
                print(target_id, band, int_time, " matched ", filename)
                with fits.open(root+filename) as hdul:
                    new_image = {}
                    new_image["int_time"] = int_time
                    new_image["target"] = target_id
                    new_image["filename"] = filename
                    new_image["data"] = hdul[0].data
                    images.append(new_image)
    # Check that all the sizes of the loaded fits objects match.
    x_sizes, y_sizes = [], []
    for image in images:
        x_sizes.append(image["data"].shape[0])
        y_sizes.append(image["data"].shape[1])
    # If all fits object dimensions match, return the existing list of objects.
    if all(x==x_sizes[0] for x in x_sizes) and all(y==y_sizes[0] for y in y_sizes):
        return(images)
    # If not, cast them into arrays of zeros with matching dimensions and copy
    # the header data over to the newly created fits objects.
    else:
        print("Imported image dimensions do not match! Framing with zeros.")
        framed_images = []
        for image in images:
            framed_image = {}
            framed_image["data"] = np.zeros((max(x_sizes),max(y_sizes)))
            framed_image["data"][:image["data"].shape[0],:image["data"].shape[1]] = image["data"]
            framed_image["int_time"] = image["int_time"]
            framed_image["target"] = image["target"]
            framed_image["filename"] = image["filename"]
            framed_images.append(framed_image)
        return(framed_images)

def average_frame(filelist, **kwargs):
    """
    Recieves a list of .fits images and returns either their mean or median as
    a HDUList object, depending upon the value of keyword argument "average=".
    Uses either the astropy or the pyraf module depending upon the value of the
    keyword argument "mode=".

    Args:
        filelist (list): list of HDUList.
        **kwargs: Arbitrary keyword arguments.
    Returns:
        average (HDUList): object containing HDUList that can be written to a
            file.
    """
    if kwargs.get("mode") == "astropy":
        if kwargs.get("average") == "mean":
            for file_object in filelist:
                average += file_object
            average = average / len(filelist)
        if kwargs.get("average") == "median":
            # Check that all the images are of the same dimensions.
            ra_pix = len(filelist[0])
            dec_pix = len(filelist[0][0])
            for file_object in filelist:
                if ra_pix != len(file_object):
                    print("Warning! Image dimensions do not match!\n")
                    break
                if dec_pix != len(file_object[0]):
                    print("Warning! Image dimensions do not match!\n")
                    break
            # Initialise empty array and populates it with median values.
            average = np.zeros((ra_pix, dec_pix))
            for i in range(ra_pix):
                for j in range(dec_pix):
                    counts = []
                    for file_object in filelist:
                        counts.append(file_object[i][j])
                    average[i][j] = median(counts)
        return average

def write_out_fits(image, filename):
    """
    Creates a header for an ndarray of reduced data and then creates a new fits
    file of this data.

    Args:
        image (ndarray): reduced data to be written to fits file.
        filename (string): name (and location) of new fits file.
    """
    hdul = fits.HDUList([fits.PrimaryHDU(image)])
    hdul.writeto(filename, overwrite=True)

def write_out_fits_2(image, filename):
    """
    Creates a header for an ndarray of reduced data and then creates a new fits
    file of this data.

    Args:
        image (ndarray): reduced data to be written to fits file.
        filename (string): name (and location) of new fits file.
    """
    hdul = fits.HDUList([fits.PrimaryHDU(image["data"])])
    hdul.writeto(filename, overwrite=True)

def normalise_flat(flat_array):
    """
    Normalises the data in a flat frame by dividing each data value by the modal
    value.

    Args:
        flat_array (ndarray): flat data to be normalised.
    Returns:
        normalised_flat (ndarray): new normalised flat data.
    """
    normalised_flat = flat_array / mode(flat_array, axis=None)[0][0]
    return normalised_flat

def max_value_centroid(image_data, **kwargs):
    """
    Receives an image array and returns the coordinates of the brightest pixel
    in that image array.

    Args:
        image_data (2darray): The image array to be searched.
    Returns:
        (x_max,y_max) (tuple): pixel coordinates of the maximum value of the
            image array.
    """
    x_max, y_max = np.where(image_data == np.amax(image_data))
    return((x_max[0], y_max[0]))

def custom_roll(array, axis=0):
    """
    Getting sum of nearest neighbours for each value in an array.

    For each value in an ndarray, sums the nearest neighbours around that
    position for one axis. Pads array with zeroes, aligns array with dstack and
    roll.

    Args:
        array (ndarray): Array to be passed over.
        axis (int): Specific axis to be rolled over. Defaults to 0.
    """
    no = 3
    array = array.T if axis==1 else array
    padding = np.zeros((no-1, array.shape[1]))
    array = np.concatenate([array, padding], axis=0)
    neighbours = np.dstack([np.roll(array, i, axis=0) for i in range(no)])
    array = neighbours.sum(2)[1:-1, :]
    array = array.T if axis==1 else array
    return array

def create_mask(image_data, **kwargs):
    # Offset image so that all values are positive
    offset_data = image_data + np.abs(np.amin(image_data))
    mask = np.empty(offset_data.shape)
    if kwargs.get("condition") == "neighbors":
        sum_of_neighbours = custom_roll(custom_roll(offset_data), axis=1) - offset_data
        mask = offset_data > sum_of_neighbours
    # Invalidate values that fall below a certain threshold (fast).
    if kwargs.get("condition") == "threshold":
        max_value = np.amax(offset_data)
        for i in range(offset_data.shape[0]):
            for j in range(offset_data.shape[1]):
                if offset_data[i,j] <= 0.67 * max_value:
                    mask[i,j] = 1
                else:
                    mask[i,j] = 0
    if "border" in list(kwargs.keys()):
        size = kwargs.get("border")
        mask[:size,:] = 1
        mask[-size:,:] = 1
        mask[:,:size] = 1
        mask[:,-size:] = 1
    return(mask)

def smooth(image_data, **kwargs):
    sigma = kwargs.get("sigma")
    smoothed_image = gaussian_filter(image_data, sigma)
    return(smoothed_image)

def weighted_mean_2D(cutout,**kwargs):
    """
    Recieves an argument of type ndarray and returns a tuple of the weighted
    mean centroid of the object contained in the cutout.

    Args:
        cutout (ndarray): portion of fits full ndarray.
    Returns:
        x_avg (int): weighted mean of x values.
        y_avg (int): weighted mean of y values.
    """
    cutout += np.abs(np.amin(cutout))
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    x_avg = np.average(range(x_sum.size), weights=x_sum)
    y_avg = np.average(range(y_sum.size), weights=y_sum)
    if kwargs.get("floor") == True:
        return((int(np.floor(x_avg)), int(np.floor(y_avg))))
    else:
        return((x_avg, y_avg))

def hybrid_centroid(image_data, **kwargs):
    """
    Recieves an array of image data and returns the pixel coordinates of the
    centroid of the brightest star in the frame. Makes an initial guess at the
    position of the star by finding the maximum value in the array, then
    performs a weighted mean in two dimensions about the guess for finer accuracy.

    Args:
        image_data (2darray): array of image data containing reference star.
        size (int): the radius of the reference star, in pixels. Used to create
            cutout of appropriate size.
    Returns:
        (x_avg,y_avg) (tuple): pixel coordinates of the centroid of the
            brightest star in the image array.
    """
    size = kwargs.get("size")
    # Attempt to invalidate pixels which may confuse the initial guess.
    if kwargs.get("filter") == "mask":
        mask_array = create_mask(image_data, condition="neighbors", border=size)
        masked_data = np.ma.array(image_data, mask=mask_array)
        x_max, y_max = max_value_centroid(masked_data)
    # Attempt to smooth out pixels which may confuse the initial guess.
    elif kwargs.get("filter") == "gaussian":
        x_max, y_max = max_value_centroid(smooth(image_data, sigma=0.25*size))
    # A hybrid method for aligning very faint images.
    elif kwargs.get("filter") == "combined":
        smoothed_data = smooth(image_data, sigma=0.25*size)
        mask_array = create_mask(smoothed_data, condition="neighbors", border=100)
        masked_data = np.ma.array(smoothed_data, mask=mask_array)
        x_max, y_max = max_value_centroid(masked_data)
    # Get the maximum value of the cutout as an initial guess.
    else:
        x_max, y_max = max_value_centroid(image_data)
    # Create a smaller cutout around the initial guess.
    cutout = np.array(image_data[
                x_max-size if x_max > size else 0:x_max+size if x_max+size < image_data.shape[0] else image_data.shape[0],
                y_max-size if y_max > size else 0:y_max+size if y_max+size < image_data.shape[1] else image_data.shape[1]
                ])
    # Get the mean weighted average of the smaller cutout.
    x_new, y_new = weighted_mean_2D(cutout, floor=True)
    # Map the centroid back to coordinates of original cutout.
    x_avg = x_max - size + x_new
    y_avg = y_max - size + y_new
    return((x_avg, y_avg))

def align(images, **kwargs):
    """
    Recieves a list of image arrays containing a common object to use for
    alignment of the image stack. Returns a list of image arrays of different
    size, aligned, and with zero borders where the image has been shifted.

    Args:
        images (list of dict): Frames to be aligned.

    Returns:
        aligned_images (list of dict): new frames that have been aligned and can
            be stacked.
    """
    # Which centroiding function to use.
    centroid_func = kwargs.get("centroid")
    # Boolean, whether or not to mask images for hot pixels on the detector.
    filter = kwargs.get("filter")
    # Find the centroid of the reference star in each image.
    x_centroids, y_centroids = [], []
    print("---Beginning Alignment---")
    counter = 0
    for image in images:
        counter += 1
        print("---Finding Centre {} of {}".format(counter, len(images)), end="\r")
        centroid = centroid_func(image["data"], size=50, filter=filter)
        x_centroids.append(centroid[0])
        y_centroids.append(centroid[1])
        image["XCENT"] = centroid[0]
        image["YCENT"] = centroid[1]

        # if counter == 1:
        #
        #     fig1 = plt.imshow(smooth(image["data"], sigma=3), origin="lower", cmap="viridis", norm=LogNorm())
        #     plt.scatter(centroid[1], centroid[0], s=1, c="red", marker="o")
        #     import matplotlib.patches as patches
        #     rect = patches.Rectangle((centroid[1]-30,centroid[0]-30),60,60,linewidth=1,edgecolor="r",facecolor="none")
        #     fig1.axes.add_patch(rect)
        #     fig1.axes.get_xaxis().set_visible(False)
        #     fig1.axes.get_yaxis().set_visible(False)
        #     plt.savefig("r_hybrid_centroid_full.jpeg", bbox_inches="tight", pad_inches=0, dpi=1000)
        #
        #     small_image = np.array(image["data"][centroid[0]-30:centroid[0]+30,centroid[1]-30:centroid[1]+30])
        #     fig2 = plt.imshow(smooth(small_image, sigma=3), origin="lower", cmap="viridis", norm=LogNorm())
        #     plt.scatter(30,30, s=1, c="red", marker="o")
        #     fig2.axes.get_xaxis().set_visible(False)
        #     fig2.axes.get_yaxis().set_visible(False)
        #     plt.savefig("r_hybrid_centroid_small.jpeg", bbox_inches="tight", pad_inches=0, dpi=1000)
        #
        #     small_image = np.array(image["data"][centroid[0]-30:centroid[0]+30,centroid[1]-30:centroid[1]+30])
        #     margin = [np.sum(small_image, axis=1)]
        #     fig2 = plt.imshow(margin, origin="lower", cmap="viridis", norm=LogNorm())
        #     fig2.axes.get_xaxis().set_visible(False)
        #     fig2.axes.get_yaxis().set_visible(False)
        #     plt.savefig("r_max_margin_1.jpeg", bbox_inches="tight", pad_inches=0, dpi=1000)

    print()
    max_pos_x, max_pos_y = max(x_centroids), max(y_centroids)
    min_pos_x, min_pos_y = min(x_centroids), min(y_centroids)
    max_dif_x, max_dif_y = max_pos_x - min_pos_x, max_pos_y - min_pos_y
    # Create new stack of aligned images using the centroid in each frame.
    aligned_images = []
    counter = 0
    for image in images:
        counter += 1
        print("---Aligning Image {} of {}".format(counter, len(images)), end="\r")
        # Determine region in which to cast the image.
        disp_x, disp_y = max_pos_x - image["XCENT"], max_pos_y - image["YCENT"]
        imsize_x, imsize_y = image["data"].shape[0], image["data"].shape[1]
        # Create new array containing aligned image data.
        aligned_image_data = np.zeros((imsize_x+max_dif_x, imsize_y+max_dif_y), dtype=int)
        aligned_image_data[disp_x:disp_x+imsize_x,disp_y:disp_y+imsize_y] = image["data"]
        # Create new dictionary value of frames per pixel column (one everywhere for unstacked image).
        frame_count = np.zeros((imsize_x+max_dif_x, imsize_y+max_dif_y), dtype=int)
        frame_count[disp_x:disp_x+imsize_x,disp_y:disp_y+imsize_y] = 1
        # Create new image dictionary and copy over header data from image.
        aligned_image = {"target"      : image["target"],
                         "filename"    : image["filename"],
                         "int_time"    : image["int_time"],
                         "data"        : aligned_image_data,
                         "frame_count" : frame_count
                         }
        # Add the new aligned image dictionary to a list to be returned.
        aligned_images.append(aligned_image)
    print("---Alignment Complete---")
    return(aligned_images)

def stack(aligned_image_stack, **kwargs):
    """
    Receives a list of aligned images and returns their summation along the axis
    of the list.

    Args:
        aligned_image_stack (list of dict): aligned frames ready to be stacked.
    Returns:
        stacked_image (dict): new combined single frame.
    """

    print("---Stacking Images---")

    # Check that the aligned images to be stacked have matching dimensions.
    for image in aligned_image_stack:
        if image["data"].shape != aligned_image_stack[0]["data"].shape:
            print("Aligned image dimensions do not match!")
            break

    # Initialise an empty array into which aligned images are stacked.
    stacked_image_data = np.zeros(aligned_image_stack[0]["data"].shape)
    counter = 0

    if kwargs.get("correct_exposure") == True:

        # Also initialise an empty array into which frames per pixel column are summed.
        total_frame_count = np.zeros(aligned_image_stack[0]["data"].shape)

        for image in aligned_image_stack:
            counter += 1
            print("---Stacking Image {} of {}".format(counter, len(aligned_image_stack)), end="\r")
            # Correct the image data for exposure and sum up the fluxes.
            stacked_image_data += image["data"] / image["int_time"]
            # Sum up the number of overlapping frames per pixel column.
            total_frame_count += image["frame_count"]
        # Replace zeros with ones to prevent zero division error (dodgy!)
        total_frame_count[np.where(total_frame_count == 0)] = 1
        # Average each pixel column over the number of frames in that column.
        stacked_image_data = np.floor(stacked_image_data / total_frame_count)
        stacked_image_data.astype(dtype=int, copy=False)
        # Create new dictionary containg the average flux (counts/second/pixel).
        exposure_corrected_stack = {"data" : stacked_image_data}

        print("---Stacking Complete---")

        return(exposure_corrected_stack)
    else:
        for image in aligned_image_stack:
            stacked_image_data += image["data"]
        stacked_image = {"data" : stacked_image_data}
        return(stacked_image)

def rgb(image_r, image_g, image_b):
    """
        Recieves three arrays of equal size. Maps these values to RGB values
        using the Lupton algorithm and displays the resulting image.
        # TODO: Retrieve the source for this algorithm.
    """
    from astropy.visualization import make_lupton_rgb
    rgb_image = make_lupton_rgb(image_r, image_g, image_b, Q=10, stretch=1000.)
    plt.imshow(rgb_image, norm=LogNorm())
    plt.show()
    plt.axis("off")
    fig = plt.imshow(rgb_image, norm=LogNorm())
    fig.axes.get_xaxis().set_visible(False)
    fig.axes.get_yaxis().set_visible(False)
    plt.savefig("false_colour.jpeg", bbox_inches="tight", pad_inches=0, dpi=1000)

def reduce_raws(raw_list, master_dark_frame, master_flat_frame, dir):
    """
    Reduces raw images into science images. This function loops through a
    list of raws. Each raw image is dark subtracted and then flat divided.

    Args:
        raw_list (list): Raw ndarray objects.
        master_dark_frame (dict): Dark ndarrays.
        master_flat_frame (dict): Flat ndarrays.
        dir (directory): Location of .fits files to be reduced.
    Returns:
        science_list (list): Reduced ndarray objects.
    """
    #: list of ndarray: Empty list for reduced images.
    science_list = {}
    for raw in raw_list:
        print("Reducing {} of {} images.".format(len(science_list), len(raw_list)), end="\r"),
        with fits.open(dir / raw["filename"]) as hdul:
            #: ndarray: Dark subtracted image data
            ds_data = np.subtract(hdul[0].data, master_dark_frame[raw["integration_time"]])
            science_list[raw["filename"]] = np.divide(ds_data, master_flat_frame[raw["band"]])
    print("\nDone!")
    return science_list

def trim(filename):
    hdul = fits.open(filename)
    hdul[0].data = np.array(hdul[0].data[50:-50,50:-50])
    hdul.writeto(filename, overwrite=True)

def get_zero_points(input_airmass):

    # bd62_standard_mags = {"r":9.332, "g":9.872, "u":11.44}
    # bd25_standard_mags = {"r":9.929, "g":9.450, "u":9.023}

    for band in ["r","g","u"]:
        counts_and_errs = np.loadtxt("standard_stars/std_{}.csv".format(band))
        airmasses = np.array(counts_and_errs[:,3])
        zero_points = counts_and_errs[:,0] + 2.5 * np.log10(counts_and_errs[:,1])
        zero_point_errs = counts_and_errs[:,2] * 2.5 / counts_and_errs[:,1] / np.log(10)

        gradient = np.sum((airmasses-np.mean(airmasses))*(zero_points-np.mean(zero_points))) / np.sum((airmasses-np.mean(airmasses))**2)
        intercept = np.mean(zero_points) - gradient * np.mean(airmasses)
        zero_point = gradient * input_airmass + intercept

        if band == "r": zpr = zero_point
        elif band == "g": zpg = zero_point
        elif band == "u": zpu = zero_point

        # plt.plot(airmasses, zero_points, 'o')
        # plt.plot(airmasses, gradient*airmasses+intercept, '-')
        # plt.title('Zero Point vs Airmass ({}-band)'.format(band))
        # plt.show()

    return(zpr, zpg, zpu)

def correct_pleiades(p_data):
    """
    Converts from Johnson/Cousins system to SDSS system, using transformations
    from Jester et al. (2005): https://sdss3.org/dr8/algorithms/sdssUBVRITransform.php

    Also de-reddens pleiades data using literature values of absorption from NEDS:
    http://ned.ipac.caltech.edu/cgi-bin/nph-objsearch?objname=MESSIER%2045&extend=no&out_csys=Equatorial&out_equinox=J2000.0&obj_sort=RA+or+Longitude&of=pre_text&zv_breaker=30000.0&list_limit=5&img_stamp=YES

    [0] : g-r
    [1] : u-g
    [2] : r
    """
    # Convert colour excess from Johnson U-B, B-V to Sloan u-g, g-r.
    #: B-V --> g-r
    p_data[:,0] = 1.02 * p_data[:,0] - 0.22
    #: U-B --> u-g
    p_data[:,1] = 1.28 * p_data[:,1] + 1.13
    #: V --> r
    p_data[:,2] = p_data[:,2] - 0.46 * p_data[:,0] + 0.11
    # De-redden the converted data using transformations from NED.
    p_data[:,0] = p_data[:,0] - 1.009 + 0.787
    p_data[:,1] = p_data[:,1] - 0.787 + 0.544
    p_data[:,2] = p_data[:,2] - 0.544
    return(p_data)

def get_mag(flux, flux_error, zero_point):
    """
    Recieves the flux of an object, the error on that flux, and an instrumental
    zero point and returns a zero point corrected magnitude and magnitude error
    for the object.
    """
    mag = -2.5 * np.log10(flux) + zero_point
    mag_err = (-2.5/flux/np.log(10)) * flux_error
    return(mag, mag_err)

def load_cat(filename, zpr, zpg, zpu):
    """
    Loads in a catalogue output by Source Extractor and returns a numpy arrays
    of calatogue fluxes and their errors after zero point correction.
    """
    catalog = np.loadtxt(filename)
    r_mag, r_err = get_mag(catalog[:,5], catalog[:,6], zpr)
    g_mag, g_err = get_mag(catalog[:,3], catalog[:,4], zpg)
    u_mag, u_err = get_mag(catalog[:,7], catalog[:,8], zpu)
    return(r_mag, r_err, g_mag, g_err, u_mag, u_err)

def write_cat(r_mag, g_mag, u_mag, filename):
    """
    Method for creating a new catalog file from an ndarray.

    Takes an ndarray and writes it out to a new catalog file with a ".cat"
    extension. Creates a header for this file with the indexes included.

    Args:
        catalog (ndarray): Merged catalog data.
        filename (str): Filename of the new catalog. Do not include the
            extension, this is added automatically.
    """
    catalog = np.stack((r_mag, g_mag, u_mag),axis=1)
    #: str: New header text for the output file.
    header_txt = "\n".join(["[0] : NUMBER",
                            "[1] : ALPHAPEAK_J2000",
                            "[2] : DELTAPEAK_J2000",
                            "[3] : FLUX_APER_G",
                            "[4] : FLUXERR_APER_G",
                            "[5] : FLUX_APER_R",
                            "[6] : FLUXERR_APER_R",
                            "[7] : FLUX_APER_U",
                            "[8] : FLUXERR_APER_U"])
    np.savetxt("cat/{}.cat".format(filename), catalog, header=header_txt)

def polynomial(x, coeffs):
    """
    Returns the corresponding y value for x, given the coefficients for an
    nth-order polynomial as a list descending in order.
    """
    # Hard code for 4th order, better to use general case.
    # y = A*x**4 + B*x**3 + C*x**2 + D*x**1 + E*x**0)
    y = 0.0
    for n in range(len(coeffs)):
        y += coeffs[n] * x ** n
    return(y)

def get_r(red_x, red_y, hyp_x, hyp_y, func, coeffs):
    # Slope of the DE-reddening vector.
    slope = (hyp_y - red_y) / (hyp_x - red_x)
    x_vals = np.linspace(red_x, hyp_x, 1000)
    y_val_vec = red_y - slope * (red_x - x_vals)
    y_val_cur = polynomial(x_vals, coeffs)
    y_diffs = abs(y_val_vec - y_val_cur)
    index = np.where(y_diffs == np.amin(y_diffs))[0]
    x_int, y_int = x_vals[index], y_val_cur[index]
    r = ((hyp_x-x_int)**2 + (hyp_y-y_int)**2)**0.5
    return(r)

def get_chi_squ(x, y, func, coeffs, error):
    """
    Returns the chi-squared value for a given x,y dataset and a function handle
    that returns the predicted predicted values.
    """
    chi_squ_tot = np.sum(((y - func(x, coeffs)) / error)**2)
    return(chi_squ_tot)

def get_sumsqu_res(x, y, func, coeffs):
    """
    Returns the sum of squared residuals for a given x,y dataset and a function
    handle that returns the predicted predicted values.
    """
    sum_squ_res = np.sum((y - func(x, coeffs))**2)
    return(sum_squ_res)

def remove_outliers(x, lower_bound, upper_bound):
    """
    Returns the indices of a reduced data set within a given lower and upper bound.
    """
    indices = np.where((x >= lower_bound) & (x <= upper_bound))[0]
    return(indices)

def minimiser(array):
    """
    Returns the index of the minimum value in a 1d-array.
    """
    index = np.where(array==np.amin(array))[0]
    return(index)

def cardelli_a(x):
    y = x - 1.82
    spam = (1 + 0.17699 * y
            - 0.50447 * y**2
            -0.02427 * y**3
            + 0.72085 * y**4
            + 0.01979 * y**5
            - 0.77530 * y**6
            + 0.32999 * y**7
            )
    return(spam)

def cardelli_b(x):
    y = x - 1.82
    eggs = (1.41338 * y
            + 2.28305 * y**2
            + 1.07233 * y**3
            - 5.38434 * y**4
            - 0.62251 * y**5
            + 5.30260 * y**6
            - 2.09002 * y**7
            )
    return(eggs)

def cardelli_const(not_gamma):
    R_v = 3.1 # Ratio of selective to total extinction
    const = cardelli_a(1/not_gamma) + cardelli_b(1/not_gamma) / R_v
    return(const)

def get_cardelli_slope(c_constants):
    return((c_constants["u"]-c_constants["g"])/(c_constants["g"]-c_constants["r"]))

def get_spectral_types(sources):
    pass

def plot_diagram(plts, **kwargs):
    """
    Method for plotting a HR diagram. Uses matplotlib to create a HR diagram
    of the magnitudes/colors. The plot is a simple scatter plot. Saves the
    plot to the plots output folder.

    Args:
        plts (dictionary): Multiple plots to be put on the same axis. Key should
            be the label for the plot. Each tuple should contain: (x values, y
            values, and the marker type for the plot.)
    """
    #: fig, ax objects: New figure and axis created by matplotlib.
    fig, ax = plt.subplots()
    plt_names = []
    for plt_name, plot in plts.items():
        ax.plot(plot[0], plot[1], plot[2], markersize=0.75)
        plt_names.append(plt_name)
    ax.set(
           xlabel=kwargs.get("x_label"),
           ylabel=kwargs.get("y_label"),
           title=kwargs.get("sup_title"),
           #: Invert the y axis for the plot.
           ylim=ax.get_ylim()[::-1]
          )
    if kwargs.get("legend")==True:
        ax.legend(plt_names)
    plt.draw()
    plt.show()
    if kwargs.get("filename")!=None:
        fig.savefig("plots/{}.jpeg".format(kwargs.get("filename")), dpi=1000)
