#! /usr/bin/env python

"""Reduction program

This program is used to reduce fits files so that they are ready to be used for
science. The raw files are stored in the "dat/" folder. The program utilises the
configparser module in order to read a "config.ini" file. This file contains a
target ID, standard star ID, and observing bands. The program sorts the data
into lists of specific integration times and bands. The files are then reduced.
The program outputs reduced fits files.

"""
import math as m
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import mode
from pathlib import Path
from astropy.io import fits
from os import walk
import configparser

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
    # Why is this function defined inside another? ~Wheaton
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
            elif "attempt" in filename:
                pass
            elif target_id in filename and band in filename:
                int_time = filename.split("_")[2][:-1]
                print(target_id, band, int_time, " matched ", filename)
                hdul = fits.open(root + filename)
                hdul[0].header["EXPTIME"] = int_time
                hdul[0].header["TARGTID"] = target_id
                images.append(fits.open(root + filename))
    return(images)

def median(list):
    """
    DEPRECATED
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
    if kwargs.get("mode") == "pyraf":
        # lol, don"t do this.
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
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    x_avg = np.average(range(x_sum.size), weights=x_sum)
    y_avg = np.average(range(y_sum.size), weights=y_sum)
    # plt.imshow(cutout)
    # plt.scatter(x_avg, y_avg, s=2, c='red', marker='o')
    # plt.show()
    if kwargs.get("floor") == True:
        return((int(np.floor(x_avg)), int(np.floor(y_avg))))
    else:
        return((x_avg, y_avg))

def max_value_centroid(cutout, **kwargs):
    """
    Returns the coordinates of the brightest pixel in a cutout of an array.
    """
    cutout += np.abs(np.min(cutout))
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    x_max = np.where(x_sum == max(x_sum))[0]
    y_max = np.where(y_sum == max(y_sum))[0]
    # plt.imshow(cutout)
    # plt.scatter(x_max, y_max, s=2, c='red', marker='o')
    # plt.show()
    return((int(np.floor(x_max)), int(np.floor(y_max))))

def threshold_centroid(cutout, **kwargs):
    """
    Returns the weighted mean centroid of values above 2/3 the maximum value in
    a cutout of an array.
    """
    x_sum = np.sum(cutout, axis=0)
    y_sum = np.sum(cutout, axis=1)
    x_fwh = np.where(x_sum >= 0.67 * max(x_sum))
    y_fwh = np.where(y_sum >= 0.67 * max(y_sum))
    x_avg = np.average(x_fwh[0], weights=x_sum[x_fwh])
    y_avg = np.average(y_fwh[0], weights=y_sum[y_fwh])
    # plt.imshow(cutout)
    # plt.scatter(x_avg, y_avg, s=2, c='red', marker='o')
    # plt.show()
    return((int(np.floor(x_avg)), int(np.floor(y_avg))))

def hybrid_centroid(cutout, **kwargs):
    # Get the maximum value of the cutout as an initial guess.
    x_max, y_max = max_value_centroid(cutout)
    # Create a smaller cutout around the initial guess.
    obj_size = 50
    reduced_cutout = cutout[x_max-obj_size:x_max+obj_size,y_max-obj_size:y_max+obj_size]
    # Get the mean weighted average of the smaller cutout.
    x_new, y_new = weighted_mean_2D(reduced_cutout, floor=True)
    # Map the centroid back to coordinates of original cutout.
    x_avg = x_max - obj_size + x_new
    y_avg = y_max - obj_size + y_new
    # plt.imshow(cutout)
    # plt.scatter(x_avg, y_avg, s=2, c='red', marker='o')
    # plt.show()
    return((x_avg, y_avg))

def align(image_stack, **kwargs):
    """
    Recieves a list of image arrays and some "cutout" range containing a common
    object to use for alignment of the image stack.

    Returns a list of image arrays of different size, aligned, and with zero
    borders where the image has been shifted.

    Args:
        image_stack (list): frames to be aligned.
        centroid (func_handle): function handle for arbitrary centroiding
        which returns a tuple.
    Returns:
        aligned_image_stack (list): new frames that have been aligned and can be
            stacked.
    """
    x, y, dx, dy = kwargs.get("cutout")
    centroid = kwargs.get("centroid")
    # Get lists of all the x and y centroids.
    x_centroids, y_centroids = [], []
    for image in image_stack:
        x_centroids.append(centroid(image[0].data[x:x+dx,y:y+dy],floor=True)[0])
        y_centroids.append(centroid(image[0].data[x:x+dx,y:y+dy],floor=True)[1])
    x_ref, y_ref = min(x_centroids), min(y_centroids)
    x_max_offset = max(x_centroids) - min(x_centroids)
    y_max_offset = max(y_centroids) - min(y_centroids)
    # Create new list of image arrays with offset.
    aligned_image_stack = []
    for image in image_stack:
        aligned_image = np.zeros((image[0].data.shape[0]+x_max_offset, image[0].data.shape[1]+y_max_offset))
        x_image_offset = centroid(image[0].data[x:x+dx,y:y+dy],floor=True)[0] - x_ref
        y_image_offset = centroid(image[0].data[x:x+dx,y:y+dy],floor=True)[1] - y_ref
        aligned_image[x_image_offset:x_image_offset+image[0].data.shape[0],y_image_offset:y_image_offset+image[0].data.shape[1]] = image[0].data
        aligned_image_stack.append(aligned_image)
    return(aligned_image_stack)

def stack(aligned_image_stack):
    """
    Receives a list of aligned images and returns their sum.

    Args:
        aligned_image_stack (list): aligned frames ready to be stacked.
    Returns:
        stacked_image (ndarray): new combined single frame.
    """
    # Check that the aligned images to be stacked have matching dimensions.
    for image in aligned_image_stack:
        if image.shape != aligned_image_stack[0].shape:
            print("Aligned image dimensions do not match!")
            break
    # If all dimensions match, initialise an empty array with those dimensions
    # into which aligned images are stacked.
    stacked_image = np.zeros(aligned_image_stack[0].shape)
    for image in aligned_image_stack:
        stacked_image += image
    return(stacked_image)

def main():
    """
    Grabs sorted lists of files from get_lists function. Creates a list of
    possible integration times from the dark files. Combines darks and flats
    into master HDUList objects. These objects are stored in dictionaries that
    have the integration times or bands as the keys.
    """
    #: path obj: Various folder locations.
    data_folder = Path("dat/")
    science_folder = Path("sci/")
    temp_folder = Path("tmp/")
    #: list of dict: Lists contain dicts with filename (str) and other keys.
    raw_dark_list, raw_flat_list, raw_target_list, raw_std_star_list = get_lists(data_folder)
    #: list of str: Contains possible integration times. No duplicates.
    possible_int_times = list(set(sub["integration_time"] for sub in raw_dark_list))
    #: list of str: Possible bands. No duplicates.
    possible_bands = list(set(sub["band"] for sub in raw_flat_list))
    #: dict of ndarray: Contains master dark objects and integration_times.
    print("Creating dark frames..."),
    master_dark_frame = {}
    for pos_int_time in possible_int_times:
        #: list of ndarray: Contains data of dark files.
        sorted_dark_list = []
        for dark in raw_dark_list:
            if dark["integration_time"] == pos_int_time:
                with fits.open(data_folder / dark["filename"]) as hdul:
                    #: ndarray: image data
                    data = hdul[0].data
                    sorted_dark_list.append(data)
        master_dark_frame[pos_int_time] = np.floor(np.median(sorted_dark_list, 0))
    print("Done!")
    #: dict of ndarray: Master flat objects, bands, and integration times.
    print("Creating flat frames..."),
    master_flat_frame = {}
    for pos_band in possible_bands:
        #: list of ndarray: Contains flat objects.
        sorted_flat_list = []
        for flat in raw_flat_list:
            if flat["band"] == pos_band:
                with fits.open(data_folder / flat["filename"]) as hdul:
                    #: ndarray: dark subtracted image data
                    data = np.subtract(hdul[0].data, master_dark_frame[flat["integration_time"]])
                    sorted_flat_list.append(data)
        master_flat_frame[pos_band] = normalise_flat(np.floor(np.median(sorted_flat_list, 0)))
    print("Done!")

    def reduce_raws(raw_list):
        """
        Reduces raw images into science images. This function loops through a
        list of raws. Each raw image is dark subtracted and then flat divided.

        Args:
            raw_list (list): Raw ndarray objects.
    , "g", "r"]    Returns:
            science_list (list): Reduced ndarray objects.
        """
        #: list of ndarray: Empty list for reduced images.
        science_list = {}
        for raw in raw_list:
            print("Reducing " + str(len(science_list)) +  " of " + str(len(raw_list)) + " images.", end="\r"),
            with fits.open(data_folder / raw["filename"]) as hdul:
                #: ndarray: Dark subtracted image data
                ds_data = np.subtract(hdul[0].data, master_dark_frame[raw["integration_time"]])
                science_list[raw["filename"]] = np.divide(ds_data, master_flat_frame[raw["band"]])
        print("\nDone!")
        return science_list

    #: dictionary of ndarray objects: Reduced target image list.
    reduced_target_list = reduce_raws(raw_target_list)
    #: dictionary of ndarray objects: Reduced standard star image list.
    reduced_std_star_list = reduce_raws(raw_std_star_list)

    print("Writing out dark frames..."),
    for key, value in master_dark_frame.items():
        write_out_fits(value, "tmp/"+"dark_"+str(key))
    print("Done!")
    print("Writing out flat frames..."),
    for key, value in master_flat_frame.items():
        write_out_fits(value, "tmp/"+"flat_"+str(key))
    print("Done!")

    counter = 0
    total = len(reduced_target_list)
    for key, value in reduced_target_list.items():
        counter += 1
        print('Writing out ' + str(counter) +  ' of ' + str(total) + ' target images.', end='\r'),
        write_out_fits(value, "sci/"+str(key))
    print("\nDone!")

    counter = 0
    total = len(reduced_std_star_list)
    for key, value in reduced_std_star_list.items():
        counter += 1
        print("Writing out " + str(counter) +  " of " + str(total) + " standard star images...", end="\r"),
        write_out_fits(value, "sci/"+str(key))
    print("\nDone!")

def test_main():

    # x, y, dx, dy = 350, 1450, 350, 350plt.imshow(cutout)
    x, y, dx, dy = 0, 0, 3351, 2531
    # the radius, in pixels, of the reference star
    obj_size = 50

    for target in ["m52"]:
        for band in ["r","g","u"]:
            unaligned_images = load_fits(path="sci/", target=target, band=band)
            aligned_images = align(unaligned_images, cutout=(x,y,dx,dy), centroid=hybrid_centroid)
            stacked_image = stack(aligned_images)
            write_out_fits(stacked_image, "sta/" + target + "_" + band + "_stacked.fits")

    # Create a list of reduced science frames for alignment.
    # Align images.
    # aligned_science_list = align_images(science_list)
    # Stack the frames and write out.
    # stacked_science_frame = stack_images(aligned_science_list)
    # write_out_fits(stacked_science_frame, filename)

    # frame size: 3352x2532px

# main()
test_main()
