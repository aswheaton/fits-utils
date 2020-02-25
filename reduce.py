"""Script for reducing images.

"""
from fits_utils import *

def main():
    """
    Grabs sorted lists of files from get_lists function. Creates a list of
    possible integration times from the dark files. Combines darks and flats
    into master HDUList objects. These objects are stored in dictionaries that
    have the integration times or bands as the keys.
    """
    # Generate the config file.
    # gen_config()
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

    #: dictionary of ndarray objects: Reduced target image list.
    reduced_target_list = reduce_raws(raw_target_list, master_dark_frame, master_flat_frame, data_folder)
    #: dictionary of ndarray objects: Reduced standard star image list.
    reduced_std_star_list = reduce_raws(raw_std_star_list, master_dark_frame, master_flat_frame, data_folder)

    print("Writing out dark frames..."),
    for key, value in master_dark_frame.items():
        write_out_fits(value, "tmp/dark_{}".format(key))
    print("Done!")
    print("Writing out flat frames..."),
    for key, value in master_flat_frame.items():
        write_out_fits(value, "tmp/flat_{}".format(key))
    print("Done!")

    counter = 0
    total = len(reduced_target_list)
    for key, value in reduced_target_list.items():
        counter += 1
        print('Writing out {} of {} target images.'.format(counter, total), end='\r'),
        write_out_fits(value, "sci/{}".format(key))
    print("\nDone!")

    counter = 0
    total = len(reduced_std_star_list)
    for key, value in reduced_std_star_list.items():
        counter += 1
        print("Writing out {} of {} standard star images...".format(counter, total), end="\r"),
        write_out_fits(value, "sci/{}".format(key))
    print("\nDone!")

if __name__ == '__main__':
    main()
