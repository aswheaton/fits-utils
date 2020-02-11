"""
Script for combining the three individual catalogs of g, r, and u band stars
into one catalog.

"""

import numpy as np

def read_catalog(filename):
    """Function for reading in the data from a catalog.

    Reads in .cat files as ndarrays. Rearranges the ndarrays so that the columns
    are in a more suitable order.

    Args:
        filename (str): Filename and location of the specific catalog.

    Returns:
        catalog (ndarray): Ndarray of the catalog data.

    """
    catalog = np.loadtxt(filename)
    permutation = [0, 3, 4, 1, 2]
    index = np.empty_like(permutation)
    index[permutation] = np.arange(len(permutation))
    catalog[:] = catalog[:, index]
    return catalog

def create_catalog(catalog, filename):
    """Method for creating a new catalog file from an ndarray.

    Takes an ndarray and writes it out to a new catalog file with a '.cat'
    extension. Creates a header for this file with the indexes included.

    Args:
        catalog (ndarray): Merged catalog data.
        filename (str): Filename of the new catalog. Do not include the
            extension, this is added automatically.

    """
    #: str: New header text for the output file.
    header_txt = '\n'.join(['[0] : NUMBER',
                            '[1] : ALPHAPEAK_J2000',
                            '[2] : DELTAPEAK_J2000',
                            '[3] : FLUX_APER_G',
                            '[4] : FLUXERR_APER_G',
                            '[5] : FLUX_APER_R',
                            '[6] : FLUXERR_APER_R',
                            '[7] : FLUX_APER_U',
                            '[8] : FLUXERR_APER_U'])
    np.savetxt('catalogs/{}.cat'.format(filename), catalog, header=header_txt)

def match_sources(catalog_1, catalog_2):
    """Function for matching sources in two catalogs.

    Takes two catalogs. For each catalog, matches sources based on their RA and
    DEC. Returns a new catalog.

    Args:
        catalog_1 (ndarray): First catalog to be matched.
        catalog_2 (ndarray): Second catalog to be matched.

    Returns:
        new_catalog (ndarray): New catalog containing the RA, DEC, fluxes, and
            flux errors of the two older catalogs.

    """
    new_catalog = []
    for source in catalog_1:
        index = np.where((np.isclose(catalog_2[:,1],source[1])) &
                         (np.isclose(catalog_2[:,2],source[2])))[0]
        if np.size(index) != 0:
            new_catalog.append(np.append(source, catalog_2[index[0], 3:5]))
    new_catalog = np.vstack(new_catalog)
    return new_catalog

def main():
    #: tuple: Catalog bands. Should just be g, r, u.
    bands = ('g', 'r', 'u')
    #: dict: For storing the raw catalog ndarrays.
    catalog = {}
    for band in bands:
        #: Store the catalogs.
        catalog[band] = read_catalog('catalogs/{}.cat'.format(band))
    #: ndarray: New merged catalog for all bands.
    new_catalog = match_sources(match_sources(catalog['g'], catalog['r']),
                                catalog['u'])
    create_catalog(new_catalog, 'combined')

if __name__ == '__main__':
    main()
