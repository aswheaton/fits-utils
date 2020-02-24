"""Script for removing intervening stars from a star catalog.

Want to remove intervening stars from a catalog of stars so that only members
of the star cluster remain.

Theory:
    Stars in the cluster should be of similar age and therefore should have
    similar brightness when compared with intervening stars. Intervening stars
    should be much brighter than stars in the cluster.

Method:
    Read in the catalog of stars as a numpy array and plot a distribution of
    star brightness. Majority of stars should be in the cluster, so remove any
    stars which are outside of x standard deviations from the distribution.

    Create a new catalog of stars with outliers removed. This should only
    contain cluster members.

Catalog:
    [0]: NUMBER : Runing object number.
    [1]: FLUX_APER : Flux vector within circular aperture.
    [2]: FLUXERR_APER : RMS error for aperture.
    [3]: X_IMAGE : X-pos.
    [4]: Y_IMAGE : Y-pos.

"""
import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy import stats

def import_catalog(filename):
    """Function for reading in catalog of objects.

    Reads in catalog of objects as an ndarray.

    Args:
        filename (str): Name of the catalog.

    Returns:
        catalog (ndarray): Ndarray of objects. Contains various parameters, only
            care about object number and brightness.

    """
    #: ndarray: Catalog of stars. Contains fluxes etc.
    catalog = np.loadtxt(filename)
    return catalog

def plot_dist(flux, mu, sigma):
    """Plots a distribution of fluxes.

    Args:
        flux (ndarray): Flux values.
        mu (flt): Mean flux value.
        sigma (flt): Standard deviation of the fluxes.

    """
    #: ndarrays: Histogram and bin edges from data.
    sorted_flux = np.sort(flux)
    fit = stats.norm.pdf(sorted_flux, mu, sigma)
    plt.plot(sorted_flux, fit, '-o')
    xlines = [mu, mu + sigma, mu - sigma]
    for line in xlines:
        plt.axvline(x=line)
    # plt.hist(flux)
    plt.xlabel('Flux')
    plt.ylabel('Frequency')
    plt.title('Distribution of Catalog Fluxes')
    plt.show()

def remove_outliers(catalog, k):
    """Removes outliers from catalog.

    Finds the standard deviation of the object fluxes. Recieves a user defined
    standard deviation limit and removes any objects that are outside of this
    limit.

    Args:
        catalog (ndarray): Object catalog. Contains IDs and fluxes as well as
            other uninteresting values.
        k (int): Standard deviation limit for the fluxes. Data that is outside
            this range is not included in the new catalog.

    Returns:
        new_catalog (ndarray): Catalog with outliers removed.

    """
    #: flt: Mean and standard deviation of the data.
    mu, sigma = np.mean(catalog[:,1]), np.std(catalog[:,1])
    plot_dist(catalog[:,1], mu, sigma)
    #: ndarray: New array which excludes outliers.
    new_catalog = catalog[np.logical_and(
                                         catalog[:,1] <= mu + k*sigma,
                                         catalog[:,1] >= mu - k*sigma
                                        )]
    return new_catalog

def make_vals(catalog):
    flux = catalog[:,1]
    u = zp['u'] - 2.5 * np.log(flux)/np.log(10)
    pass

def main():
    """Main method.

    """
    #: str: User defined filename from command line input.
    filename = sys.argv[1]
    #: ndarray: Imported catalog of stars.
    catalog = import_catalog('catalogs/{}.cat'.format(filename))
    #: ndarray: New catalog of stars, with outliers removed.
    new_catalog = remove_outliers(remove_outliers(catalog, 1), 2)


if __name__ == '__main__':
    main()
