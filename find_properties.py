import numpy as np
import matplotlib.pyplot as plt
from fits_utils import *
from scipy import optimize

solar_lum= 3.828E+26
solar_mass= 1.99E+30
h=6.63E-34
c=3E+8
cen_wav_r= 658E-9
dist_pl_pc=150


def absolute_magnitude(apparent_mag, distance):
    absolute_magnitude = apparent_mag - 5.0 * np.log10(distance / 10.0)
    return(absolute_magnitude)

def mag_to_flux(m,zp):
    counts= 10**((zp-m)/2.5)
    energy_photon=h*c/cen_wav_r
    flux=counts*energy_photon
    return(flux)

def remove_outliers(g_r,r):
    index= np.where((g_r<-0.2) & (g_r>-0.7))[0]
    return(g_r[index],r[index])

def get_distance_2(g_r,r,param_pl):
    #USING BIN METHOD
    cut_off=np.mean(g_r)
    index_bin_1= np.where((g_r>np.amin(g_r))&(g_r<cut_off))[0]
    index_bin_2= np.where((g_r>cut_off)&(g_r<np.amax(g_r)))[0]
    bin_1_g_r= g_r[index_bin_1]
    bin_2_g_r= g_r[index_bin_2]
    bin_1_r=r[index_bin_1]
    bin_2_r=r[index_bin_2]
    mean_bin_1_r=np.mean(bin_1_r)
    mean_bin_2_r=np.mean(bin_2_r)
    mean_bin_1_g_r=np.mean(bin_1_g_r) #NOT CENTRAL VALUE OF BIN
    mean_bin_2_g_r=np.mean(bin_2_g_r) #NOT CENTRAL VALUE OF BIN
    bin_1_r_pl=polynomial(mean_bin_1_g_r,*param_pl)
    bin_2_r_pl=polynomial(mean_bin_2_g_r,*param_pl)
    r_bin_1_diff=abs(bin_1_r_pl-mean_bin_1_r)
    r_bin_2_diff=abs(bin_2_r_pl- mean_bin_2_r)
    dist_mod=(r_bin_1_diff+r_bin_2_diff)/2
    print(dist_mod)
    exponent=((dist_mod/5)+1)
    distance= 10**(exponent)
    #err_dist_mod=
    #err_distance=(2**((dist_mod/5)+1)*np.log(10)*(5**(dist_mod/5)))*err_dist_mod
    return(distance)


def distance_modulus(params1, params2):
    """
        Returns the average distance between two polynomial curves over a
        particular range. (Hard coded to 0.2 < x < 0.7 for now.)

        Args:
            params1: list of polynomial coefficients for a curve in ascending order
            params2: list of polynomial coefficients for a second curve in ascending order
        Returns:
            average_distance: float of average distance between the two curves
    """

    x_range = np.linspace(0.2, 0.7, 1000)
    average_distance = np.sum(abs(polynomial(x_range,params1)-polynomial(x_range, params2))) / 1000
    return(average_distance)

def get_age(flux, distance):
    """
        James write a docstring please <3
        Method to determine the age of a source based on the flux. Flux is
        converted to a luminosity using the equation L = 4 * pi * (d**2) * F

    """
    # Get the age of a star from the mass luminosity relationship.
    # Constant formula: (7.875e-04 * c**2 * solar mass) / (solar luminosity)**(1/3)
    luminosity = flux * 4.0 * np.pi * distance**2
    age = 1.941499183e35 * luminosity**(-2/3)
    return(age)

def get_errors_distance(err_r,err_g,cov_pl,cov_m52,param_pl,param_m52,g_r_m52,dist_mod):
    err_g_r=np.sqrt(err_r**2 +err_g**2)
    dm52_dgr=(4*param_m52[0]*(g_r_m52**3) +3*param_m52[1]*(g_r_m52**2) + 2*param_m52[2]*(g_r_m52) +param_m52[3])
    err_r_pl=np.sqrt((g_r_pl**8)*cov_pl[0,0] +(g_r_pl**6)*cov_pl[1,1] +(g_r_pl**4)*cov_pl[2,2] +(g_r_pl**2)*cov_pl[3,3] +cov_pl[4,4])
    err_r_m52=np.sqrt((g_r_m52**8)*cov_m52[0,0] +(g_r_m52**6)*cov_m52[1,1] +(g_r_m52**4)*cov_m52[2,2] +(g_r_m52**2)*cov_m52[3,3] +cov_m52[4,4] +(dm52_dgr**2)*(err_g_r**2))
    err_r_diff=np.sqrt(err_r_pl**2 + err_r_m52**2)
    err_dist_mod=err_r_diff
    err_distance=(2**((dist_mod/5)+1)*np.log(10)*(5**(dist_mod/5)))*err_dist_mod
    return(err_distance)

def get_errors_age(err_zp_r,err_r,zpr,r_mag,err_distance,flux,distance,lum,mass):
    dc_dr=(2**(1-(2*r_mag/5)+(2*zpr/5)))*(np.log(10))*(5**(-1-(2*r_mag/5)+(2*zpr/5)))
    err_flux=dc_dr*np.sqrt(err_zp_r**2 +err_r**2)
    err_lum=np.sqrt((4*err_flux*np.pi*distance**2)+(8*err_distance*distance*np.pi*flux)**2)
    err_mass=err_lum*((solar_lum**-(1/3))*(1/(3*lum**(2/3))))
    err_age= (10**10)*-2.5*(mass**-3.5)*err_mass
    return(err_age)

def main():

    # Catalogue directory.
    cat_dir = "cat/cumulative_trim/"
    # Central wavelength of filters, in micrometers.
    r_lambda, g_lambda, u_lambda = 0.6231, 0.4770, 0.3543
    # Zero points from zero-point-calculator.
    zpr, zpg, zpu = get_zero_points(1.04)
    # Load in the pleiades data and limit it to reasonable values for fitting.
    pleiades_data = np.loadtxt("pleiades/pleiades_johnson.txt")
    pleiades_data = correct_pleiades(pleiades_data)
    reduced_indices = np.where((pleiades_data[:,0] > 0.2) & (pleiades_data[:,0] < 0.7))[0]
    reduced_gr_pleiades = pleiades_data[reduced_indices,0]
    reduced_r_pleiades = absolute_magnitude(pleiades_data[reduced_indices,2],dist_pl_pc)
    # Load in the cluster data and limit it to reasonable values for fitting.
    gr_r_catalog = np.loadtxt(cat_dir+"de_reddened_gr_r.cat")
    gr_excess, r_mag = gr_r_catalog[:,0], gr_r_catalog[:,1]
    reduced_gr_excess,reduced_r_mag= remove_outliers(gr_excess,r_mag)

    # param_pleiades, cov_pl = get_fit(polynomial, reduced_gr_pleiades, reduced_r_pleiades)
    # params, covariance = get_fit(polynomial, reduced_gr_excess, reduced_r_mag)

    # Fit a fourth order polynomial to both datasets.
    params_pleiades, cov_pleiades = np.polyfit(reduced_gr_pleiades, reduced_r_pleiades, deg=4, cov=True)
    params_pleiades = np.flip(params_pleiades)
    params, cov = np.polyfit(reduced_gr_excess, reduced_r_mag, deg=4, cov=True)
    params = np.flip(params)

    # Calculate the distance modulus between the two fits and
    mean_distance = distance_modulus(params_pleiades, params)
    distance_parsecs = 10.0**((mean_distance / 5) + 1)
    distance_meters= 3.08567782e16 * distance_parsecs

    min_r_flux = mag_to_flux(np.min(reduced_r_mag), zpr)
    cluster_age = get_age(min_r_flux, distance_meters)

    print("Distance to the cluster: {} pc.".format(distance_parsecs))
    print("Age of the cluster: {} mya.".format(cluster_age / 1000000.0))

    #err_distance=get_errors_distance(err_r,err_g,cov_pl,cov_m52,param_pl,param_m52,reduced_gr_excess,dist_mod)
    #err_age=get_errors_age(err_zp_r,err_r,zpr,r_min,err_distance,r_min_flux,distance,r_min_lum,r_min_mass)
    #NOTE: This code doesnt currently possess a err_zp_r, err_r, err_g need to be added from updated catalogue & from zp calc

    x_range=np.linspace(np.amin(reduced_gr_excess),np.amax(reduced_gr_excess), 1000)

    dict = {"Pleiades"     : (reduced_gr_pleiades, reduced_r_pleiades, 'o'),
            "M52"          : (reduced_gr_excess, reduced_r_mag, 'o'),
            "Pleiades Fit" : (x_range, polynomial(x_range, params_pleiades), '-'),
            "M52 Fit"      : (x_range, polynomial(x_range, params), '-')
            }
    plot_diagram(dict, xlabel="G-R Colour", ylabel="R Magnitude")

main()
