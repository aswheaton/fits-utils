import numpy as np
import matplotlib.pyplot as plt
from fits_utils import *
from scipy import optimize

# def remove_outliers(g_r,r):
#     index= np.where((g_r<-0.2) & (g_r>-0.7))[0]
#     return(g_r[index],r[index])
#
# def get_distance_2(g_r,r,param_pl):
#     #USING BIN METHOD
#     cut_off=np.mean(g_r)
#     index_bin_1= np.where((g_r>np.amin(g_r))&(g_r<cut_off))[0]
#     index_bin_2= np.where((g_r>cut_off)&(g_r<np.amax(g_r)))[0]
#     bin_1_g_r= g_r[index_bin_1]
#     bin_2_g_r= g_r[index_bin_2]
#     bin_1_r=r[index_bin_1]
#     bin_2_r=r[index_bin_2]
#     mean_bin_1_r=np.mean(bin_1_r)
#     mean_bin_2_r=np.mean(bin_2_r)
#     mean_bin_1_g_r=np.mean(bin_1_g_r) #NOT CENTRAL VALUE OF BIN
#     mean_bin_2_g_r=np.mean(bin_2_g_r) #NOT CENTRAL VALUE OF BIN
#     bin_1_r_pl=polynomial(mean_bin_1_g_r,*param_pl)
#     bin_2_r_pl=polynomial(mean_bin_2_g_r,*param_pl)
#     r_bin_1_diff=abs(bin_1_r_pl-mean_bin_1_r)
#     r_bin_2_diff=abs(bin_2_r_pl- mean_bin_2_r)
#     dist_mod=(r_bin_1_diff+r_bin_2_diff)/2
#     print(dist_mod)
#     exponent=((dist_mod/5)+1)
#     distance= 10**(exponent)
#     #err_dist_mod=
#     #err_distance=(2**((dist_mod/5)+1)*np.log(10)*(5**(dist_mod/5)))*err_dist_mod
#     return(distance)

# def get_errors_distance(err_r,err_g,cov_pl,cov_m52,param_pl,param_m52,g_r_m52,dist_mod):
#     err_g_r=np.sqrt(err_r**2 +err_g**2)
#     dm52_dgr=(4*param_m52[0]*(g_r_m52**3) +3*param_m52[1]*(g_r_m52**2) + 2*param_m52[2]*(g_r_m52) +param_m52[3])
#     err_r_pl=np.sqrt((g_r_pl**8)*cov_pl[0,0] +(g_r_pl**6)*cov_pl[1,1] +(g_r_pl**4)*cov_pl[2,2] +(g_r_pl**2)*cov_pl[3,3] +cov_pl[4,4])
#     err_r_m52=np.sqrt((g_r_m52**8)*cov_m52[0,0] +(g_r_m52**6)*cov_m52[1,1] +(g_r_m52**4)*cov_m52[2,2] +(g_r_m52**2)*cov_m52[3,3] +cov_m52[4,4] +(dm52_dgr**2)*(err_g_r**2))
#     err_r_diff=np.sqrt(err_r_pl**2 + err_r_m52**2)
#     err_dist_mod=err_r_diff
#     err_distance=(2**((dist_mod/5)+1)*np.log(10)*(5**(dist_mod/5)))*err_dist_mod
#     return(err_distance)
#
# def get_errors_age(err_zp_r,err_r,zpr,r_mag,err_distance,flux,distance,lum,mass):
#     dc_dr=(2**(1-(2*r_mag/5)+(2*zpr/5)))*(np.log(10))*(5**(-1-(2*r_mag/5)+(2*zpr/5)))
#     err_flux=dc_dr*np.sqrt(err_zp_r**2 +err_r**2)
#     err_lum=np.sqrt((4*err_flux*np.pi*distance**2)+(8*err_distance*distance*np.pi*flux)**2)
#     err_mass=err_lum*((solar_lum**-(1/3))*(1/(3*lum**(2/3))))
#     err_age= (10**10)*-2.5*(mass**-3.5)*err_mass
#     return(err_age)

solar_lum = 3.828E+26
solar_mass = 1.989E+30
h = 6.62607015E-34
c = 3E+8
lambda_r = 6.231E-7
dist_pl_pc = 150
AU_meters = 1.496E+11
meters_per_parsec = 3.08567782E+16
lower_colour, upper_colour = -0.6, -0.2

def absolute_magnitude(apparent_mag, distance):
    absolute_magnitude = apparent_mag - 5.0 * np.log10(distance / 10.0)
    return(absolute_magnitude)

def mag_to_flux(magnitude):
    counts_per_sec = 10**(- magnitude / 2.5)
    energy_photon = h * c / lambda_r
    flux = counts_per_sec * energy_photon
    return(flux)

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

    x_range = np.linspace(lower_colour, upper_colour, 1000)
    average_distance = np.mean(abs(polynomial(x_range,params1)-polynomial(x_range, params2)))
    return(average_distance)

def legacy_get_age(flux,distance):
    r_min_lum=flux*(4*np.pi*(distance**2))
    r_min_mass= (r_min_lum/solar_lum)**(1/3)
    age= (r_min_mass**-2.5)* (10**10)
    return(age)

def new_get_age(flux, distance):
    """
        James write a docstring please <3
        Method to determine the age of a source based on the flux. Flux is
        converted to a luminosity using the equation L = 4 * pi * (d**2) * F
    """
    # Get the age of a star from the mass luminosity relationship.
    # Constant formula: (7.875e-04 * c**2 * solar mass) / (solar luminosity)**(1/3)
    luminosity = flux * 4.0 * np.pi * distance**2
    age = 1E+10 / (luminosity / solar_lum)**(3/4)
    return(age)

def get_age_2_electric_boogaloo(flux, distance):
    sol_flux_r = mag_to_flux(4.65)
    sol_lum_r = sol_flux_r * 4.0 * np.pi * (10*meters_per_parsec)**2
    # sol_flux_r = mag_to_flux(-26.93)
    # sol_lum_r = sol_flux_r * 4.0 * np.pi * AU_meters**2
    luminosity = flux * 4.0 * np.pi * distance**2
    print("Solar luminosity r: ", sol_lum_r)
    print("Cluster member luminosity r: ", luminosity)
    age = 9E+9 * (luminosity/sol_lum_r)**(-0.6875)
    return(age)

def main():

    # Catalogue directory.
    cat_dir = "cat/cumulative_trim/"
    # Load in the pleiades data and limit it to reasonable values for fitting.
    pleiades_data = np.loadtxt("pleiades/pleiades_johnson.txt")
    pleiades_data = correct_pleiades(pleiades_data)
    reduced_indices = np.where((pleiades_data[:,0] > lower_colour) & (pleiades_data[:,0] < upper_colour))[0]
    reduced_gr_pleiades = pleiades_data[reduced_indices,0]
    reduced_r_pleiades = absolute_magnitude(pleiades_data[reduced_indices,2],dist_pl_pc)
    # Load in the cluster data and limit it to reasonable values for fitting.
    gr_r_catalog = np.loadtxt(cat_dir+"de_reddened_gr_r.cat")
    gr_excess, r_mag = gr_r_catalog[:,0], gr_r_catalog[:,1]
    reduced_indices = np.where((gr_excess > lower_colour) & (gr_excess < upper_colour))[0]
    reduced_gr_excess = gr_excess[reduced_indices]
    reduced_r_mag = r_mag[reduced_indices]

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
    distance_meters= meters_per_parsec * distance_parsecs

    print("Distance Modulus: {}".format(mean_distance))
    print("Max Member Apparent Magnitude: {}".format(np.min(r_mag)))
    print("Corresponding flux: {}".format(mag_to_flux(np.min(r_mag))))
    print("Max Member Absolute Magnitude: {}".format(absolute_magnitude(np.min(r_mag),distance_parsecs)))

    max_r_flux = mag_to_flux(np.min(r_mag))
    cluster_age = get_age_2_electric_boogaloo(max_r_flux, distance_meters)

    # distances = np.linspace(8.0,20.0,1000)
    # plt.plot(distances, new_get_age(max_r_flux, distances)/1000000, 'b-')
    # plt.plot(distances, legacy_get_age(max_r_flux, distances), 'r-')
    # plt.show()

    print("Distance to the cluster: {} pc.".format(distance_parsecs))
    print("Age of the cluster: {} myrs.".format(cluster_age/1000000))

    #err_distance=get_errors_distance(err_r,err_g,cov_pl,cov_m52,param_pl,param_m52,reduced_gr_excess,dist_mod)
    #err_age=get_errors_age(err_zp_r,err_r,zpr,r_min,err_distance,r_min_flux,distance,r_min_lum,r_min_mass)
    #NOTE: This code doesnt currently possess a err_zp_r, err_r, err_g need to be added from updated catalogue & from zp calc

    x_range=np.linspace(lower_colour,upper_colour, 1000)

    dict = {"Pleiades"     : (reduced_gr_pleiades, reduced_r_pleiades, 'o'),
            "M52"          : (reduced_gr_excess, reduced_r_mag, 'o'),
            "Pleiades Fit" : (x_range, polynomial(x_range, params_pleiades), '-'),
            "M52 Fit"      : (x_range, polynomial(x_range, params), '-')
            }
    plot_diagram(dict, x_label="G-R Colour", y_label="R Magnitude", legend=True)

main()
