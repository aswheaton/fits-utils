from fits_utils import *

# def get_errors_distance(err_r,err_g,cov_pl,cov_m52,param_pl,param_m52,g_r_m52,dist_mod):
#     err_g_r=np.sqt(err_r**2 +err_g**2)
#     dm52_dgr=(4*param_m52[0]*(g_r_m52**3) +3*param_m52[1]*(g_r_m52**2) + 2*param_m52[2]*(g_r_m52) +param_m52[3])
#     err_r_pl=np.sqt((g_r_pl**8)*cov_pl[0,0] +(g_r_pl**6)*cov_pl[1,1] +(g_r_pl**4)*cov_pl[2,2] +(g_r_pl**2)*cov_pl[3,3] +cov_pl[4,4])
#     err_r_m52=np.sqt((g_r_m52**8)*cov_m52[0,0] +(g_r_m52**6)*cov_m52[1,1] +(g_r_m52**4)*cov_m52[2,2] +(g_r_m52**2)*cov_m52[3,3] +cov_m52[4,4] +(dm52_dgr**2)*(err_g_r**2))
#     err_r_diff=np.sqt(err_r_pl**2 + err_r_m52**2)
#     err_dist_mod=err_r_diff
#     err_distance=(2**((dist_mod/5)+1)*np.log(10)*(5**(dist_mod/5)))*err_dist_mod
#     return(err_distance)
#
# def get_errors_age(err_zp_r,err_r,zpr,r_mag,err_distance,flux,distance,lum,mass):
#     dc_dr=(2**(1-(2*r_mag/5)+(2*zpr/5)))*(np.log(10))*(5**(-1-(2*r_mag/5)+(2*zpr/5)))
#     err_flux=dc_dr*np.sqt(err_zp_r**2 +err_r**2)
#     err_lum=np.sqt((4*err_flux*np.pi*distance**2)+(8*err_distance*distance*np.pi*flux)**2)
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

def absolute_magnitude(apparent_mag, apparent_mag_err_sq, distance,
                       distance_err_sq):
    """
        Returns the absolute magnitude from an apparent magnitude and known
        distance.

        Args:
            apparent_mag (flt): Apparent magnitude of the object.
            apparent_mag_err (flt): Uncertainty on the apparent magnitude.
            distance (flt): Distance in parsecs.
            distance_err_sq (flt): Uncertainty on the distance in parsecs.
        Returns:
            absolute_magnitude (flt): Absolute magnitude of the object.
            absolute_magnitude_err_sq (flt): Error on the absolute magnitude.
    """
    absolute_magnitude = apparent_mag - 5.0 * np.log10(distance / 10.0)
    absolute_magnitude_err_sq = (apparent_mag_err_sq -
                                 (5/(np.log(10)*distance))**2 * distance_err_sq)
    return(absolute_magnitude, absolute_magnitude_err_sq)

def mag_to_flux(magnitude, magnitude_err_sq):
    counts_per_sec = 10**(- magnitude / 2.5)
    counts_per_sec_err_sq = (((-1/2.5) * np.log(10) * 10**(-magnitude/2.5))**2 *
                             magnitude_err_sq)
    energy_photon = h * c / lambda_r
    flux = counts_per_sec * energy_photon
    flux_err_sq = energy_photon**2 * counts_per_sec_err_sq
    return(flux, flux_err_sq)

def distance_modulus(params1, params1_errs_sq, params2, params2_errs_sq):
    """
        Returns the average distance between two polynomial curves over a
        particular range. (Hard coded to 0.2 < x < 0.7 for now.)

        Args:
            params1: list of polynomial coefficients for a curve in ascending
                order.
            params2: list of polynomial coefficients for a second curve in
                ascending order.
        Returns:
            average_distance: float of average distance between the two curves
    """

    x_range = np.linspace(lower_colour, upper_colour, 1000)
    fit_1, fit_1_err_sq = polynomial(x_range, x_sq_err=0, params1,
                                     params1_errs_sq)
    fit_2, fit_2_err_sq = polynomial(x_range, x_sq_err=0, params2,
                                     params2_errs_sq)
    average_distance = np.mean(abs(fit_1-fit_2))
    average_distance_err_sq = np.mean(fit_1_err_sq + fit_2_err_sq)
    return(average_distance, average_distance_err_sq)

def get_age_2_electric_boogaloo(flux, flux_err_sq, distance, distance_err_sq):
    # Calculate solar luminiosty using absolute magnitudes.
    sol_flux_r, sol_flux_r_err_sq = mag_to_flux(4.65, 0.0)
    sol_lum_r = sol_flux_r * 4.0 * np.pi * (10*meters_per_parsec)**2
    sol_lum_r_err_sq = np.abs(4.0 * np.pi * (10*meters_per_parsec)**2
                              * sol_flux_r_err_sq)
    # Calculate solar luminosity using apparent magnitudes.
    # sol_flux_r = mag_to_flux(-26.93)
    # sol_lum_r = sol_flux_r * 4.0 * np.pi * AU_meters**2
    luminosity = flux * 4.0 * np.pi * distance**2
    luminosity_err_sq = np.abs((8 * np.pi * distance * flux)**2 *
                                distance_err_sq + (4 * np.pi * distance**2)**2 *
                                flux_err_sq)
    print("Cluster member luminosity r: {} +- {}".format(luminosity,
                                                         luminosity_err_sq**0.5)
         )
    age = 9E+9 * (luminosity/sol_lum_r)**(-0.6875)
    age_err = np.abs(-6.1875E+9 * (1/sol_lum_r) *
                     (luminosity/sol_lum_r)**(-1.6875) * luminosity_err_sq**0.5)
    return(age, age_err)

def main():

    # Catalogue directory.
    cat_dir = "cat/cumulative_trim/"
    # Load in the pleiades data and limit it to reasonable values for fitting.
    pleiades_data = np.loadtxt("pleiades/pleiades_johnson.txt")
    pleiades_data = correct_pleiades(pleiades_data)
    reduced_indices = np.where((pleiades_data[:,0] > lower_colour) &
                               (pleiades_data[:,0] < upper_colour))[0]
    reduced_gr_pleiades = pleiades_data[reduced_indices,0]
    reduced_r_pleiades = absolute_magnitude(pleiades_data[reduced_indices,2],
                                            pleiades_data[reduced_indices,3],
                                            dist_pl_pc, dist_pl_pc_err_sq)
    # Load in the cluster data and limit it to reasonable values for fitting.
    gr_r_catalog = np.loadtxt(cat_dir+"de_reddened_gr_r.cat")
    gr_excess, r_mag = gr_r_catalog[:,0], gr_r_catalog[:,2]
    gr_excess_err_sq, r_mag_err_sq = gr_r_catalog[:,1]**2, gr_r_catalog[:,3]**2
    reduced_indices = np.where((gr_excess > lower_colour) &
                               (gr_excess < upper_colour))[0]
    reduced_gr_excess = gr_excess[reduced_indices]
    reduced_gr_excess_err_sq = reduced_gr_exc
    reduced_r_mag = r_mag[reduced_indices]

    # Fit a fourth order polynomial to both datasets.
    params_pleiades, cov_pleiades = np.polyfit(reduced_gr_pleiades,
                                               reduced_r_pleiades, deg=4,
                                               cov=True)
    params_pleiades = np.flip(params_pleiades)
    cov_pleiades_sq = cov_pleiades**2
    params, cov = np.polyfit(reduced_gr_excess, reduced_r_mag, deg=4, cov=True)
    params = np.flip(params)
    cov_sq = cov**2

    # Calculate the distance modulus between the two fits and
    mean_distance, mean_distance_err_sq = distance_modulus(params_pleiades,
                                                           cov_pleiades_sq,
                                                           params, cov_sq)
    distance_parsecs = 10.0**((mean_distance / 5) + 1)
    distance_parsecs_err_sq = ((np.log(10) / 5 * 10**(mean_distance / 5))**2 *
                               mean_distance_err_sq)
    distance_meters = meters_per_parsec * distance_parsecs
    distance_meters_err_sq = distance_parsecs_err_sq * meters_per_parsec

    # print("Distance Modulus: {}".format(mean_distance))
    # print("Max Member Apparent Magnitude: {}".format(np.min(r_mag)))
    # print("Corresponding flux: {}".format(mag_to_flux(np.min(r_mag))))
    # print("Max Member Absolute Magnitude: {}".format(absolute_magnitude(np.min(r_mag),distance_parsecs)))

    brightest_index = np.where(r_mag == np.min(r_mag))[0]
    max_r_flux, max_r_flux_err_sq = mag_to_flux(r_mag[brightest_index],
                                                r_mag_err_sq[brightest_index])
    cluster_age, cluster_age_err = get_age(max_r_flux, max_r_flux_err_sq,
                                           distance_meters,
                                           distance_meters_err_sq)

    # distances = np.linspace(8.0,20.0,1000)
    # plt.plot(distances, new_get_age(max_r_flux, distances)/1000000, 'b-')
    # plt.plot(distances, legacy_get_age(max_r_flux, distances), 'r-')
    # plt.show()

    print("Distance to the cluster: {} pc.".format(distance_parsecs))
    print("Age of the cluster: {} myrs.".format(cluster_age/1000000))

    #err_distance=get_errors_distance(err_r,err_g,cov_pl,cov_m52,param_pl,param_m52,cor_g_r_m52,dist_mod)
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
