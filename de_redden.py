from fits_utils import *
import matplotlib.pyplot as plt

def main():

    # Central wavelength of filters, in micrometers.
    r_lambda, g_lambda, u_lambda = 0.6231, 0.4770, 0.3543
    # Zero points from zero-point-calculator.
    zpr, zpg, zpu = get_zero_points(1.04) # Rory's ZP: 30.236, 29.719, 27.075
    print("zpr = {}, zpg = {}, zpu = {}".format(str(zpr)[:5], str(zpg)[:5], str(zpu)[:5]))
    # Get the zero point corrected catalogue and error.
    r_mag, r_err, g_mag, g_err, u_mag, u_err = load_cat("cat/ugr.cat", zpr, zpg, zpu)
    # error = (g_err**2 + r_err**2 + u_err**2)**0.5

    pleiades_data = np.loadtxt("pleiades/pleiades_johnson.txt")
    pleiades_data = correct_pleiades(pleiades_data)
    indices = np.where(pleiades_data[:,0] < 0.5)[0]
    reduced_data = pleiades_data[indices,:]
    # pleiades_coeffs = np.flip(np.polyfit(pleiades_data[:,0], pleiades_data[:,1], 4),axis=0)
    pleiades_coeffs = np.flip(np.polyfit(reduced_data[:,0], reduced_data[:,1], 4),axis=0)
    # Calculate the colour excess.
    gr_excess = g_mag - r_mag # x-axis variable
    ug_excess = u_mag - g_mag # y-axis variable
    # Calculate error on colour excess.
    gr_excess_err = g_err + r_err
    ug_excess_err = u_err + g_err

    # Calculate the slope of the reddening vector according to Cardelli et al.
    cardelli_consts = {"r" : cardelli_const(r_lambda),
                       "g" : cardelli_const(g_lambda),
                       "u" : cardelli_const(u_lambda)
                      }
    cardelli_slope = get_cardelli_slope(cardelli_consts)
    # 2d Array containing the various reddening vector magnitudes, x and y
    # components of reddening vector and the chi-squared value associated with
    # the particular reddening vector magnitude.
    params_and_fit = []

    mp_x = (1/gr_excess.shape[0])*np.sum(gr_excess)
    mp_y = (1/ug_excess.shape[0])*np.sum(ug_excess)
    y_cept = mp_y - cardelli_slope*mp_x

    # Iterate over reasonable range of values for the reddening vector magnitude.
    for red_vec_mag in np.linspace(0.6, 1.75, 1000):
        # Separate reddening vector into components in colour-colour space.
        red_vec_x = (red_vec_mag**2 / (1 + cardelli_slope**2))**0.5
        red_vec_y = (red_vec_mag**2 / (1 + cardelli_slope**-2))**0.5
        # Remove outliers which might skew the chi-squared minimisation.
        reduced_indices = remove_outliers(gr_excess, 0.2, 0.7)
        gr_excess_shifted = gr_excess[reduced_indices] - red_vec_x
        ug_excess_shifted = ug_excess[reduced_indices] - red_vec_y
        # Get the minimum chi-squared value for the particular magnitude.
        chi_squ = get_chi_squ(gr_excess_shifted, ug_excess_shifted, polynomial,
                                 pleiades_coeffs, ug_excess_err[reduced_indices]
                                )
        params_and_fit.append([red_vec_mag, red_vec_x, red_vec_y, chi_squ])

    # Dirty data type change so minimisation can be performed.
    params_and_fit = np.array(params_and_fit)
    row = minimiser(params_and_fit[:,3])
    best_red_vec_mag, best_red_vec_x, best_red_vec_y, chi_squ_min = params_and_fit[row,:][0]

    g_abs = best_red_vec_y * ((g_lambda/u_lambda) - 1)**-1
    u_abs = best_red_vec_y + g_abs
    r_abs = g_abs - best_red_vec_x

    v_abs_g = g_abs / cardelli_consts["g"]
    v_abs_u = u_abs / cardelli_consts["u"]
    v_abs_r = r_abs / cardelli_consts["r"]

    print("Cardelli slope: {}\nMagnitide: {}\nx-comp: {}\ny-comp: {}\nChi-Squ: {}".format(cardelli_slope, best_red_vec_mag, best_red_vec_x, best_red_vec_y, chi_squ_min))
    print("A_g = {}\nA_u = {}\nA_r = {}".format(g_abs,u_abs,r_abs))
    print("A_v = {} (from A_g)\nA_v = {} (from A_u)\nA_v = {} (from A_r)".format(v_abs_g,v_abs_u,v_abs_r))

    de_reddened_gr_excess = gr_excess - best_red_vec_x
    de_reddened_ug_excess = ug_excess - best_red_vec_y

    new_catalogue = np.column_stack((u_mag-u_abs, g_mag-g_abs, r_mag-r_abs))
    np.savetxt("cat/de_reddened_ugr.cat", new_catalogue)

    # Plot chi-sqaured as a function of reddening vector magnitude and the
    # reddened data alongside the de-reddened data and the Pleiades data.
    plt.plot(params_and_fit[:,0], params_and_fit[:,3])
    dict = {}
    dict["M52 Uncorrected"] = (gr_excess,ug_excess,"o")
    dict["M52 De-reddened"] = (de_reddened_gr_excess,de_reddened_ug_excess,"o")
    value_range = np.linspace(-0.8, 1.0, 1000)
    dict["Pleiades Data"] = (pleiades_data[:,0], pleiades_data[:,1], "o")
    dict["Pleiades Fit"] = (value_range, polynomial(value_range, pleiades_coeffs), "-")
    dict["Cardelli Slope"] = (value_range, cardelli_slope*value_range + y_cept, "-")
    plot_diagram(dict, x_label="Colour:(g-r)", y_label="Colour:(u-g)",
                 sup_title="M52\nColour-Colour Diagram",
                 legend=True, filename="M52_Colour-Colour_Diagram"
                )

    # Load in the larger g and r catalogue of objects which are invisible in u.
    catalog = np.loadtxt("cat/gr.cat")
    r_mag, r_err = get_mag(catalog[:,5], catalog[:,6], zpr)
    g_mag, g_err = get_mag(catalog[:,3], catalog[:,4], zpg)
    # De-redden the larger catalogue with newly found r_abs and g_abs values.
    de_reddened_r_mag = r_mag - r_abs
    de_reddened_g_mag = g_mag - g_abs
    # Calculate the de-reddened colour excess.
    de_reddened_gr_excess = de_reddened_g_mag - de_reddened_r_mag
    # Write the corrected catalogue out.
    de_reddened_gr_r = np.column_stack((de_reddened_gr_excess, de_reddened_r_mag))
    np.savetxt("cat/de_reddened_gr_r.cat", de_reddened_gr_r)
    # Plot the de-reddened diagram.
    dict = {"M52 g vs. g-r":(de_reddened_gr_excess,de_reddened_g_mag,'o'),
            "Pleiades g vs. g-r":(pleiades_data[:,0], pleiades_data[:,0]+pleiades_data[:,2], 'o')
           }
    plot_diagram(dict, x_label="Colour:(g-r)", y_label="Magnitude: g",
                 sup_title="M52\nColour-Colour Diagram",
                 legend=True, filename="M52_Colour-Colour_Diagram"
                )

main()
