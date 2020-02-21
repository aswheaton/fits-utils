from fits_utils import *

def main():

    # Central wavelength of filters, in Angstroms.
    r_lambda, g_lambda, u_lambda = 6231, 4770, 3543
    # Zero points from zero-point-calculator.
    zpr, zpg, zpu = 30.236, 29.719, 27.075
    # Get the zero point corrected catalogue and error.
    r_mag, r_err, g_mag, g_err, u_mag, u_err = load_cat("cat/ugr.cat", zpr, zpg, zpu)
    error = (g_err**2 + r_err**2 + u_err**2)**0.5
    # Pleaides fit values in E(u-g) vs. E(g-r) colour-colour space.
    A, B, C, D, E = 1.283, -0.107, -2.748, 8.477, -3.141
    pleiades_coeffs = [A, B, C, D, E]

    # Calculate the colour excess.
    gr_excess = g_mag - r_mag # x-axis variable
    ug_excess = u_mag - g_mag # y-axis variable
    # Calculate error on colour excess.
    gr_excess_err = g_err + r_err
    ug_excess_err = u_err + g_err

    # Calculate the slope of the reddening vector according to Cardelli et al.
    cardelli_consts = {'r' : cardelli_const(r_lambda),
                       'g' : cardelli_const(g_lambda),
                       'u' : cardelli_const(u_lambda)
                      }
    cardelli_slope = get_cardelli_slope(cardelli_consts)
    # 2d Array containing the various reddening vector magnitudes, x and y
    # components of reddening vector and the chi-squared value associated with
    # the particular reddening vector magnitude.
    params_and_fit = []
    # Iterate over reasonable range of values for the reddening vector magnitude.
    for red_vec_mag in np.linspace(-2.0,  2.0, 1000):
        red_vec_x = red_vec_mag**2 / (1 + red_vec_mag**2)
        red_vec_y = red_vec_mag**2 / (1 + red_vec_mag**-2)
        gr_excess_shifted = gr_excess - red_vec_x
        ug_excess_shifted = ug_excess - red_vec_y
        chi_squ = get_chi_squ(gr_excess, ug_excess, gr_excess_shifted,
                              ug_excess_shifted, polynomial, pleiades_coeffs,
                              error
                              )
        params_and_fit.append([red_vec_mag, red_vec_x, red_vec_y, chi_squ])

    # Dirty data type change so minimisation can be performed.
    params_and_fit = np.array(params_and_fit)
    row = minimiser(params_and_fit[:,3])
    best_red_vec_mag, best_red_vec_x, best_red_vec_y, chi_squ_min = params_and_fit[row,:]

    de_reddened_gr_excess = gr_excess + best_red_vec_x
    de_reddened_ug_excess = ug_excess + best_red_vec_y

    g_abs = best_red_vec_y * (1 - ((g_lambda/u_lambda) - 1)**-1)
    u_abs = best_red_vec_y + g_abs
    r_abs = g_abs - best_red_vec_x

    de_reddened_r_mag = r_mag + r_abs
    de_reddened_g_mag = g_mag + g_abs
    de_reddened_u_mag = u_mag + u_abs

    # All dereddening calculations complete, write out the catalogue.
    write_cat(de_reddened_r_mag, de_reddened_g_mag, de_reddened_u_mag,
              "de_reddened_combined")

    print("Best chi-squared: {}, for reddening vector with components {}E(g-r)+{}E(u-g)".format(chi_squ_min, best_red_vec_x, best_red_vec_y))

    dict = {}
    dict["M52 Uncorrected"] = (gr_excess,ug_excess,"o")
    dict["M52 De-reddened"] = (de_reddened_gr_excess,de_reddened_ug_excess,"o")
    dict["Pleiades Fit"] = (np.linspace(-0.6,0.6,1000),polynomial(np.linspace(-0.6,0.6,1000), pleiades_coeffs),"-")

    plot_diagram(dict, x_label="Colour:(g-r)", y_label="Colour:(u-g)",
                 sup_title="M52\nColour-Colour Diagram",
                 legend=True, filename="M52_Colour-Colour_Diagram"
                )

main()
