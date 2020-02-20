from fits_utils import *

def main():

    # Central wavelength of filters, in Angstroms.
    r_lambda, g_lambda, u_lambda = 6231, 4770, 3543
    # Zero points from zero-point-calculator.
    zpr, zpg, zpu = 30.236, 29.719, 27.075
    # Get the zero point corrected catalogue and error.
    r_mag, r_err, g_mag, g_err, u_mag, u_err = load_cat("cat/ugr.cat", zpr, zpg, zpu)
    error = (g_err**2 + r_err**2 + u_err**2)**0.5

# Begin first reddening vector determination, which calculates the de-reddened
# colour excess in the u-g and g-r, and the absorption in the g-band.

    # Pleaides fit values in E(u-g) vs. E(g-r) colour-colour space.
    A, B, C, D, E = 1.283, -0.107, -2.748, 8.477, -3.141
    pleiades_coeffs = [A, B, C, D, E]

    # Calculate the colour excess.
    gr_excess = g_mag - r_mag # x-axis variable
    ug_excess = u_mag - g_mag # y-axis variable
    # Calculate error on colour excess.
    gr_excess_err = g_err + r_err
    ug_excess_err = u_err + g_err


    # 2d Array containing x,y components of reddening vector and the
    # corresponding chi-squared value for that shift.
    params_and_fit = np.zeros((10000,3))
    row = 0
    for red_vec_x in np.linspace(-2.0, -0.01, 100):
        for red_vec_y in np.linspace(-2.0, -0.01, 100):
            gr_excess_shifted = gr_excess + red_vec_x
            ug_excess_shifted = ug_excess + red_vec_y
            chi_squ = get_chi_squ(gr_excess, ug_excess, gr_excess_shifted,
                                  ug_excess_shifted, polynomial, pleiades_coeffs,
                                  error
                                  )
            params_and_fit[row,0] = red_vec_x
            params_and_fit[row,1] = red_vec_y
            params_and_fit[row,2] = chi_squ
            row += 1

    best_row = minimiser_2(params_and_fit[:,2])
    best_fit = np.array(params_and_fit[best_row,:])
    best_red_vec_x, best_red_vec_y, chi_squ_min = best_fit[0,0], best_fit[0,1], best_fit[0,2]

    print("Best chi-squared: {}, for reddening vector with components {}E(g-r)+{}E(u-g)".format(chi_squ_min, best_red_vec_x, best_red_vec_y))
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

    # dict = {}
    # dict["M52 Uncorrected"] = (gr_excess,ug_excess,"o")
    # dict["M52 De-reddened"] = (de_reddened_gr_excess,de_reddened_ug_excess,"o")
    # dict["Pleiades Fit"] = (np.linspace(-0.6,0.6,1000),polynomial(np.linspace(-0.6,0.6,1000), pleiades_coeffs),"-")
    #
    # plot_diagram(dict, x_label="Colour:(g-r)", y_label="Colour:(u-g)",
    #              sup_title="M52\nColour-Colour Diagram",
    #              legend=True, filename="M52_Colour-Colour_Diagram"
    #             )
main()
