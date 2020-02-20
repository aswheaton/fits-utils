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

    u_abs = best_red_vec_y * (1 - ((u_lambda/g_lambda) - 1)**-1)
    de_reddened_u_mag = u_mag + u_abs

# Begin second reddening vector determination, which calculates the de-reddened
# colour excess in the g-u and u-r, and the absorption in the g-band.

    # Pleaides fit values in E(g-r) vs. E(r-u) colour-colour space.
    A, B, C, D, E = 1.283, -0.107, -2.748, 8.477, -3.141
    pleiades_coeffs = [A, B, C, D, E]

    # Calculate the colour excess.
    ru_excess = r_mag - u_mag # x-axis variable
    gr_excess = g_mag - r_mag # y-axis variable
    # Calculate error on colour excess.
    ru_excess_err = r_err + u_err
    gr_excess_err = g_err + r_err

    # 2d Array containing x,y components of reddening vector and the
    # corresponding chi-squared value for that shift.
    params_and_fit = np.zeros((10000,3))
    row = 0
    for red_vec_x in np.linspace(-2.0, -0.01, 100):
        for red_vec_y in np.linspace(-2.0, -0.01, 100):
            ru_excess_shifted = ru_excess + red_vec_x
            gr_excess_shifted = gr_excess + red_vec_y
            chi_squ = get_chi_squ(ru_excess, gr_excess, ru_excess_shifted,
                                  gr_excess_shifted, polynomial, pleiades_coeffs,
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
    de_reddened_ru_excess = ru_excess + best_red_vec_x
    de_reddened_gr_excess = gr_excess + best_red_vec_y

    g_abs = best_red_vec_y * (1 - ((g_lambda/r_lambda) - 1)**-1)
    de_reddened_g_mag = g_mag + g_abs

# Begin third reddening vector determination, which calculates the de-reddened
# colour excess in the r-u and u-g, and the absorption in the g-band.

    # Pleaides fit values in E(r-u) vs. E(u-g) colour-colour space.
    A, B, C, D, E = 1.283, -0.107, -2.748, 8.477, -3.141
    pleiades_coeffs = [A, B, C, D, E]

    # Calculate the colour excess.
    ug_excess = u_mag - g_mag # x-axis variable
    ru_excess = r_mag - u_mag # y-axis variable
    # Calculate error on colour excess.
    ug_excess_err = u_err + g_err
    ru_excess_err = r_err + u_err

    # 2d Array containing x,y components of reddening vector and the
    # corresponding chi-squared value for that shift.
    params_and_fit = np.zeros((10000,3))
    row = 0
    for red_vec_x in np.linspace(-2.0, -0.01, 100):
        for red_vec_y in np.linspace(-2.0, -0.01, 100):
            ug_excess_shifted = ug_excess + red_vec_x
            ru_excess_shifted = ru_excess + red_vec_y
            chi_squ = get_chi_squ(ug_excess, ru_excess, ug_excess_shifted,
                                  ru_excess_shifted, polynomial, pleiades_coeffs,
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
    de_reddened_ug_excess = ug_excess + best_red_vec_x
    de_reddened_ru_excess = ru_excess + best_red_vec_y

    r_abs = best_red_vec_y * (1 - ((r_lambda/u_lambda) - 1)**-1)
    de_reddened_r_mag = r_mag + r_abs

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
