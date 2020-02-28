from fits_utils import *
import matplotlib.pyplot as plt

def main():

    cat_dir = "cat/cumulative_trim/"
    # Central wavelength of filters, in micrometers.
    r_lambda, g_lambda, u_lambda = 0.6231, 0.4770, 0.3543
    # Zero points from zero-point-calculator.
    zpr, zpg, zpu = get_zero_points(1.04)
    print("zpr = {}, zpg = {}, zpu = {}".format(str(zpr)[:5], str(zpg)[:5], str(zpu)[:5]))
    # Get the zero point corrected catalogue and error.
    r_mag, r_err, g_mag, g_err, u_mag, u_err = load_cat(cat_dir+"ugr.cat", zpr, zpg, zpu)

    pleiades_data = np.loadtxt("pleiades/pleiades_johnson.txt")
    pleiades_data = correct_pleiades(pleiades_data)
    indices = np.where(pleiades_data[:,0] < 0.5)[0]
    reduced_data = pleiades_data[indices,:]
    pleiades_coeffs, cov = np.polyfit(reduced_data[:,0], reduced_data[:,1], deg=4, cov=True)
    pleiades_coeffs = np.flip(pleiades_coeffs)
    squ_pleiades_cov = np.flip(np.diag(cov))
    # Calculate the colour excess.
    gr_excess = g_mag - r_mag # x-axis variable
    ug_excess = u_mag - g_mag # y-axis variable
    # Calculate error on colour excess.
    squ_gr_excess_err = g_err**2 + r_err**2
    squ_ug_excess_err = u_err**2 + g_err**2

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

    # Remove outliers which might skew the chi-squared minimisation.
    reduced_indices = remove_outliers(gr_excess, 0.2, 0.7)
    # Iterate over reasonable range of values for the reddening vector magnitude.
    start, stop, steps = 0.6, 1.5, 1000
    for red_vec_mag in np.linspace(start, stop, steps):
        squ_err_mag = ((stop - start) / steps)**2
        # Separate reddening vector into components in colour-colour space.
        red_vec_x = (red_vec_mag**2 / (1 + cardelli_slope**2))**0.5
        red_vec_y = (red_vec_mag**2 / (1 + cardelli_slope**-2))**0.5
        squ_err_x = (1.0 + cardelli_slope**2)**-1 * 2 * red_vec_mag * squ_err_mag
        squ_err_y = (1.0 + cardelli_slope**-2)**-1 * 2 * red_vec_mag * squ_err_mag
        # Move positions in colour-colour space along the direction of the reddening vector.
        gr_excess_shifted = gr_excess[reduced_indices] - red_vec_x
        ug_excess_shifted = ug_excess[reduced_indices] - red_vec_y
        squ_gr_shifted_err = squ_gr_excess_err[reduced_indices] + squ_err_x
        squ_ug_shifted_err = squ_ug_excess_err[reduced_indices] + squ_err_y
        # Get the minimum chi-squared value for the particular magnitude.
        chi_squ = get_chi_squ(gr_excess_shifted, squ_gr_shifted_err,
                              ug_excess_shifted, squ_ug_shifted_err,
                              polynomial, pleiades_coeffs, squ_pleiades_cov
                              )
        params_and_fit.append([red_vec_mag, squ_err_mag,
                               red_vec_x, squ_err_x,
                               red_vec_y, squ_err_y,
                               chi_squ])

    # Dirty data type change so minimisation can be performed.
    params_and_fit = np.array(params_and_fit)
    row = minimiser(params_and_fit[:,6])
    best_red_vec_mag, squ_err_mag, best_red_vec_x, squ_err_x, best_red_vec_y, squ_err_y, chi_squ_min = params_and_fit[row,:][0]

    g_abs = best_red_vec_y * ((g_lambda/u_lambda) - 1)**-1
    u_abs = best_red_vec_y + g_abs
    r_abs = g_abs - best_red_vec_x
    squ_err_g = ((g_lambda/u_lambda) - 1)**-2 * squ_err_y
    squ_err_u = squ_err_y + squ_err_g
    squ_err_r = squ_err_x + squ_err_g

    v_abs_g = g_abs / cardelli_consts["g"]
    v_abs_u = u_abs / cardelli_consts["u"]
    v_abs_r = r_abs / cardelli_consts["r"]
    squ_v_abs_g_err = squ_err_g / cardelli_consts["g"]**2
    squ_v_abs_u_err = squ_err_u / cardelli_consts["u"]**2
    squ_v_abs_r_err = squ_err_r / cardelli_consts["r"]**2

    print("Cardelli slope: {}".format(cardelli_slope))
    print("Magnitude: {} +/- {}".format(best_red_vec_mag, squ_err_mag**0.5))
    print("x-comp: {}\ny-comp: {}".format(best_red_vec_x, best_red_vec_y))
    print("Chi-Squ: {}".format(chi_squ_min))
    print("A_g = {} +/- {}".format(g_abs, squ_err_g**0.5))
    print("A_u = {} +/- {}".format(u_abs, squ_err_u**0.5))
    print("A_r = {} +/- {}".format(r_abs, squ_err_r**0.5))
    print("A_v = {} +/- {} (from A_g)".format(v_abs_g, squ_v_abs_g_err**0.5))
    print("A_v = {} +/- {} (from A_u)".format(v_abs_u, squ_v_abs_u_err**0.5))
    print("A_v = {} +/- {} (from A_r)".format(v_abs_r, squ_v_abs_r_err**0.5))

    # Create de-reddened ugr.cat catalogue with errors.
    new_catalogue = np.column_stack((u_mag - u_abs, u_err**2 + squ_err_u,
                                     g_mag - g_abs, g_err**2 + squ_err_g,
                                     r_mag - r_abs, r_err**2 + squ_err_r
                                     ))
    header_txt = "\n".join(["[0] : DE_REDDENED_U_MAG",
                            "[1] : SQUARE_ERR_U_MAG",
                            "[2] : DE_REDDENED_G_MAG",
                            "[3] : SQUARE_ERR_U_MAG",
                            "[4] : DE_REDDENED R_MAG",
                            "[5] : SQUARE_ERR_U_MAG",
                            ])
    np.savetxt(cat_dir+"de_red_ugr.cat", new_catalogue, header=header_txt )

    de_red_gr_excess = gr_excess - best_red_vec_x
    de_red_ug_excess = ug_excess - best_red_vec_y
    squ_dered_gr_err = squ_err_g + squ_err_r + squ_err_x
    squ_dered_ug_err = squ_err_u + squ_err_g + squ_err_y

    # Plot chi-sqaured as a function of reddening vector magnitude and the
    # reddened data alongside the de-reddened data and the Pleiades data.
    colour_range = np.linspace(-0.8, 1.0, 1000)
    pleiades_curve, _ = polynomial(colour_range, 0.0, pleiades_coeffs, squ_pleiades_cov)
    dict = {"M52 Uncorrected" : (gr_excess,ug_excess,"o"),
            "M52 De-reddened" : (de_red_gr_excess,de_red_ug_excess,"o"),
            "Pleiades Data"   : (pleiades_data[:,0], pleiades_data[:,1], "o"),
            "Pleiades Fit"    : (colour_range, pleiades_curve, "-"),
            "Cardelli Slope"  : (colour_range, cardelli_slope*colour_range+y_cept, "-")
            }
    plot_diagram(dict, x_label="Colour:(g-r)", y_label="Colour:(u-g)",
                 sup_title="M52\nColour-Colour Diagram",
                 legend=True, filename="M52_Colour_Colour_Diagram"
                )
    plt.plot(params_and_fit[:,0], params_and_fit[:,6])

    # Load in the larger g and r catalogue of objects which are invisible in u.
    catalog = np.loadtxt(cat_dir+"gr.cat")
    r_mag, r_err = get_mag(catalog[:,5], catalog[:,6], zpr, 600)
    g_mag, g_err = get_mag(catalog[:,3], catalog[:,4], zpg, 600)
    # De-redden the larger catalogue with newly found r_abs and g_abs values.
    de_red_r_mag = r_mag - r_abs
    de_red_g_mag = g_mag - g_abs
    squ_dered_r_err = r_err**2 + squ_err_r
    squ_dered_g_err = g_err**2 + squ_err_g
    # Calculate the de-reddened colour excess.
    de_red_gr_excess = de_red_g_mag - de_red_r_mag
    squ_dered_gr_err = squ_dered_r_err + squ_dered_g_err

    # Write the corrected catalogue out.
    new_catalogue = np.column_stack((de_red_gr_excess, squ_dered_gr_err,
                                     de_red_r_mag, squ_dered_g_err
                                     ))
    header_txt = "\n".join(["[0] : DE_REDDENED_GR_COLOR",
                            "[1] : SQUARE_ERR_GR_COLOR",
                            "[2] : DE_REDDENED_R_MAG",
                            "[3] : SQUARE_ERR_R_MAG"
                            ])
    np.savetxt(cat_dir+"de_red_gr_r.cat", new_catalogue, header=header_txt)

    # Plot the de-reddened magnitude vs. colour diagram.
    dict = {"M52 r vs. g-r"      : (de_red_gr_excess, de_red_r_mag,"o"),
            "Pleiades r vs. g-r" : (pleiades_data[:,0], pleiades_data[:,2], "o")
           }
    plot_diagram(dict, x_label="Colour:(g-r)", y_label="Magnitude: g",
                 sup_title="M52\nColour-Magnitude Diagram",
                 legend=True, filename="M52_Colour_Magnitude_Diagram"
                )

main()
