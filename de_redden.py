from fits_utils import *
import matplotlib.pyplot as plts

def correct_pleiades(p_data):
    #: Deredden data in johnson band.
    # p_data[:,1] = p_data[:,1] - 0.653
    # p_data[:,2] = p_data[:,2] - 0.863 + 0.653
    # p_data[:,3] = p_data[:,3] - 1.032 + 0.863
    #: ndarray: Pleiades data in sloan bands.
    p_data[:,2] = 1.28*p_data[:,2] + 1.13
    p_data[:,3] = 1.02*p_data[:,3] - 0.22
    p_data[:,2] = p_data[:,2] - 0.787 + 0.544
    p_data[:,3] = p_data[:,3] - 1.009 + 0.787

    header_txt = '[0]:No\t[1]:V\t[2]:g-r\t[3]:u-g'
    np.savetxt('pleiades/SDSS_dereddened.txt', p_data, header=header_txt)

def main():

    # Central wavelength of filters, in micrometers.
    r_lambda, g_lambda, u_lambda = 0.6231, 0.4770, 0.3543
    # Zero points from zero-point-calculator.
    zpr, zpg, zpu = 30.236, 29.719, 27.075
    # Get the zero point corrected catalogue and error.
    r_mag, r_err, g_mag, g_err, u_mag, u_err = load_cat("cat/ugr.cat", zpr, zpg, zpu)
    error = (g_err**2 + r_err**2 + u_err**2)**0.5

# Begin first reddening vector determination, which calculates the de-reddened
# colour excess in the u-g and g-r, and the absorption in the g-band.
    #: ndarray: Pleiades data.
    dr_p_data = np.loadtxt('pleiades/pleiades_johnson.txt')
    correct_pleiades(dr_p_data)
    p_data = np.loadtxt('pleiades/SDSS_dereddened.txt')
    #: ndarray: Deredden Pleiades data using extinction coeffs from NED.
    pleiades_coeffs = np.polyfit(p_data[:,2], p_data[:,3], 4)
    print(pleiades_coeffs)

    # Pleaides fit values in E(u-g) vs. E(g-r) colour-colour space.
    # A, B, C, D, E = 1.283, -0.107, -2.748, 8.477, -3.141
    # pleiades_coeffs = [A, B, C, D, E]

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
    for red_vec_mag in np.linspace(0.0,  3.0, 1000):
        red_vec_x = red_vec_mag**2 / (1 + cardelli_slope**2)
        red_vec_y = red_vec_mag**2 / (1 + cardelli_slope**-2)
        gr_excess_shifted = gr_excess - red_vec_x
        ug_excess_shifted = ug_excess - red_vec_y
        chi_squ = get_chi_squ(gr_excess, ug_excess, gr_excess_shifted,
                              ug_excess_shifted, polynomial, pleiades_coeffs,
                              error
                              )
        params_and_fit.append([red_vec_mag, red_vec_x, red_vec_y, chi_squ])

    # Dirty data type change so minimisation can be performed.
    params_and_fit = np.array(params_and_fit)
    plt.plot(params_and_fit[:,0], params_and_fit[:,3])
    row = minimiser(params_and_fit[:,3])
    best_red_vec_mag, best_red_vec_x, best_red_vec_y, chi_squ_min = params_and_fit[row,:][0]

    g_abs = best_red_vec_y * ((g_lambda/u_lambda) - 1)**-1
    u_abs = best_red_vec_y + g_abs
    r_abs = g_abs - best_red_vec_x

    v_abs_g = g_abs / cardelli_consts['g']
    v_abs_u = u_abs / cardelli_consts['u']
    v_abs_r = r_abs / cardelli_consts['r']

    print("Cardelli slope: {}\nMagnitide: {}\nx-comp: {}\ny-comp: {}\nChi-Squ: {}".format(cardelli_slope, best_red_vec_mag, best_red_vec_x, best_red_vec_y, chi_squ_min))
    print("A_g = {}\nA_u = {}\nA_r = {}".format(g_abs,u_abs,r_abs))
    print("A_v = {} (from A_g)\nA_v = {} (from A_u)\nA_v = {} (from A_r)".format(v_abs_g,v_abs_u,v_abs_r))

    de_reddened_gr_excess = gr_excess - best_red_vec_x
    de_reddened_ug_excess = ug_excess - best_red_vec_y

    de_reddened_r_mag = r_mag + r_abs
    de_reddened_g_mag = g_mag + g_abs
    de_reddened_u_mag = u_mag + u_abs

    # All dereddening calculations complete, write out the catalogue.
    write_cat(de_reddened_r_mag, de_reddened_g_mag, de_reddened_u_mag,
              "de_reddened_combined")

    dict = {}
    dict["M52 Uncorrected"] = (gr_excess,ug_excess,"o")
    dict["M52 De-reddened"] = (de_reddened_gr_excess,de_reddened_ug_excess,"o")
    value_range = np.linspace(-0.6, 0.6, 1000)
    dict["Pleiades Fit"] = (value_range,
                            np.poly1d(pleiades_coeffs)(value_range), "-")
    # dict["Pleiades Fit"] = (np.linspace(-0.6,0.6,1000),polynomial(np.linspace(-0.6,0.6,1000), pleiades_coeffs),"-")

    plot_diagram(dict, x_label="Colour:(g-r)", y_label="Colour:(u-g)",
                 sup_title="M52\nColour-Colour Diagram",
                 legend=True, filename="M52_Colour-Colour_Diagram"
                )

main()
