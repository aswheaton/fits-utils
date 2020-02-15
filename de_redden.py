from fits_utils import *

def main():

    # Pleaides fit values.
    A, B, C, D, E = 1.283, -0.107, -2.748, 8.477, -3.141
    pleiades_coeffs = [A, B, C, D, E]
    # Zero points from zero-point-calculator.
    zpr, zpg, zpu = 30.0, 30.0, 27.0
    # Get the zero point corrected catalog and error.
    r_mag, r_err, g_mag, g_err, u_mag, u_err = load_cat('cat/combined.cat', zpr, zpg, zpu)
    sqr_err = g_err**2 + r_err**2 + u_err**2
    # Calculate the colour excess.
    gr_excess = g_mag - r_mag
    ug_excess = u_mag - g_mag
    # Calculate error on colour excess.
    gr_excess_err = g_err + r_err
    ug_excess_err = u_err + g_err

    plt.scatter(ug_excess,gr_excess)
    plt.title("Reddened U-V vs. G-R Colour Excess")
    plt.show()

    # 2d Array containing x,y components of reddening vector and the
    # corresponding chi-squared value for that shift.
    params_and_fit = np.zeros((10000,3))
    row = 0
    for red_vec_x in np.linspace(-10.0, 10.0, 100):
        for red_vec_y in np.linspace(-10.0, 10.0, 100):
            gr_excess_shifted = gr_excess + red_vec_x
            ug_excess_shifted = ug_excess + red_vec_y
            chi_squ = get_chi_squ(gr_excess_shifted, ug_excess_shifted,
                                  ug_excess_err, polynomial, pleiades_coeffs
                                 )
            params_and_fit[row,0] = red_vec_x
            params_and_fit[row,1] = red_vec_y
            params_and_fit[row,2] = chi_squ
            row += 1

    best_row = minimiser_2(params_and_fit[:,2])
    best_fit = np.array(params_and_fit[best_row,:])
    best_red_vec_x, best_red_vec_y, chi_squ_min = best_fit[0,0], best_fit[0,1], best_fit[0,2]

    print('Best chi-squared: {}, for reddening vector with components {}, {}'.format(chi_squ_min, best_red_vec_x, best_red_vec_y))

    de_reddened_gr_excess = gr_excess + best_red_vec_x
    de_reddened_ug_excess = ug_excess + best_red_vec_y

    plt.scatter(gr_excess, ug_excess)
    plt.scatter(de_reddened_gr_excess, de_reddened_ug_excess)
    plt.scatter(gr_excess, polynomial(gr_excess, pleiades_coeffs))
    plt.title("Dereddened U-G vs. G-R Colour Excess")
    plt.show()

    # m= 0.9199548
    # ms=(m**2)+1
    # delrelx= (1/ms)**0.5
    # delrely=(m)/((ms)**0.5)
    # a= -6.9834437936715
    # b=-0.36562274114452
    # c=1.3158896161408
    # k2 = m+b
    # ll= int(input("Please enter a lower limit for the reddening vector: "))
    # ul= int(input("Please enter an upper limit for the reddening vector: "))
    # dr= int(input("Please enter a step-size: "))

    # lambda_fit = []
    # for i in np.linspace(ll, ul, num=(ul-ll)/dr):
    #     delta_x, delta_y = i * delrelx, i* delrely
    #     chi_sqr = 0.0
    #     correct_ug_excess = ug_excess + delta_y
    #     correct_gr_excess = gr_excess - delta_x
    #     k1 = (c+(correct_ug_excess)) - (m*correct_gr_excess)
    #     #: Only care about negative solution.
    #     ex_x = (k2-((k2**2)+(4*a*k1))**0.5)/(2*a)
    #     ex_y = a*(ex_x**2) + b*ex_x + c
    #     r_sq = (correct_gr_excess - ex_x)**2 + (correct_ug_excess - ex_y)**2
    #     chi_sq = r_sq/sqr_err
    #     fit = np.array([np.abs(i), np.sum(chi_sq)])
    #     lambda_fit.append(fit)
    # lambda_fit = np.array(lambda_fit)
    # best_fit = minimiser(lambda_fit)
    # print('Best a lambda: {}'.format(best_fit))
    # de_red_ug = ug_excess + delrely * best_fit
    # de_red_gr = gr_excess - delrelx * best_fit

main()
