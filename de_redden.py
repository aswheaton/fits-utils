from fits_utils import *

def main():
    m= 0.9199548
    ms=(m**2)+1
    delrelx= (1/ms)**0.5
    delrely=(m)/((ms)**0.5)
    a= 5
    b=0.1
    c=15
    k2 = m+b
    ll= 0#int(input("Please enter a lower limit for the reddening vector: "))
    ul= 10#int(input("Please enter an upper limit for the reddening vector: "))
    dr= 0.1#int(input("Please enter a step-size: "))
    zpu= 1#float(input("What is your instrumental zero point in the u band?: "))
    zpg= 2#float(input("What is your instrumental zero point in the g band?: "))
    zpr= 3#float(input("What is your instrumental zero point in the r band?: "))
    catalog = np.loadtxt('cat/combined.cat')
    g_mag, g_err = get_mag(catalog[:,3], catalog[:,4], zpg)
    r_mag, r_err = get_mag(catalog[:,5], catalog[:,6], zpr)
    u_mag, u_err = get_mag(catalog[:,7], catalog[:,8], zpu)
    color_ug = u_mag - g_mag
    color_gr = g_mag - r_mag
    sqr_err = g_err**2 + r_err**2 + u_err**2

    lambda_fit = []
    for i in np.linspace(ll, ul, num=(ul-ll)/dr):
        delta_x, delta_y = i * delrelx, i* delrely
        chi_sqr = 0.0
        correct_color_ug = color_ug + delta_y
        correct_color_gr = color_gr - delta_x
        k1 = (c+(correct_color_ug)) - (m*correct_color_gr)
        #: Only care about negative solution.
        ex_x = (k2-((k2**2)+(4*a*k1))**0.5)/(2*a)
        ex_y = a*(ex_x**2) + b*ex_x + c
        r_sq = (correct_color_gr - ex_x)**2 + (correct_color_ug - ex_y)**2
        chi_sq = r_sq/sqr_err
        fit = np.array([np.abs(i), np.sum(chi_sq)])
        lambda_fit.append(fit)
    lambda_fit = np.array(lambda_fit)
    best_fit = minimiser(lambda_fit)
    print('Best a lambda: {}'.format(best_fit))
    de_red_ug = color_ug + delrely * best_fit
    de_red_gr = color_gr - delrelx * best_fit

main()
