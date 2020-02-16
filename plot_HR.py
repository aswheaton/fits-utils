from fits_utils import *

def main():
    z_points = {'u' : 27.075, 'g' : 29.719, 'r' : 30.236}
    catalog = np.loadtxt('cat/combined.cat')
    g_mag, g_err = get_mag(catalog[:,3], catalog[:,4], z_points['g'])
    r_mag, r_err = get_mag(catalog[:,5], catalog[:,6], z_points['r'])
    u_mag, u_err = get_mag(catalog[:,7], catalog[:,8], z_points['u'])
    color_ug = u_mag - g_mag
    color_gr = g_mag - r_mag

    pleiades_data = np.loadtxt('pleiades/pleiadescutoff.txt')
    pleiades_gr, pleiades_ug = pleiades_data[:,0], pleiades_data[:,1]
    pleiades_fit = np.polyfit(pleiades_gr, pleiades_ug, 4)
    trendline = np.poly1d(pleiades_fit)
    
    plts = {'Pleiades Data' : (pleiades_gr, pleiades_ug, 'o'),
            'M52 Uncorrected' : (color_gr, color_ug, 'o'),
            'Pleiades Best Fit Line' : (np.sort(pleiades_gr),
                                        trendline(np.sort(pleiades_gr)), '-')}
    plot_diagram(plts, x_label='Color: g - r', y_label='Color: g - u',
                 sup_title='Color-Color Diagram \n M52', legend=True,
                 filename='M52_color_color_uncorrected')

    # plot_HR(color_gr, color_ug, 'Color: G-R', 'Color: U-G', 'M52')

if __name__ == '__main__':
    main()
