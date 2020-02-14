def main():
    catalog = np.loadtxt('cat/combined.cat')
    g_mag, g_err = get_mag(catalog[:,3], catalog[:,4], zpg)
    r_mag, r_err = get_mag(catalog[:,5], catalog[:,6], zpr)
    u_mag, u_err = get_mag(catalog[:,7], catalog[:,8], zpu)
    color_ug = u_mag - g_mag
    color_gr = g_mag - r_mag
    plot_HR()

if __name__ == '__main__':
    main()
