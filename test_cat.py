from fits_utils import *

# Zero points from zero-point-calculator.
zpr, zpg, zpu = get_zero_points(1.04) # Rory's ZP: 30.236, 29.719, 27.075
print("zpr = {}, zpg = {}, zpu = {}".format(str(zpr)[:5], str(zpg)[:5], str(zpu)[:5]))
# Get the zero point corrected catalogue and error.
exp_cor_r_mag, _, exp_cor_g_mag, _, exp_cor_u_mag, _ = load_cat("cat/exp_corr/ugr.cat", zpr, zpg, zpu)
raw_r_mag, _, raw_g_mag, _, raw_u_mag, _ = load_cat("cat/ugr.cat", zpr, zpg, zpu)

bins = np.linspace(0, 10000, 100)

exp_cor_r_counts = 10 ** ( (exp_cor_r_mag - zpr) / -2.5)
raw_r_counts = 10 ** ( (raw_r_mag - zpr) / -2.5)
plt.hist([exp_cor_r_counts, raw_r_counts/600], bins, label=['Exposure Corrected c/s', 'Raw Stack Counts / 600s'])
plt.title('Counts/Second in r')
plt.legend(loc='upper right')
plt.show()

exp_cor_g_counts = 10 ** ( (exp_cor_g_mag - zpg) / -2.5)
raw_g_counts = 10 ** ( (raw_g_mag - zpg) / -2.5)
plt.hist([exp_cor_g_counts, raw_g_counts/600], bins, label=['Exposure Corrected c/s', 'Raw Stack Counts / 600s'])
plt.title('Counts/Second in g')
plt.legend(loc='upper right')
plt.show()

bins = np.linspace(0, 5000, 100)

exp_cor_u_counts = 10 ** ( (exp_cor_u_mag - zpu) / -2.5)
raw_u_counts = 10 ** ( (raw_u_mag - zpu) / -2.5)
plt.hist([exp_cor_u_counts, raw_u_counts/600], bins, label=['Exposure Corrected c/s', 'Raw Stack Counts / 600s'])
plt.title('Counts/Second in u')
plt.legend(loc='upper right')
plt.show()
