import numpy as np
import matplotlib.pyplot as plt
from fits_utils import *
from scipy import optimize

zpu=27.075
zpg=29.719
zpr= 22.9
solar_lum= 3.828E+26
solar_mass= 1.99E+30
h=6.63E-34
c=3E+8
cen_wav_r= 658E-9
dist_pl_pc=150


def get_abs_mag(app_mag):
    return(app_mag-5*np.log10(dist_pl_pc/10))

def get_fit(f,x,y):
    param,cov= optimize.curve_fit(f,x,y,p0=None)
    return(param)

def polynomial(x,a,b,c,d,e):
    return(a*x**4 +b*x**3 +c*x**2 +d*x +e)

def mag_convert(m,zp):
    counts= 10**((zp-m)/2.5)
    energy_photon=h*c/cen_wav_r
    flux=counts*energy_photon
    return(flux)

def remove_outlie(g_r,r):
    index= np.where((g_r<-0.2) & (g_r>-0.7))[0]
    return(g_r[index],r[index])

def get_distance_2(g_r,r,param_pl):
    cut_off=np.mean(g_r)
    index_bin_1= np.where((g_r>np.amin(g_r))&(g_r<cut_off))[0]
    index_bin_2= np.where((g_r>cut_off)&(g_r<np.amax(g_r)))[0]
    bin_1_g_r= g_r[index_bin_1]
    bin_2_g_r= g_r[index_bin_2]
    bin_1_r=r[index_bin_1]
    bin_2_r=r[index_bin_2]
    mean_bin_1_r=np.mean(bin_1_r)
    mean_bin_2_r=np.mean(bin_2_r)
    mean_bin_1_g_r=np.mean(bin_1_g_r)
    mean_bin_2_g_r=np.mean(bin_2_g_r)
    bin_1_r_pl=polynomial(mean_bin_1_g_r,*param_pl)
    bin_2_r_pl=polynomial(mean_bin_2_g_r,*param_pl)
    r_bin_1_diff=abs(bin_1_r_pl-mean_bin_1_r)
    r_bin_2_diff=abs(bin_2_r_pl- mean_bin_2_r)
    dist_mod=(r_bin_1_diff+r_bin_2_diff)/2
    print(dist_mod)
    exponent=((dist_mod/5)+1)
    distance= 10**(exponent)
    err_dist_mod=
    err_distance=(2**((dist_mod/5)+1)*np.log(10)*(5**(dist_mod/5)))*err_dist_mod
    return(distance)


def get_distance(param,param_pl,color_gr,g_r_pl,cov_pl,cov_m52):
    start= np.amin(color_gr)
    stop=np.amax(color_gr)
    shift_list=[]
    pl_1=((g_r_pl**4)*err_a)**2
    pl_2=((g_r_pl**3)*err_b)**2
    pl_3=((g_r_pl**2)*err_c)**2
    pl_4=((g_r_pl)*err_d)**2
    pl_5=(err_e)**2
    pl_6=([4*param_pl[0]*(g_r_pl**3) + 3*param_pl[1]*(g_r_pl)**2 + 2*param_pl[2]*(g_r_pl) + param_pl[3]]*err_g_r_pl)**2
    for x in np.linspace(start,stop,num=1000):
        y_pl=polynomial(x,*param_pl)
        y_pl_err=np.sqrt(pl_1 +pl_2 +pl_3 +pl_4 +pl_5 +pl_6)        
        y_m52= polynomial(x,*param)
        y_diff=abs(y_pl -y_m52)
        shift_list.append(y_diff)

    dist_mod=np.mean(shift_list)
    print(dist_mod)
    err_dist_mod=
    err_distance=(2**((dist_mod/5)+1)*np.log(10)*(5**(dist_mod/5)))*err_dist_mod
    exponent=((dist_mod/5)+1)
    distance= 10**(exponent)
    return(distance)

def main():
    data=np.loadtxt("pleiades/pleiades_johnson.txt")
    data=correct_pleiades(data)
    catalog = np.loadtxt('cat/de_reddened_gr_r.cat')
    r_mag = catalog[:,1]
    color_gr = catalog[:,0]
    g_r_pl=data[:,0]
    r_pl=data[:,2]
    r_pl_abs= get_abs_mag(r_pl)
    cor_g_r_pl,cor_r_abs_pl= remove_outlie(g_r_pl,r_pl_abs)
    cor_g_r_m52,cor_r_m52= remove_outlie(color_gr,r_mag)
    param_pl,cov_pl= get_fit(polynomial,cor_g_r_pl,cor_r_abs_pl)
    param_m52,cov_m52= get_fit(polynomial,cor_g_r_m52,cor_r_m52)
    dist_pc=get_distance(param_m52,param_pl,cor_g_r_m52,cor_g_r_pl,cov_pl,cov_m52)
    #dist_pc=get_distance_2(cor_g_r_m52,cor_r_m52,param_pl)
    dist_rounded=np.round(dist_pc,decimals=3)
    print("The distance to Messier 52 is: " +str(dist_rounded)+" parsecs.")
    distance=(3.08567782E+16)*dist_pc
    r_min=np.min(cor_r_m52)
    r_min_flux=mag_convert(r_min,zpr)
    r_min_lum= r_min_flux*(4*np.pi*(distance**2))
    r_min_mass= (r_min_lum/solar_lum)**(1/3.5)
    age= (r_min_mass**-2.5)* (10**10)
    age_mil= age/1000000
    print("The age of Messier 52 is: "+str(age_mil)+" Million years old.")


    x=np.linspace(np.amin(cor_g_r_m52),np.amax(cor_g_r_m52),1000)

    plt.plot(cor_g_r_pl,cor_r_abs_pl,'o',label="Pleiades literature")
    plt.plot(cor_g_r_m52,cor_r_m52,'o',label="M52 Data(Uncorrected)")
    plt.plot(x,polynomial(x,*param_pl),'-',label="Pleiades")
    plt.plot(x,polynomial(x,*param_m52),'-',label="M52")
    plt.gca().invert_yaxis()
    plt.xlabel("G-R Colour")
    plt.ylabel("R Magnitude")
    plt.legend()
    plt.show()

main()
