import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize

zpu=27.075
zpg=29.719
zpr=30.236
solar_lum= 3.828E+26
solar_mass= 1.99E+30
h=6.63E-34
c=3E+8
cen_wav_r= 658E-9
dist_pl_pc=150

def get_mag(flux, flux_err, zpoint):
    mag = zpoint - (2.5*np.log(flux)/np.log(10))
    mag_err = (-2.5/np.log(10))*(flux_err/flux)
    return mag, mag_err

def sys_change(b_v,v):
    g_r= 1.02*(b_v) -0.22
    r= v- 0.42*(b_v)+0.11
    r_abs=r-5*np.log10(15)
    return(g_r,r_abs)

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

def get_distance(param,param_pl,color_gr,g_r_pl):
    start= np.min(color_gr)
    stop=0.215
    shift_list=[]
    for x in np.linspace(start,stop,num=1000):
        y_pl=polynomial(x,*param_pl)
        y_m52= polynomial(x,*param)
        y_diff=abs(y_pl -y_m52)
        shift_list.append(y_diff)

    dist_mod= np.mean(shift_list)
    #print(dist_mod)
    exponent=((dist_mod/5)+1)
    distance= 10**(exponent)
    return(distance)

def main():
    data=np.loadtxt("pl_data_UBV.txt",usecols=[2,3])
    catalog = np.loadtxt('combined.cat')
    r_flux= catalog[:,5]
    g_mag, g_err = get_mag(catalog[:,3], catalog[:,4], zpg)
    r_mag, r_err = get_mag(catalog[:,5], catalog[:,6], zpr)
    u_mag, u_err = get_mag(catalog[:,7], catalog[:,8], zpu)
    color_gr = g_mag - r_mag
    b_v= data[:,1]
    v= data[:,0]
    index= np.where(b_v<0.215)[0]
    b_v_cut=b_v[index]
    v_cut=v[index]
    g_r_pl, r_pl= sys_change(b_v_cut,v_cut)
    param_pl= get_fit(polynomial,g_r_pl,r_pl)
    param_m52= get_fit(polynomial,color_gr,r_mag)
    dist_pc=get_distance(param_m52,param_pl,color_gr,g_r_pl)
    dist_rounded=np.round(dist_pc,decimals=3)
    print("The distance to Messier 52 is: " +str(dist_rounded)+" parsecs.")
    distance=(3.08567782E+16)*dist_pc
    r_min=np.min(r_mag)
    r_min_flux=mag_convert(r_min,zpr)
    r_min_lum= r_min_flux*(4*np.pi*(distance**2))
    r_min_mass= (r_min_lum/solar_lum)**(1/3.5)
    age= (r_min_mass**-2.5)* (10**10)
    age_mil= age/1000000
    print("The age of Messier 52 is: "+str(age)+" Million years old.")


    x=np.linspace(-0.7,0.1,1000)

    plt.plot(g_r_pl,r_pl,'o',label="Pleiades literature")
    plt.plot(color_gr,r_mag,'o',label="M52 Data(Uncorrected)")
    plt.plot(x,polynomial(x,*param_pl),'-',label="Pleiades")
    plt.plot(x,polynomial(x,*param_m52),'-',label="M52")
    plt.gca().invert_yaxis()
    plt.xlabel("G-R Colour")
    plt.ylabel("R Magnitude")
    plt.legend()
    plt.show()

main()
