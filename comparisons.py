import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd

data = pd.read_csv("PS_2023.02.01_18.28.47.csv",comment='#')

def names(dataframe):
    for p in dataframe.keys():
        print(format(p))

parameters = ['pl_name',
            'pl_bmassj','pl_radj','st_mass','st_rad','pl_eqt',
            'pl_orbper','pl_orbsmax','pl_orbeccen','pl_orbincl',
            'st_met','st_teff','st_age','pl_trandep','st_spectype']
parameter_high = ['pl_name',
                  'pl_bmassjerr1','pl_radjerr1','st_masserr1','st_raderr1','pl_eqterr1',
                  'pl_orbpererr1','pl_orbsmaxerr1','pl_orbeccenerr1','pl_orbinclerr1',
                 'st_metfeerr1','st_tefferr1','st_ageerr1','pl_trandeperr1']
parameter_low = ['pl_name',
                 'pl_bmassjerr2','pl_radjerr2','st_masserr2','st_raderr2','pl_eqterr2',
                  'pl_orbpererr2','pl_orbsmaxerr2','pl_orbeccenerr2','pl_orbinclerr2',
                 'st_metfeerr2','st_tefferr2','st_ageerr2','pl_trandeperr2']

det_type=data.loc[:,'discoverymethod']
trans = data[det_type == 'Transit']
rv = data[det_type == 'Radial Velocity']
img = data[det_type == 'Imaging']
ttv=data[det_type =='Transit Timing Variations']
lens=data[det_type =='Microlensing'] 
# timing includes TTVs, Eclipse Timing Variations, Pulsar Timing, and Pulsar Timing Variations,
timing=(data["discoverymethod"]== "Pulsar Timing")|(data["discoverymethod"]== "Pulsation Timing Variations")|(data["discoverymethod"]== "Eclipse Timing Variations")|(data["discoverymethod"]== "Transit Timing Variations")

#orbital params
per=data['pl_orbper'] #days
sma=data['pl_orbsmax'] #AU
ecc=data['pl_orbeccen']
inc=data['pl_orbincl']
# planet params
mp_earth=data['pl_bmasse'] # M_earth
mp=data['pl_bmassj'] # M_jupiter
rp=data['pl_radj'] # R_jupiter
# stellar params
mstar=data['st_mass'] # M_sun
rstar=data['st_rad'] # R_sun
stype=data['st_spectype']
feh=data['st_met']
teff=data['st_teff']

MEARTH = 5.9742e24 #kg 
REARTH = 6.3781e6 #km
MJUP = 1.8986e27 #kg
RJUP = 6.9911e11 #km

# Solar System params
ssmass = np.array([0.3302, 4.8685, 5.9736, 0.64185, 1898.6, 568.46, 86.832, 102.43]) * 10**24 #kg
ssper = np.array([88.0,224.7,365.2,687.0,4331,10747,30589,59800]) #days
ssecc= [0.205,0.007,0.017,0.094,0.049,0.057,0.046,0.011]
ssrad = [2439.7,6051.8,6371.00,3389.50,69911,58232,25362,24622] #km
ssrad_pm = [1.0,1.0,0.01,0.2,6,6,7,19]

multi=(data["sy_pnum"]> 1) & (data["discoverymethod"]== "Transit")
single=(data["sy_pnum"]== 1) & (data["discoverymethod"]== "Transit")
other=(data["discoverymethod"]!= "Transit")# (data["pl_discmethod"]== "Transit Timing Variations")
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 12
        }
def plot_per_vs_mp():
fig,ax = plt.subplots(figsize=(8,6))
ax.set_xscale('log') 
ax.set_yscale('log')
plt.scatter(data[other].pl_orbper,data[other].pl_bmasse,label="Other Disc. Method",marker=".",color='gray',alpha=0.45)
plt.scatter(data[single].pl_orbper,data[single].pl_bmasse,label="Single Transiting System",marker=".",color="black")
plt.scatter(data[multi].pl_orbper,data[multi].pl_bmasse,label="Multi-Transiting System",marker="o",color='white',edgecolors='black')
plt.scatter(ssper,ssmass/MEARTH,marker="o",color="red")
plt.text(300,0.6,"Venus",fontdict=font)
plt.text(480,0.88,"Earth",fontdict=font)
plt.text(120,0.045,"Mercury",fontdict=font)
plt.text(990,0.095,"Mars",fontdict=font)
plt.text(5200,270,"Jupiter",fontdict=font)
plt.text(14000,80,"Saturn",fontdict=font)
plt.text(2800,13,"Uranus",fontdict=font)
plt.text(81000,15,"Neptune",fontdict=font)
plt.xlabel("Orbital Period [days]",fontsize='large')
plt.ylabel("Mass or $M\sin{i}$ [M$_{Earth}$]",fontsize='large')
plt.tick_params(axis='both',which='both',direction='in',top=True,right=True)
ax.minorticks_on()
plt.legend(loc=4,fontsize='medium')
plt.show()

def plot_mp_vs_rp():
fig,ax = plt.subplots(figsize=(8,6))
ax.set_xscale('log') 
ax.set_yscale('log')
plt.scatter(data[other].pl_radj,data[other].pl_bmassj,label="Other Disc. Method",marker=".",color='gray',alpha=0.45)
plt.scatter(data[single].pl_radj,data[single].pl_bmassj,label="Single Transiting System",marker=".",color="black")
plt.scatter(data[multi].pl_radj,data[multi].pl_bmassj,label="Multi-Transiting 


def system_conditions(max_a,min_mp,max_mstar,max_p):
    acut = (sma < max_a) 
    mcuts = (mp > min_mp) & (mstar < max_mstar)
    pcut = (per < max_p) #very hot Jupiters < 3 day?
    full_cut = acut & mcuts
    
    #print("# of Exoplanets with sma <",max_a,'AU, mstar <',max_mstar,'M_sun, mp >',min_mp,'M_Jup ==',len(name[full_cut]))
    plt.loglog(per[pcut],mp[pcut],'.',alpha=0.75,label='Confirmed Exolanets')
    # plt.loglog(per[mcuts],mp[mcuts],'o',alpha=0.75,label='Exolanets -- mcuts')
    plt.loglog(2.30,0.779,'+',color='red',label='HAT-P-68b',markersize=15)
    plt.xlabel("Orbital Period [days]")
    plt.ylabel("Planet Mass [$M_{Jup}$]")
    plt.xlim(right=10)
    plt.legend()
    plt.show()
    #plt.savefig('hatp68_in_context.png')
    #print(name[full_cut])

#system_conditions(0.05,0.25,0.75,20)

# sol = pt.pl_bmassj*0.000954588
# ratio = sol/pt.st_mass
# print(pt.pl_name[ratio>min_depth])
# #print(pt.pl_name[delta>0.16])


def jup_to_sol(radius):
	#sol = mass * 9.55E-4 #solar masses
    sol = radius * 0.10049 #solar radii
    return print(sol)


def plot_TEPs(x,y):
    plt.loglog(pt.pl_orbper,pt.pl_bmassj,'o',alpha=0.75,label='TEP')
    # plt.loglog(rv.pl_orbper,rv.pl_bmassj,'v',alpha=0.75,label='Doppler')
# # plt.loglog(im.pl_orbper[cuts],im.pl_bmassj[cuts],'s',alpha=0.5,label='Imaging')
# # plt.loglog(ml.pl_orbper[cuts],ml.pl_bmassj[cuts],'d',alpha=0.5,label='Microlensing')
    plt.loglog(2.30,0.779,'+',label='HAT-P-68b',markersize=15)
    plt.xlabel("Orbital Period [days]")
    plt.ylabel("Planet Mass [$M_{Jup}$]")
    plt.legend()
    plt.show()

def plot_depth():
    plt.plot(pt.pl_orbper,delta,'.',alpha=0.75,label='TEP')
    plt.plot(2.30,(0.11014/0.69)**2,'+',label='HAT-P-68b',markersize=15)
    plt.xlabel("Orbital Period [days]")
    plt.ylabel("Transit Depth $\delta$")
    plt.legend()
    #plt.ylim(0.1, 0) 
    plt.show()

def planet_distribution(p_max):
    bin_size = 0.5
    period_bins = np.arange(0, p_max + bin_size, bin_size/2)
    period_range = (0,p_max) #HJ
    mp_bins = np.arange(0,2 + bin_size, bin_size/4)
    mp_range = (0,2)
    rp_bins = np.arange(0,2 + bin_size, bin_size/4)
    rp_range = (0,2)
    ecc_bins = np.arange(0,0.6+bin_size,bin_size/4)
    ecc_range = (0,1)

    #period-mass
    pm_hist,period_bins,mp_bins=np.histogram2d(per,mp,bins=(period_bins,mp_bins),range=(period_range,mp_range))
    pm_hist = pm_hist.T
    #period-radius
    pr_hist,period_bins,rp_bins=np.histogram2d(per,rp,bins=(period_bins,rp_bins),range=(period_range,rp_range))
    pr_hist = pr_hist.T    
    #period-eccentricity
    pe_hist,period_bins,ecc_bins=np.histogram2d(per,ecc,bins=(period_bins,ecc_bins),range=(period_range,ecc_range))
    pe_hist = pe_hist.T

    fig = plt.figure(figsize=(6,8))
    fig.subplots_adjust(wspace=0.25, left=0.1, right=0.95, bottom=0.07, top=0.95)
    #plt.suptitle("Hot Jupiter Distributions")
    plt.subplot(311, title = '$P-M_{p}$ Distribution')
    #plt.xlabel("Orbital Period [days]")
    plt.ylabel("Planet Mass [$M_{Jup}$]")
    im1 = plt.pcolormesh(*np.meshgrid(period_bins, mp_bins), pm_hist)
    plt.colorbar(im1)
    plt.subplot(312, title = '$P-R_{p}$ Distribution')
    #plt.ylabel("Sky Projected Obliquity [deg]")
    #plt.xlabel("Orbital Period [days]")
    plt.ylabel("Planet Radius [$R_{Jup}$]")
    im2 = plt.pcolormesh(*np.meshgrid(period_bins, rp_bins), pr_hist)
    plt.colorbar(im2)
    plt.subplot(313, title = '$P-e$ Distribution')
    #plt.ylabel("Sky Projected Obliquity [deg]")
    plt.xlabel("Orbital Period [days]")
    plt.ylabel("Planet Eccentricity")
    im3 = plt.pcolormesh(*np.meshgrid(period_bins, ecc_bins), pe_hist)
    plt.colorbar(im3)
    #plt.show()
#planet_distribution(10)

def star_distribution(p_max):
    bin_size = 0.5
    period_bins = np.arange(0, p_max + bin_size, bin_size/2)
    period_range = (0,p_max) #HJ
    feh_bins = np.arange(-2,2 + bin_size, bin_size/4)
    feh_range = (-2,2)
    rp_bins = np.arange(0,2 + bin_size, bin_size/4)
    rp_range = (0,2)
    obl_bins = np.arange(0,180+bin_size,bin_size*50)
    obl_range = (0,180)
    #period-feh
    pf_hist,period_bins,feh_bins=np.histogram2d(per,feh,bins=(period_bins,feh_bins),range=(period_range,feh_range))
    pf_hist = pf_hist.T 
    #period-obliquity
    po_hist,period_bins,obl_bins=np.histogram2d(per,mp,bins=(period_bins,obl_bins),range=(period_range,obl_range))
    po_hist = po_hist.T
    
    fig = plt.figure(figsize=(9,6))
    fig.subplots_adjust(wspace=0.25, left=0.1, right=0.95, bottom=0.07, top=0.95)
    plt.subplot(121, title = 'a.')
    plt.xlabel("Orbital Period [days]")
    plt.ylabel("Stellar Metallicity [dex]")
    im1 = plt.pcolormesh(*np.meshgrid(period_bins, feh_bins), pf_hist)
    plt.colorbar(im1)
    plt.subplot(122, title = 'b.')
    #plt.ylabel("Sky Projected Obliquity [deg]")
    plt.xlabel("Orbital Period [days]")
    plt.ylabel("Sky Projected Angle [degrees]")
    im2 = plt.pcolormesh(*np.meshgrid(period_bins, obl_bins), po_hist)
    plt.colorbar(im2)
    #plt.show()
#star_distribution(10)
