import numpy as np
import h5py
import pandas as pd
import os
import sys
import scipy
import math
import matplotlib.pyplot as plt
script_R = 3.36*10**3
T = 10 ##Kelvin
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.989*10.0**30 ## kg 
G = 6.674*10.0**(-11)  ## m^3*kg^-1*s^-2
a = 0.975471932310942752606143078377
rmax = (G*M)/(2*c_s**2)
outdir = './outputdata'

'''
index = sys.argv[1]
file_path = f'./outputdata/particle_dat{index}.csv'

try:
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        simdat = pd.read_csv(file_path)
    else:
        print(f"Skipping: {file_path} is missing or empty.")
        sys.exit(0)
except pd.errors.EmptyDataError:
    print(f"Skipping: {file_path} is empty or has no columns.")
    sys.exit(0)
except Exception as e:
    print(f"Error reading {file_path}: {e}")
    sys.exit(0)
'''

v_x = np.logspace(-12, np.log10(2), 10000)
v_integral = np.loadtxt('/Users/richardanderson/workdir/planetesimal_project/ShuScripts/shuInt.txt')

IntHelper = scipy.interpolate.interp1d(v_x, v_integral, kind = 'cubic')

def cubic_meters_to_cubic_au(volume_m3):
    AU_in_meters = 1.495978707e11
    return volume_m3 / (AU_in_meters ** 3)

def get_x_integral(x):
    if x < 10.0**(-12):
        return 0
    if x > 2:
        return IntHelper(2) + 2*(x-2)
    return IntHelper(x)

def Mtot(r, t):
    r0 = c_s*t
    centralMass = .975502*c_s**2*r0/G
    Mcalc =  centralMass + r0**3/(G*t**2)*get_x_integral(r/r0)
    return min(Mcalc, M)

def getr_max(t):
    return scipy.optimize.brentq(lambda r: (Mp(t)+getMenc(r, t)-2*10**30), .4, 20000)

def getPhi(r, t):
    rmax = (G*M)/(2*c_s**2)
    phi_max = -G*M/rmax
    r_meters = r
    if r >= rmax:
        return -G*M/r_meters
    else:
        return phi_max - G*trapIntegrateLog(lambda rp: (Mtot(rp, t))/rp**2, r_meters, rmax, 10000)

def get_bmax(R, v0, t):
    phi = getPhi(R, t)
    return R*math.sqrt(1-2*phi/v0**2)

def getcrossect(b):
    return np.pi * b**2

def getn0(tf, v0, b):
    return 1/(tf * v0 * getcrossect(b))

'''
Ncapture = 0

for particle in simdat['particle'].unique():
    if any(energy < 0 for energy in simdat['energy']):
        Ncapture +=1
    V = 0
    #b = get_bmax(rmax, v0, t_f)
    V += Ncapture / getn0(t_f, v0, b)
'''

'''
def getV(file):
    dat = pd.read_csv(f'{fpath}')
    Ncapture = 0
    for p in dat['particle'].unique():
        particle_data = dat[dat['particle'] == p]
        for e in particle_data['energy']:
            if e < 0:
                Ncapture += 0
                break
'''

capturefrac=[]
v_init=[]
captV=[]
for f in os.listdir(outdir):
    Ncapture = 0
    Ntot = 0
    V=0
    if f.endswith('.csv'):
        fpath = os.path.join(outdir, f)
        try:
            dat = pd.read_csv(f'{fpath}')

            for p in dat['particle'].unique():
                Ntot+=1
                pdata = dat[dat['particle'] == p]

                if any(pdata['energy'] < 0):
                    #print(f"File: {f}, Particle: {p}, Negative Energy: {e}")
                    Ncapture+=1
                
                tf = 350000*3.156e7
                b = pdata['impact parameter'].unique()
                v0 = pdata['v0'].unique()
                #b = get_bmax(rmax, v0, tf)
                V += Ncapture/getn0(tf, v0, b) 

        except pd.errors.EmptyDataError:
            #print(f"Skipping empty file: {fpath}")
            print(" ")
        except Exception as e:
            print(f"Error processing file {fpath}: {e}")
        
        if Ncapture !=0:
            print(f"File: {f}, Ntot: {Ntot}, Ncapture: {Ncapture}, Vcapture: {V}")
            captV.append(V)
            v_init.append(v0)

        if Ntot !=0:
            capturefrac.append(Ncapture/Ntot)

v_init = np.array(v_init)
captV = np.array(captV)
plt.yscale('log')
plt.xscale('log')
plt.ylim(1e13, 1e16)
plt.xlim(1e-1, 1e1)
plt.xlabel(r'$v_0/c_s$')
plt.ylabel(r'Capture Volume $[AU^3]$')
plt.scatter(v_init/c_s, cubic_meters_to_cubic_au(captV))
plt.show()

'''
plt.yscale('log')
plt.hist(capturefrac, bins='auto')
plt.show()
'''
