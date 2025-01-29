import numpy as np
import h5py
import pandas as pd
import os
import sys
import scipy
import math
script_R = 3.36*10**3
T = 10 ##Kelvin
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.989*10.0**30 ## kg 
G = 6.674*10.0**(-11)  ## m^3*kg^-1*s^-2
a = 0.975471932310942752606143078377

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

v_x = np.logspace(-12, np.log10(2), 10000)
v_integral = np.loadtxt('/home/rbanderson/projectfiles/bminbmaxShu/shuInt.txt')

IntHelper = scipy.interpolate.interp1d(v_x, v_integral, kind = 'cubic')

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

Ncapture = 0

for particle in simdat['particle'].unique():
    if any(energy < 0 for energy in simdat['energy']):
        Ncapture +=1
    V = 0
    t_f = 350000 * 3.154e7
    for v0 in simdat['v0'].unique():
        b = get_bmax(c_s*t_f, v0, t_f)
        V += Ncapture / getn0(t_f, v0, b)

if V != 0:
    with open('capturevoloutput.txt', 'a') as f:
        f.write(f"{v0} {V}\n")
