import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
import pandas as pd
from scipy.optimize import brentq
import h5py
import math
import rebound
import sys
import time
import pandas as pd
import os
import matplotlib.pyplot as plt
script_R = 3.36*10**3
T = 10 ##Kelvin
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.989*10.0**30 ## kg 
G = 6.674*10.0**(-11)  ## m^3*kg^-1*s^-2
a = 0.975471932310942752606143078377
tf = 35000*3.156e7
#^ For R_fin = Rsun/10
class SplitArray:
    def __init__(self, array, m):
        self.array = array
        self.m = m
        self.n = len(array)
        self.chunk_size = self.n // m
        self.remainder = self.n % m

    def __getitem__(self, index):
        if index >= self.m or index < -self.m:
            raise IndexError("Split index out of range")

        if index < 0:
            index += self.m

        start = index * self.chunk_size + min(index, self.remainder)
        end = start + self.chunk_size
        if index < self.remainder:
            end += 1
        return self.array[start:end]

    def __len__(self):
        return self.m

Ncapt = {}

for f in os.listdir('.'):
    if f.endswith('.csv'):
        dat = pd.read_csv(f)
        v_0 = float(dat['v0'].unique()[0])
        Ncapt[f'{v_0}'] = 0
        for p in dat['particle'].unique():
            pdat = dat[dat['particle'] == p]
            if any(pdat['energy'] < 0):
                Ncapt[f'{v_0}'] +=1

print({k: v for k, v in Ncapt.items() if v > 0})

def trapIntegrateLog(f, xmin, xmax, N):
    s = np.logspace(np.log10(xmin), np.log10(xmax), N)
    fvec = np.zeros(N)
    m = 0
    for i in range(0, N):
        fvec[i] = f(s[i])
    for i in range(0, N-1):
        deltax = s[i+1] - s[i]
        av = (fvec[i+1] + fvec[i])/2
        m = m + deltax*av
    return m

def trapIntegrateLinear(f, xmin, xmax, N):
    s = np.linspace(xmin, xmax, N)
    fvec = np.zeros(N)
    m = 0
    for i in range(0, N):
        fvec[i] = f(s[i])
    for i in range(0, N-1):
        deltax = s[i+1] - s[i]
        av = (fvec[i+1] + fvec[i])/2
        m = m + deltax*av
    return m

v_x = np.logspace(-12, np.log10(2), 10000)

v_integral = np.loadtxt('/Users/richardanderson/workdir/planetesimal_project/shuInt.txt')

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

def getPhi(r, t):
    rmax = (G*M)/(2*c_s**2)
    phi_max = -G*M/rmax
    r_meters = r
    if r >= rmax:
        return -G*M/r
    else:
        return phi_max - G*trapIntegrateLog(lambda rp: (Mtot(rp, t))/rp**2, r, rmax, 10000)

def get_bmax(R, v0, t):
    phi = getPhi(R, t)
    return R*math.sqrt(1-2*phi/v0**2)

def getcaptV(v0, tf, Nsim):
    rmax = (G*M)/(2*c_s**2)
    b = get_bmax(rmax, v0, tf)
    Ncapture = Ncapt[f'{v0}']
    return Ncapture*tf*v0*np.pi*b**2/Nsim

def plotmean(x, y, nbins):
    xbins = SplitArray(x, nbins)
    
    def yinbin(bin):
        yi = np.linspace((bin)*(len(x)/len(xbins)),(bin+1)*(len(x)/len(xbins)), int(len(x)/len(xbins)+1))
        yvec = np.array([y[int(yii)] for yii in yi])
        return yvec

    ymean = [np.mean(yinbin(bindex)) for bindex in range(len(xbins)-1)]
    xmid = [(max(xbins[bindex])+min(xbins[bindex]))/2 for bindex in range(len(xbins)-1)]
    plt.plot(xmid, ymean, label = 'Richard Code Capture Volume')

keys_sorted = sorted(Ncapt.keys(), key=float)

captVolvec = np.array([
    [float(k), getcaptV(float(k), tf, 10) / (1.496e11**3)]
    for k in keys_sorted
])

#plotmean(captVolvec[:, 0]/c_s, captVolvec[:, 1], 100)

plt.scatter(captVolvec[:, 0], captVolvec[:, 1])

'''
data = np.loadtxt("../silsbeecapturevol.txt", skiprows=1)

v0_loaded = data[:, 0]
Vvol_loaded = data[:, 1]

plt.plot(v0_loaded/c_s, Vvol_loaded, label = 'Dr. Silsbee Capture Volume')
'''

'''
plt.yscale('log')
plt.xscale('log')
plt.legend()
plt.xlabel(r'$v_0/c_s$')
plt.ylabel(r'$V_{capt}\ [AU^3]$')
plt.show()
'''
