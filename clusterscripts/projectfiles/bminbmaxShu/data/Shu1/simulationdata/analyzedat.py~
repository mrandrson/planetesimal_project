import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy
import math
import sys
import os
script_R = 3.36*10**3
T = 10 ##Kelvin
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.9891*10.0**30
G = 6.67*10.0**(-11)  ## m^3*kg^-1*s^-2
a = 0.975471932310942752606143078377
rmax = 1.974342*10.0**15
file_path = sys.argv[1]

filename = os.path.splitext(os.path.basename(file_path))[0]

def years_to_seconds(years):
    seconds = years * 31557600
    return seconds

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
    centralMass = .9754516*c_s**2*r0/G
    Mcalc =  centralMass + r0**3/(G*t**2)*get_x_integral(r/r0)
    return min(Mcalc, M)

def getPhi(r, t):
    rmax = 1.974342*10.0**15
    phi_max = -G*M/rmax
    r_meters = r
    if r >= rmax:
        return -G*M/r_meters
    else:
        return phi_max - G*trapIntegrateLog(lambda rp: (Mtot(rp, t))/rp**2, r_meters, rmax, 1000)
    
def getEnergy(r, v, t):
    return getPhi(r, t)+0.5*v**2

df1 = pd.read_csv(f'{file_path}')

df1['energy error'] = np.nan

for particle in df1['particle_id'].unique():
    particlemask = df1['particle_id'] == particle
    filtered_df = df1[particlemask].reset_index()
    r = np.sqrt(filtered_df['x']**2+filtered_df['y']**2)
    v = np.sqrt(filtered_df['vx']**2+filtered_df['vy']**2)
    energyerror = np.array([
        filtered_df['energy'][i] - getEnergy(r[i], v[i], filtered_df['time'][i])
        for i in range(len(filtered_df))
    ])

    df1.loc[particlemask, 'energy error'] = energyerror

df1.to_csv(f'./energyerrorout/{filename}energyerror.csv', index=False)
