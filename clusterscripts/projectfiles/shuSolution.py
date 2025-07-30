import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from scipy.optimize import brentq
import h5py
import math
import rebound
import sys
import time
import os

script_R = 3.36*10**3
T = 10 ##Kelvin
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.989*10.0**30 ## kg 
G = 6.674*10.0**(-11)  ## m^3*kg^-1*s^-2
a = 0.975471932310942752606143078377

Ainput = float(sys.argv[1])
Ninput = int(sys.argv[2])

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

def takeStep(alpha, v, x, dx):
    vp = (alpha*(x-v)-2.0/x)*(x-v)/((x-v)**2-1)
    alphap = (alpha-2.0/x*(x-v))*(x-v)*alpha/((x-v)**2-1)
    vt = v + vp*dx
    alphat = alpha + alphap*dx
    xt = x + dx
    vp1 = (alphat*(xt-vt)-2.0/xt)*(xt-vt)/((xt-vt)**2-1)
    alphap1 = (alphat-2.0/xt*(xt-vt))*(xt-vt)*alphat/((xt-vt)**2-1)
    vp = 0.5*(vp + vp1)
    alphap = 0.5*(alphap + alphap1)
    return v + vp*dx, alpha + alphap*dx

def getSolution(A, Nsteps):
    xvec = np.logspace(np.log10(2), -12.0, Nsteps)
    vvec = np.zeros(Nsteps)
    alphavec = np.zeros(Nsteps)
    vvec[0] = (-(A-2)/xvec[0])-((1-A/6)*(A-2))/(xvec[0]**3)
    alphavec[0] = A/xvec[0]**2-A*(A-2)/(2*xvec[0]**4)
    for i in range(0, Nsteps-1):
        alpha = alphavec[i]
        x = xvec[i]
        v = vvec[i]
        dx = xvec[i+1]-xvec[i]
        vvec[i+1], alphavec[i+1] = takeStep(alpha, v, x, dx)
    return xvec, alphavec, vvec

x4, alpha4, v4 = getSolution(Ainput, Ninput)


alphaHelp = scipy.interpolate.interp1d(x4, alpha4)
vHelp = scipy.interpolate.interp1d(x4, v4)

with h5py.File(f"shusolutionoutput{Ainput}.h5", 'w') as hf:
    hf.create_dataset('alpha', data=alpha4)
    hf.create_dataset('x', data=x4)
    hf.create_dataset('v', data=v4)
 
print(f"Data saved to shusolutionoutput{Ainput}.h5")
