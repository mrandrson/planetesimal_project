import numpy as np
import math
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate


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

x4, alpha4, v4 = getSolution(2.0000001,100)

alphaHelp = scipy.interpolate.interp1d(x4, alpha4)
vHelp = scipy.interpolate.interp1d(x4, v4)

plt.yscale("symlog")
plt.plot(x4, alpha4, label = "$\\alpha$")
plt.plot(x4, v4, label = "$v$")
plt.legend()
plt.show()
