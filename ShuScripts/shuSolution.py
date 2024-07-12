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
    vvec[0] = (-(A-2)/xvec[0])
    alphavec[0] = A/xvec[0]**2
    for i in range(0, Nsteps-1):
        alpha = alphavec[i]
        x = xvec[i]
        v = vvec[i]
        dx = xvec[i+1]-xvec[i]
        vvec[i+1], alphavec[i+1] = takeStep(alpha, v, x, dx)
    return xvec, alphavec, vvec
'''
x1, alpha1, v1 = getSolution(2.1, 100000)
x2, alpha2, v2 = getSolution(2.01, 100000)
x3, alpha3, v3 = getSolution(2.001, 100000)
'''

x4, alpha4, v4 = getSolution(int(2+1e-6), int(1e8))

alphaHelp = scipy.interpolate.interp1d(x4, alpha4)
vHelp = scipy.interpolate.interp1d(x4, v4)
def getAlpha(x):
    return float(alphaHelp(x))

def getV(x):
    return float(vHelp(x))    
    
    
def getSlope(x):
    return (np.log(getAlpha(x*1.01)) - np.log(getAlpha(x)))/(np.log(1.01))

def analyticAlpha(x):
    c = 0.6985
    return 2/x*(1+ c/2*(1-x)**7/math.sqrt(x))

def alpha(x):
    if x<2:
        return getAlpha(x)
    else:
        return 2/x**2

def getFlux(x):
    return 4*math.pi*getV(x)*x**2*getAlpha(x)

'''
alphaAnalytical = np.zeros(1000000)
ratio = np.zeros(1000000)
for i in range(0, 1000000):
    if x4[i] > 1:
        alphaAnalytical[i] = 2.0/x4[i]**2
    else:
        alphaAnalytical[i] = analyticAlpha(x4[i])
    ratio[i] = alphaAnal[i]/alpha4[i]
'''  

fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(111)
#ax1.plot(x1, alpha1, color = 'red')
#ax1.plot(x1, -v1, color = 'red', linestyle = 'dashed')
#ax1.plot(x2, alpha2, color = 'green')
#ax1.plot(x2, -v2, color = 'green', linestyle = 'dashed')
#ax1.plot(x3, alpha3, color = 'blue')
#ax1.plot(x3, -v3, color = 'blue', linestyle = 'dashed')
#ax1.plot(x4, alpha4, color = 'violet')
#ax1.plot(x4, alphaAnal, color = 'black')
ax1.plot(x4, -v4, color = 'violet', linestyle = 'dashed')
ax1.plot(x4, ratio)
ax1.plot(x4, alpha4, color = 'red')
ax1.set_xscale("log")
ax1.set_yscale("log")

plt.show()


