import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from scipy.optimize import brentq


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

N = 100000  
x = np.zeros(N)
for i in range(0, N):
    x[i] = (1000*i)/N +.000001

eta = np.zeros(N)

eta[0] = 2/3*x[0]

xi = np.zeros(N)

xi[0] = 5/3
for i in range(1, N):
    eta[i] = eta[i-1] + (x[i]-x[i-1])*(((x[i-1]-eta[i-1])/x[i-1])*((xi[i-1]*x[i-1]*(x[i-1]-eta[i-1])-2)/((x[i-1]-eta[i-1])**2-1)))
    xi[i] = xi[i-1] + (x[i]-x[i-1])*(xi[i-1]*(((x[i-1]-eta[i-1])/x[i-1])*(xi[i-1]*x[i-1]-2*(x[i-1]-eta[i-1]))/((x[i-1]-eta[i-1])**2-1)))

etaHelper = scipy.interpolate.interp1d(x, xi)
def getEta(x):
    if x < 10.0**(-6):
        return 5.0/3
    if x > 10.0**2:
        return 8.86/x**2
    return float(etaHelper(x))

xiHelper = scipy.interpolate.interp1d(x, eta)

def getXi(x):
    if x < 1e-6:
        return (2*x)/3
    if x > 1e2:
        return 3.28
    return float(xiHelper(x))

if __name__ == "__main__":
    fig = plt.figure(figsize = (6, 8))
    ax1 = fig.add_subplot(111)
    l1, = ax1.plot(x, xi)
    ax1.set_xscale("log")
    ax1.set_xlim(10**-1, 10**2)
    ax1.set_yscale("log")
    ax1.set_ylim(10**-3, 10**1)
    ax1.text(10**0.25, 10**-0.5, "log $\eta$", fontsize=10)
    ax1.text(10**0.75, 10**0.15, "log $\u03BE$", fontsize=10)
    plt.plot(x, eta, color = 'purple')
    plt.plot(x, xi, color = 'blue')
    plt.show()

