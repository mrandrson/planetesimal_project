import shuMass as sM
import scipy
import numpy as np
import scipy.constants as sc
from LarsonSolution import trapIntegrateLinear, trapIntegrateLog
import matplotlib.pyplot as plt

def r_max(t):
    return scipy.optimize.brentq(lambda r: sM.M_enc(r, t)-2*10**30, .4, 10000000)

def Phi(r, t):
    rmax = r_max(t)
    rmax_meters = rmax*149597870690
    phi_max = -sc.G*2*10.0**30/rmax_meters
    r_meters = r*149597870690
    return phi_max - sc.G*trapIntegrateLog(lambda x: (sM.Mp(t)+sM.M_enc(x/149597870690, t))/x**2, r_meters, rmax_meters, 100)

rvec = np.logspace(1, np.log10(15000), 30)
phivec = np.zeros(30)
for i in range(0, 30):
    print(i)
    phivec[i] = Phi(rvec[i], 10000)


fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(111)
ax1.plot(rvec, -phivec)
ax1.set_xscale("log")
plt.xlabel('r(AU)', fontsize=25)
plt.ylabel('$\\phi\\left(\\frac{m}{s}\\right)^2$', fontsize=25)
plt.show()
