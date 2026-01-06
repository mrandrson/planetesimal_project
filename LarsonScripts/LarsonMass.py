##Enclosed Mass
import scipy
#from scipy.integrate import cumtrapz as ct
import LarsonSolution as LS
import scipy.constants as sc
from scipy import integrate
import math as m
import numpy as np
import matplotlib.pyplot as plt
import os
from LarsonSolution import getEta, getXi
script_R = 3.36*10**3
T = 10
cs = np.sqrt(script_R*T)
N = 600

tf = 35000*3.154e7
tmax = tf/1e8
t = np.flip(-np.logspace(np.log10(tmax), np.log10(tf), N)) ##seconds

x = LS.x
trapIntegrateLinear = LS.trapIntegrateLinear
trapIntegrateLog = LS.trapIntegrateLog
r0 = -t*np.sqrt(3.36e4)/1.496e11 ## AU

r_values = [10, 100, 1000] # AU
colors = ['blue', 'red', 'green']

fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(111)
#ax1.set_xlim(-10**13, -10)
ax1.set_xscale("symlog")
ax1.set_yscale("log")
plt.xlabel('t(yr)', fontsize=25)
plt.ylabel('M$_{enc}$(kg)', fontsize=25)

for r in r_values:
    os.system(f"rm -rf LarsonMass{r}.txt")

for i, r in enumerate(r_values):
    M_enc = np.zeros(N)
    Mshu = np.zeros(N)
    for j in range(0, N):
        integral = scipy.integrate.quad(lambda x: x**2*LS.getEta(x), 0, r/r0[j])[0]
        M_enc[j] = ((-t[j]*(script_R*T)**(3/2))/sc.G)*integral ##kg
        x = r/(-cs*t[j])
        print(f"{100*(i*j)/(3*N):.2f}")
    np.savetxt(f"LarsonMass{r}.txt", np.column_stack((t, M_enc)), fmt="%.8e", delimiter=",", header="t,M", comments="") 
    ax1.plot(t/3.154e7, M_enc, label=f'M$_{{enc}}$ (r = {r} AU)', color=colors[i])


plt.legend(fontsize=20)
plt.show()

'''
def getMenc(r, t):
    r0 = -t*np.sqrt(3.36*10.0**3*10) * (31557600/149597870691) ## AU
    if r/r0 < 1:
        integral = trapIntegrateLinear(lambda x: x**2*LS.getEta(x), 0, r/r0, 100)
    else:
        integral = trapIntegrateLinear(lambda x: x**2*LS.getEta(x), 0, 1, 100) + trapIntegrateLog(lambda x: x**2*LS.getEta(x), 1, r/r0, 100) 
    return ((-t*(31557600)*(script_R*T)**(3/2))/sc.G)*integral

def getr_max(t):
    return scipy.optimize.brentq(lambda r: getMenc(r, t)-2*10.0**30, 0, 10000000)

def getPhi(r, t):
    rmax = getr_max(t)
    rmax_meters = rmax*14959787069
    phi_max = -sc.G*2*10.0**30/rmax_meters
    r_meters = r*14959787069
    return phi_max - sc.G*trapIntegrateLog(lambda x: getMenc(x/14959787069, t)/x**2, r_meters, rmax_meters, 100)

rvec = np.logspace(1, np.log10(4144), 30)
phivec = np.zeros(30)
for i in range(0, 30):
    print(i)
    phivec[i] = getPhi(rvec[i], -10000)


fig = plt.figure(figsize=(8, 8))
ax1 = fig.add_subplot(111)
ax1.plot(rvec, -phivec)
ax1.set_xscale("log")
plt.xlabel('r(AU)', fontsize=25)
plt.ylabel('$\\phi\\left(\\frac{m}{s}\\right)^2$', fontsize=25)
plt.show()
'''
 
