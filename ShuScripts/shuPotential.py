'''
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
'''
import sys
import numpy as np
import os
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from scipy.optimize import brentq
import math
import rebound
import time
script_R = 3.36*10**3
T = 10
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.9891*10.0**30 ## kg
G = 6.67*10.0**(-11)  ## m^3*kg^-1*s^-2

def trap_integrate_log(f, xmin, xmax, N):
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

def trap_integrate_linear(f, xmin, xmax, N):
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

def years_to_seconds(years):
    return years * 365.25 * 24 * 3600

def get_b_max(potential, R, v0, t): #returns the maximum impact parameter for particles impacting the edge of the collapsing cloud at time t.
    phi = potential(R, t)
    return R*math.sqrt(1-2*phi/v0**2)

def get_params(potential, b, R, v0, t): # returns the components of the velocity vector given a potential function, and the radius R
    v = math.sqrt(v0**2 -2*potential(R, t) ) ##m/s
    vazimuthal = v0*b/R
    vr = math.sqrt(v**2 - vazimuthal**2)
    return [vr, vazimuthal]

def get_augmentation_factor_Shu(b, R, v_0, t):
    v = math.sqrt(v_0**2 -2*get_phi_Shu(R, t) ) ##m/s
    v_azimuthal = v_0*b/R
    v_r = math.sqrt(v**2 - v_azimuthal**2)
    return 1 + c_s/v_r


def getEnergy(particle, t):
    x = particle.x
    y = particle.y
    vx = particle.vx
    vy = particle.vy
    r = math.sqrt(x**2 + y**2)
    v = math.sqrt(vx**2 + vy**2)
    return 0.5*v**2 + get_phi_Shu(r, t)


def calculateFluxShu():
    return 10*10.0**(-12)*10.0**12/t_final

shu_x = np.logspace(-12, np.log10(2), 10000)

shu_integral = np.loadtxt('shuInt.txt')

shu_helper = scipy.interpolate.interp1d(shu_x, shu_integral, kind = 'cubic')
def get_shu_integral(x):
    if x < 10.0**(-12):
        return 0
    if x > 2:
        return shu_helper(2) + 2*(x-2)
    return shu_helper(x)

r_out = 1.974342*10.0**15 # radius in the Shu model with initial enclosed mass equal to one solar mass
##Is this Rmax?

def get_Shu_enclosed_mass(r, t):
    r0 = c_s*t
    centralMass = .975502*c_s**2*r0/G
    Mcalc =  centralMass + r0**3/(G*t**2)*get_shu_integral(r/r0)
    return min(Mcalc, M)

def get_Shu_little_g(r, t):
    return G*get_Shu_enclosed_mass(r, t)/r**2

def get_phi_Shu(r, t):
    if r > r_out:
        return -G*M/r
    return -G*M/(r_out) - trap_integrate_log(lambda x: get_Shu_little_g(x, t), r, r_out, 1000)

r_vals = np.logspace(13, np.log10(r_out), 10)

tsec = np.logspace(-1, np.log10(35000), 100)
t_values = np.array([years_to_seconds(t) for t in tsec])
phi_matrix = np.zeros((len(r_vals), len(t_values)))

for i, r in enumerate(r_vals):
    for j, t in enumerate(t_values):
        phi_matrix[i, j] = get_phi_Shu(r, t)

T, R = np.meshgrid(tsec, r_vals)
phi_flat = phi_matrix.flatten()
T_flat = T.flatten()
R_flat = R.flatten()

plt.figure(figsize=(10, 6))
sc = plt.scatter(T_flat, -phi_flat, c=R_flat, cmap='viridis', s=10, edgecolor='none')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('t (yr)', fontsize=14)
plt.ylabel(r'$-\phi_{\rm Shu}(r,t)$', fontsize=14)
cbar = plt.colorbar(sc, label='Radius (m)')
plt.title('Gravitational Potential Over Time with Radius as Color')
plt.tight_layout()
plt.show()
