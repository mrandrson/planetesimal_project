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

'''
bmin = float(sys.argv[1])
bmax = float(sys.argv[2])
v_0 = float(sys.argv[3])
t_final = float(sys.argv[4])
b_index = sys.argv[5]
v_index = sys.argv[6]
folder = sys.argv[7]
reps = int(sys.argv[8])
'''

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

def get_b_max(potential, R, v0, t): #returns the maximum impact parameter for particles impacting the edge of the collapsing cloud at time t.
    phi = potential(R, t)
    return R*math.sqrt(1-2*phi/v0**2)

'''
def get_t(r):
    return scipy.optimize.brentq(lambda t: get_r_max(t) - r, 0, 10000000) ##years
'''

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

r_vals = np.linspace(1e-12, r_out, 100)

import matplotlib.animation as animation

fig, ax = plt.subplots(figsize=(8, 5))
line_phi, = ax.plot([], [], label="Gravitational Potential $\phi_{Shu}(r,t)$", color='blue')
ax.set_xlim(r_vals.min(), r_vals.max())
ax.set_xlabel("Radius r (m)")
ax.set_ylabel("Gravitational Potential $\phi_{Shu}(r,t)$ (J/kg)")
ax.set_title("Shu Model Gravitational Potential Evolution Over Time")
ax.legend()
ax.grid()

t_values = np.linspace(1, 350000*3.154e7, 50)

def update(frame):
    t = t_values[frame]
    phi_vals = np.array([get_phi_Shu(r, t) for r in r_vals])
    line_phi.set_data(r_vals, phi_vals)
    y_min, y_max = phi_vals.min(), phi_vals.max()
    ax.set_ylim(y_min * 1.1, y_max * 0.9)
    ax.set_title(f"Shu Model Gravitational Potential at t = {t:.2e} s")
    return line_phi,

ani = animation.FuncAnimation(fig, update, frames=len(t_values), interval=50, blit=False)

plt.show()

