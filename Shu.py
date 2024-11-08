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

bmin = float(sys.argv[1])
bmax = float(sys.argv[2])
v_0 = float(sys.argv[3])
t_final = float(sys.argv[4])
b_index = sys.argv[5]
v_index = sys.argv[6]
folder = sys.argv[7]
reps = int(sys.argv[8])

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

##Where does this come from? 
 



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
    
    
def get_b(bmin, bmax):
    s = np.random.random()
    return math.sqrt(bmin**2 + s*(bmax**2-bmin**2))

def get_captured_number(bmin, bmax, v_0, t_final): #Nsteps is the number of individual simulations run, v0 the velocity at infinity, expected number the number of particles expected to be captured if Equation 12 is correct
    Ncaptured = 0
    sim = rebound.Simulation()
    sim.integrator = "IAS15"
   # sim.ri_ias15.min_dt = 10.0**6
    sim.G = 6.67*10.0**(-11)
    def extraForce(reb_sim):
        t = sim.t
        rad = c_s*t
        for i in range(0, len(sim.particles)):
            r = math.sqrt(sim.particles[i].x**2 + sim.particles[i].y**2) # in meters
            M_enc = get_Shu_enclosed_mass(r, t)
            a = G*M_enc/r**2
            sim.particles[i].ax = -sim.particles[i].x/r*a
            sim.particles[i].ay = -sim.particles[i].y/r*a
        return
    sim.additional_forces = extraForce
    N = 0
    Nsteps = 100
    tvec = np.logspace(np.log10(t_final/10.0**6), np.log10(t_final), Nsteps+1)
    sim.t = tvec[0]
    for i in range(0, Nsteps):
        R_max = c_s*tvec[i]
        ##Why the extra factor for v0
        bmiss = get_b_max(get_phi_Shu, R_max, 1.025*v_0, tvec[i])
        if bmiss > math.sqrt(bmin*bmax):
            N_exp = 1/t_final*(tvec[i+1] - tvec[i])*get_augmentation_factor_Shu(math.sqrt(bmin*bmax), R_max, v_0, tvec[i])
        else:
            N_exp = 0
        N_static = np.random.poisson(N_exp)
        for j in range(0, N_static):
            b = get_b(bmin, bmax)
            v = v_0*(1+.05*(np.random.random()-.5))
            if b < bmiss:
                [vxp, vyp] = get_params(get_phi_Shu, b, R_max, v, tvec[i])
                sim.add(m = 0, x = -R_max, y = 0, vx = vxp, vy = vyp)
                N = N + 1
                print(f'There are {N} particles in the simulation now')
                print('(', 6.685e-12*sim.particles[j].x, 6.685e-12*sim.particles[j].y, ')')
                print('dt:', sim.dt)
        if N > 0:
            sim.integrate(tvec[i+1])
        else:
            sim.t = tvec[i+1]
    Evec = np.zeros(N)
    for i in range(0, N):
        Evec[i] = getEnergy(sim.particles[i], tvec[Nsteps-1])
    for i in range(0, len(Evec)):
        if Evec[i] < 0:
            Ncaptured = Ncaptured + 1
    return Ncaptured
for i in range(10):
    print(get_captured_number(bmin, bmax, v_0, t_final))

'''
s = 0
for i in range(0, reps):
    s = s + get_captured_number(bmin, bmax, v_0, t_final)
    print('capture number:', s)
directory = './data/Shu' + folder
os.makedirs(directory, exist_ok=True)

np.savetxt('./data/Shu' + folder + '/'  +  'v' + v_index + 'b' + b_index+ '.txt', [s] )
'''
