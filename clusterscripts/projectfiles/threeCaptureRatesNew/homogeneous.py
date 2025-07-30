import sys
import numpy as np
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
M = 2*10.0**30 ## kg 
G = 6.67*10.0**(-11)  ## m^3*kg^-1*s^-2
alpha = 0.2

bmin = float(sys.argv[1])
bmax = float(sys.argv[2])
v_0 = float(sys.argv[3])
R_init = float(sys.argv[4])
b_index = sys.argv[5]
v_index = sys.argv[6]
folder = sys.argv[7]

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







def get_b_max(R, v0): #returns the maximum impact parameter for particles impacting the edge of the collapsing cloud at time t. 
    return R*math.sqrt(1 + 2*G*M/(R*v0**2))

def get_r_max(t):
    return (-3*alpha*math.sqrt(G*M)*t/2)**(2.0/3)

def get_t(r):
    return scipy.optimize.brentq(lambda t: get_r_max(t) - r, 0, 10000000) ##years


def get_params(b, R, v0): # returns the components of the velocity vector given a potential function, and the radius R
    v = math.sqrt(v0**2 +2*G*M/R ) ##m/s
    vazimuthal = v0*b/R
    vr = math.sqrt(v**2 - vazimuthal**2)
    return [vr, vazimuthal]

def get_reduction_factor(b, R, v0):
    v = math.sqrt(v0**2 + 2*G*M/R) ##m/s
    psi = G*M/(R*v0**2) ##dimensionless
   # print(b, R, v0)
    theta = math.acos(b/(R*math.sqrt(1+2*psi))) ##dimensionless
    vr = v*math.sin(theta) ##m/s
    vrCloud = alpha*math.sqrt(G*M/R)
    return max(1 - vrCloud/vr, 0)
    
    

def getEnergy(particle, t):
    x = particle.x
    y = particle.y
    vx = particle.vx
    vy = particle.vy 
    r = math.sqrt(x**2 + y**2)
    v = math.sqrt(vx**2 + vy**2)
    return 0.5*v**2 -G*M/r 
    
def get_b(bmin, bmax):
    s = np.random.random()
    return math.sqrt(bmin**2 + s*(bmax**2-bmin**2))


def get_captured_number(bmin, bmax, v_0, R_init): #Nsteps is the number of individual simulations run, v0 the velocity at infinity, expected number the number of particles expected to be captured if Equation 12 is correct
    Ncaptured = 0
    t_init = -math.sqrt(4*R_init**3/(9*alpha**2*G*M))
    v_init = math.sqrt(v_0**2 + 2*G*M/R_init)
    wait_time = 2*R_init/v_init
    sim = rebound.Simulation()
    sim.integrator = "IAS15"
    sim.ri_ias15.min_dt = 10.0**4
    sim.G = 6.67*10.0**(-11)
    def extraForce(reb_sim):
        if sim.t < wait_time:
            t = t_init
        else: 
            t = t_init + (sim.t - wait_time)
        for i in range(0, len(sim.particles)):
            r = math.sqrt(sim.particles[i].x**2 + sim.particles[i].y**2) # in meters
            r_edge = (-3*alpha*math.sqrt(G*M)*t/2)**(2.0/3)
            if r > r_edge:
                a = G*M/r**2 #regular force outside of sphere of mass M
            else:
                a = G*M*r/r_edge**3 
            sim.particles[i].ax = -sim.particles[i].x/r*a
            sim.particles[i].ay = -sim.particles[i].y/r*a
            #if sim.dt < 0.1:
                #print(sim.dt, r/r_edge, r_edge, r, a, t, N)
        return
    sim.additional_forces = extraForce
    N = 0
    for i in range(0, 100):
        N_exp = -1/t_init*wait_time/100
        N_static = np.random.poisson(N_exp)
        bmiss = get_b_max(R_init, v_0*1.025)
        for j in range(0, N_static):
            b = get_b(bmin, bmax)
            if b > bmiss:
                pass
            else:
                [vxp, vyp] = get_params(b, R_init, v_0)
                sim.add(m = 0, x = -R_init, y = 0, vx = vxp, vy = vyp)
                N = N + 1
        if N > 0:
            sim.integrate((i+1)*wait_time/100)
        else: 
            sim.t = (i+1)*wait_time/100
    tvec = -np.flip(np.logspace(np.log10(-t_init/10.0**9), np.log10(-t_init), 10000))
    for i in range(0, 9999):
        R_max = get_r_max(tvec[i])
        bmiss = get_b_max(R_max, v_0*1.025)
        if bmiss > math.sqrt(bmin*bmax):
            N_exp = -1/t_init*(tvec[i+1] - tvec[i])*get_reduction_factor(math.sqrt(bmin*bmax), R_max, v_0)
        else:
            N_exp = 0
        N_static = np.random.poisson(N_exp)
        for j in range(0, N_static):
            b = get_b(bmin, bmax)
            if b < bmiss:
                [vxp, vyp] = get_params(b, R_max, v_0)
                sim.add(m = 0, x = -R_max, y = 0, vx = vxp, vy = vyp)
                N = N + 1
        if N > 0:
            sim.integrate(wait_time - t_init + tvec[i+1])
        else:
            sim.t = wait_time - t_init + tvec[i+1]
    Evec = np.zeros(N)
    for i in range(0, N):
        Evec[i] = getEnergy(sim.particles[i], -0.0000001)
    for i in range(0, len(Evec)):
        if Evec[i] < 0:
            Ncaptured = Ncaptured + 1
    return Ncaptured

s = 0
for i in range(0, 300):
    s = s + get_captured_number(bmin, bmax, v_0, R_init)

np.savetxt('./data/homogeneous' + folder + '/'  +  'v' + str(v_index) + 'b' + str(b_index)+ '.txt', [s] )

