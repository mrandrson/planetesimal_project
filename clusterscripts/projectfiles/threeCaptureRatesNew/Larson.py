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

bmin = float(sys.argv[1])
bmax = float(sys.argv[2])
v_0 = float(sys.argv[3])
t_init = -float(sys.argv[4])*24*3600*365.25
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

N = 100000  
x = np.zeros(N)
for i in range(0, N):
    x[i] = (1000*i)/(N-1) +.000001

eta = np.zeros(N)

eta[0] = 2/3*x[0]

xi = np.zeros(N)

xi[0] = 1.667
for i in range(1, N):
    eta[i] = eta[i-1] + (x[i]-x[i-1])*(((x[i-1]-eta[i-1])/x[i-1])*((xi[i-1]*x[i-1]*(x[i-1]-eta[i-1])-2)/((x[i-1]-eta[i-1])**2-1)))
    xi[i] = xi[i-1] + (x[i]-x[i-1])*(xi[i-1]*(((x[i-1]-eta[i-1])/x[i-1])*(xi[i-1]*x[i-1]-2*(x[i-1]-eta[i-1]))/((x[i-1]-eta[i-1])**2-1)))

eta_Helper = scipy.interpolate.interp1d(x, xi)
def get_Eta(x):
    if x < 10.0**(-6):
        return 1.667
    if x > 10.0**3:
        return 8.785135/x**2
    return float(eta_Helper(x))
#plt.show()




v_x = np.logspace(-6, 3, 10000)
'''
v_integral = np.zeros(10000)
int01 = trap_integrate_linear(lambda x: x**2*get_Eta(x), 0, 1, 10000)
for i in range(0, 10000):
    print(i)
    if v_x[i] < 1:
        v_integral[i] = trap_integrate_linear(lambda x: x**2*get_Eta(x), 0, v_x[i], 10000)
    else:
        v_integral[i] = int01 + trap_integrate_log(lambda x: x**2*get_Eta(x), 1, v_x[i], 10000)

np.savetxt('LarsonInt.txt', v_integral)
'''
v_integral = np.loadtxt('LarsonInt.txt')

Int_Helper = scipy.interpolate.interp1d(v_x, v_integral, kind = 'cubic')
def get_x_integral(x):
    if x < 10.0**(-6):
        return 5.0/9*x**3
    if x > 1000:
        return 8.785135*x-27.423076281380418
    return float(Int_Helper(x))
        

def get_Menc_Larson(r, t): #r in meters, t in seconds
    r0 = -c_s*t
    integral = get_x_integral(r/r0) 
    return (-t*c_s**3/G)*integral

def get_r_max(t): #t in seconds
    return scipy.optimize.brentq(lambda r: get_Menc_Larson(r, t)-1.9891*10.0**30, 0, 10.0**16)

def get_phi_Larson(r, t): 
    rmax = get_r_max(t)
    phi_max = -G*M/rmax
    if r > rmax:
        return -G*M/r
    return phi_max - G*trap_integrate_log(lambda x: get_Menc_Larson(x, t)/x**2, r, rmax, 10000)


def get_rdot_Larson(t):
    h = 10**(-5)*t 
    rdot = (get_r_max(t + h) - get_r_max(t - h)) / (2 * h)
    return rdot 


def get_b_max(potential, R, v0, t): #returns the maximum impact parameter for particles impacting the edge of the collapsing cloud at time t. 
    phi = potential(R, t)
    return R*math.sqrt(1-2*phi/v0**2) 

def get_t(r):
    return scipy.optimize.brentq(lambda t: get_r_max(t) - r, 0, 10000000) ##years


def get_params(potential, b, R, v0, t): # returns the components of the velocity vector given a potential function, and the radius R
    v = math.sqrt(v0**2 -2*potential(R, t) ) ##m/s
    vazimuthal = v0*b/R
    vr = math.sqrt(v**2 - vazimuthal**2)
    return [vr, vazimuthal]

def get_reduction_factor_Larson(b, R, v0, t): #R and t in 
    v = math.sqrt(v0**2 + 2*G*M/R) ##m/s
    psi = G*M/(R*v0**2) ##dimensionless
    theta = math.acos(b/(R*math.sqrt(1+2*psi))) ##dimensionless
    vr = -v*math.sin(theta) ##m/s
    Rdot = get_rdot_Larson(t)
    return max(0, (vr - Rdot)/vr) ##m/s


def getEnergy(particle, t):
    x = particle.x
    y = particle.y
    vx = particle.vx
    vy = particle.vy 
    r = math.sqrt(x**2 + y**2)
    v = math.sqrt(vx**2 + vy**2)
    return 0.5*v**2 + get_phi_Larson(r, t) 

def calculateFlux(n0, v0, t): 
    Rmax = get_r_max(t)
    vesc = math.sqrt(G*M/Rmax)
    return math.pi*Rmax**2*n0*v0*(1+(vesc/v0)**2)*get_reduction_factor_Larson(b, Rmax, v0, t)

def get_b(bmin, bmax):
    s = np.random.random()
    return math.sqrt(bmin**2 + s*(bmax**2-bmin**2))
 
def get_captured_number(bmin, bmax, n_0, v_0, t_init): #t_init in seconds
    R_max = get_r_max(t_init) ##meters
    Ncaptured = 0
    v_init = math.sqrt(v_0**2 + 2*G*M/R_max)
    wait_time = 2*R_max/v_init
    sim = rebound.Simulation()
    sim.integrator = "IAS15"
    sim.N_active = 0
#    sim.ri_ias15.min_dt = 10.0**6
    sim.G = 6.67*10.0**(-11)
    theta = 0
    sim.t = 10000
    def extraForce(reb_sim):
        if sim.t < wait_time:
            t = t_init
        else:
            t = sim.t - wait_time + t_init
        r0 = -c_s*t
        const = -t*c_s**3/G
        for i in range(0, len(sim.particles)):
            r = math.sqrt(sim.particles[i].x**2 + sim.particles[i].y**2) # in meters
            a = min(G*get_x_integral(r/r0)*const, G*M)/r**2 
            sim.particles[i].ax = -sim.particles[i].x/r*a
            sim.particles[i].ay = -sim.particles[i].y/r*a
        return
    sim.additional_forces = extraForce
    N = 0
    for i in range(0, 100):
        #print(N, i) 
        R_max = get_r_max(t_init)
        N_exp = -5/t_init*wait_time/100
        #print(N_exp, N, i, wait_time, t_init, R_max, v_0)
        #N_static = np.random.poisson(N_exp)
        r = np.random.random()
        if r < N_exp:
           N_static = 1
        else:
           N_static = 0
        bmiss = get_b_max(get_phi_Larson, R_max, v_0*1.025, t_init)
        for j in range(0, N_static):
            b = get_b(bmin, bmax)
            v = v_0*(1+.05*(np.random.random()-.5))
            if b > bmiss:
                pass
            else:
                theta = theta + 0.1
                [vxp, vyp] = get_params(get_phi_Larson, b, R_max, v, t_init)
                sim.add(m = 0, x = -R_max*math.cos(theta), y = R_max*math.sin(theta), vx = vxp*math.cos(theta) - vxp*math.sin(theta), vy = vyp*math.cos(theta)+ vxp*math.sin(theta))
                N = N + 1
        if N > 0:
            sim.integrate((i+1)*wait_time/100)
        else: 
            sim.t = (i+1)*wait_time/100
    tvec = -np.flip(np.logspace(np.log10(-t_init/10.0**6), np.log10(-t_init), 1000))
    for i in range(0, 999):
        #print(N, i)
        R_max = get_r_max(tvec[i])
        bmiss = get_b_max(get_phi_Larson, R_max, v_0*1.025, tvec[i])
        if bmiss > math.sqrt(bmin*bmax):
            N_exp = -5/t_init*(tvec[i+1] - tvec[i])*get_reduction_factor_Larson(math.sqrt(bmin*bmax), R_max, v_0, tvec[i])
        else:
            N_exp = 0
       # N_static = np.random.poisson(N_exp)
        r = np.random.random()
        if r < N_exp:
           N_static = 1
        else:
           N_static = 0
        for j in range(0, N_static):
            theta = theta + 0.1
            b = get_b(bmin, bmax)
            v = v_0*(1+.05*(np.random.random()-.5))
            if b < bmiss:
                [vxp, vyp] = get_params(get_phi_Larson, b, R_max, v, tvec[i])
                sim.add(m = 0, x = -R_max*math.cos(theta), y = R_max*math.sin(theta), vx = vxp*math.cos(theta) - vxp*math.sin(theta), vy = vyp*math.cos(theta)+ vxp*math.sin(theta))
                N = N + 1
        if N > 0:
            sim.integrate(wait_time - t_init + tvec[i+1])
        else:
            sim.t = wait_time - t_init + tvec[i+1]
    Evec = np.zeros(N)
    for i in range(0, N):
        Evec[i] = getEnergy(sim.particles[i], -10000)
    for i in range(0, len(Evec)):
        if Evec[i] < 0:
            Ncaptured = Ncaptured + 1
    return Ncaptured

s = 0
print('ready for business')
for i in range(0, reps):
    s = s + get_captured_number(bmin, bmax, 10**(-45), v_0, t_init)

np.savetxt('./data/Larson' + folder + '/'  +  'v' + str(v_index) + 'b' + b_index + '.txt', [s] )

