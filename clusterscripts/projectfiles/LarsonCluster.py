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

v0 = float(sys.argv[1])
N_exp = int(sys.argv[2])
t_init = float(sys.argv[3])
filePath = sys.argv[4]

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
fig = plt.figure(figsize = (8, 10))
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

etaHelper = scipy.interpolate.interp1d(x, xi)
def getEta(x):
    if x < 10.0**(-6):
        return 5.0/3
    if x > 10.0**2:
        return 8.86/x**2
    return float(etaHelper(x))
#plt.show()

def au_to_meters(au):
    meters = au * 149597870700  
    return meters
def meters_to_au(meters):
    au = meters / 149597870700
    return au
def years_to_seconds(years):
    seconds = years * 31557600 
    return seconds
def seconds_to_years(seconds):
    years = seconds / 31557600
    return years




v_x = np.logspace(-6, 4, 10000)
'''
v_integral = np.zeros(10000)
int01 = trapIntegrateLinear(lambda x: x**2*getEta(x), 0, 1, 10000)
for i in range(0, 10000):
    print(i)
    if v_x[i] < 1:
        v_integral[i] = trapIntegrateLinear(lambda x: x**2*getEta(x), 0, v_x[i], 10000)
    else:
        v_integral[i] = int01 + trapIntegrateLog(lambda x: x**2*getEta(x), 1, v_x[i], 10000)

np.savetxt('vInt.txt', v_integral)
'''
v_integral = np.loadtxt('vInt.txt')

IntHelper = scipy.interpolate.interp1d(v_x, v_integral, kind = 'cubic')
def get_x_integral(x):
    if x < 10.0**(-6):
        return 5.0/9*x**3
    if x > 100:
        return 8.86*x-35.252838
    return IntHelper(x)
        

def getMenc(r, t):
    r0 = meters_to_au(c_s*years_to_seconds(-t)) ## AU
    integral = get_x_integral(r/r0) 
    return (-years_to_seconds(t)*c_s**3/G)*integral

def getr_max(t):
    return scipy.optimize.brentq(lambda r: getMenc(r, t)-2*10.0**30, 0, 10000000)

def getPhi(r, t): #r in AU, t in years
    rmax = getr_max(t)
    rmax_meters = au_to_meters(rmax)
    phi_max = -G*M/rmax_meters
    r_meters = au_to_meters(r)
    if r > rmax:
        return -G*M/r_meters
    return phi_max - G*trapIntegrateLog(lambda x: getMenc(meters_to_au(x), t)/x**2, r_meters, rmax_meters, 10000)


def get_rdot(t):
    h = 10**(-4) ##years
    rdot = (getr_max(t + h) - getr_max(t - h)) / (2 * h)
    return rdot ##AU/years




def get_bmax(v0, t): #returns the maximum impact parameter for particles impacting the edge of the collapsing cloud at time t
    R = au_to_meters(getr_max(t)) ## - meters
    return R*math.sqrt(1+2*G*M/(R*v0**2)) #meters

def get_t(r):
    return scipy.optimize.brentq(lambda t: getr_max(t) - r, 0, 10000000) ##years


def get_params(b, v0, R): # returns the components of the velocity vector as the particle enters the mass distribution
    v = math.sqrt(v0**2 + 2*G*M/R) ##m/s
    psi = G*M/(R*v0**2) ##dimensionless
    theta = math.acos(b/(R*math.sqrt(1+2*psi))) ##dimensionless
    vr = v*math.sin(theta) ## m/s
    vazimuthal = v*math.cos(theta) ##m/s
    return [vr, vazimuthal]

def get_reduction_factor(b, v0, R, t): #returns the fraction by which the flux is reduced due to the contraction of the sphere relative to the case of a static sphere.  Please check this!
    v = math.sqrt(v0**2 + 2*G*M/R) ##m/s
    psi = G*M/(R*v0**2) ##dimensionless
    theta = math.acos(b/(R*math.sqrt(1+2*psi))) ##dimensionless
    vr = v*math.sin(theta) ##m/s
    Rdot = get_rdot(t)*(au_to_meters(1)/seconds_to_years(1)) ##m/s
    return max(0, (vr - Rdot)/vr) ##m/s

def getEnergy(particle, t):
    x = particle.x
    y = particle.y
    vx = particle.vx
    vy = particle.vy 
    r = math.sqrt(x**2 + y**2)
    v = math.sqrt(vx**2 + vy**2)
    return 0.5*v**2 + getPhi(meters_to_au(r), t) ##(m/s)^2

def calculateFlux(n0, v0, t_0): #This gives the expected flux of particles onto the sphere in terms of their number density and velocity at infinity.  Make sure you understand where this comes from
    b = get_bmax(v0, t_0)##meters
    return math.pi*b**2*n0*v0 ###particles/s
 
 
def get_n0(expectedNumber, v0, Rmax): #gives the number density of planetesimals given the number expected to be captured, via inversion of equation 12
    if v0 > 5*c_s:
        return get_n0(expectedNumber, 5*c_s, Rmax)
    f = v0/c_s
    return 3*B*c_s*expectedNumber/(math.pi*v0*Rmax**3*(1 + 4*B/(3*f**2))) ###m^-3

###New
def calculate_acceleration(r, t):
    h = 10**(-4)
    a = -((getPhi(r+h, t)-getPhi(r-h, t))/(2*h))
    return meters_to_au(a)
###



def get_captured_number(v0, expectedNumber, t_init): #Nsteps is the number of individual simulations run, v0 the velocity at infinity, expected number the number of particles expected to be captured if Equation 12 is correct
    Rmax = au_to_meters(getr_max(t_init)) ##meters
    n0 = get_n0(expectedNumber, v0, Rmax) ##m^-3
    Ncaptured = 0
    v_init = math.sqrt(v0**2 + 2*G*M/Rmax)
    wait_time = 2*Rmax/v_init
    sim = rebound.Simulation()
    sim.integrator = "IAS15"
    sim.ri_ias15.min_dt = 10.0**6
    sim.G = 6.67*10.0**(-11)  #m^3kg^-1s^-2
    def extraForce(reb_sim):
        if sim.t < wait_time:
            t = t_init
        else:
            t = seconds_to_years(sim.t - wait_time) + t_init
        rad = au_to_meters(getr_max(t)) #rad in AU
        r0 = max(.001, meters_to_au(years_to_seconds(-t)*c_s)) #in AU
        const = ((years_to_seconds(-t)*c_s**3)/G)
        for i in range(0, len(sim.particles)):
            r = math.sqrt(sim.particles[i].x**2 + sim.particles[i].y**2) # in meters
            if r > rad:
                a = G*M/r**2 #regular force outside of sphere of mass M
            else:
                a = G*get_x_integral(meters_to_au(r)/r0)*const/r**2 
            sim.particles[i].ax = -sim.particles[i].x/r*a
            sim.particles[i].ay = -sim.particles[i].y/r*a
        return
    sim.additional_forces = extraForce
    N = 0
    bmax = get_bmax(v0, t_init)
    for i in range(0, 100):
        N_static = np.random.poisson(calculateFlux(n0, v0, t_init)*wait_time/100)
        for j in range(0, N_static):
            b = math.sqrt(bmax**2*np.random.random()) # Note that the impact parameter b is distributed such that b^2 is uniform; try to derive this.
            [vxp, vyp] = get_params(b, v0, Rmax)
            sim.add(m = 0, x = -Rmax, y = 0, vx = vxp, vy = vyp)
            N = N + 1
        if N > 0:
            sim.integrate((i+1)*wait_time/100)
        else: 
            sim.t = (i+1)*wait_time/100
    for i in range(0, 100):
        N_static = np.random.poisson(calculateFlux(n0, v0, t_init)*(-t_init/100))
        bmax = get_bmax(v0, t_init*(1-i/100))
        xvec = np.zeros(len(sim.particles))
        yvec = np.zeros(len(sim.particles))
        for k in range(0, len(sim.particles)):
            xvec[k] = sim.particles[k].x
            yvec[k] = sim.particles[k].y
        np.savetxt('./data/positions/x' + str(i+1) + '.txt', xvec)
        np.savetxt('./data/positions/y' + str(i+1) + '.txt', yvec)
        for j in range(0, N_static):
            b = math.sqrt(bmax**2*np.random.random()) # Note that the impact parameter b is distributed such that b^2 is uniform; try to derive this.
            Rmax = au_to_meters(getr_max(t_init*(1-i/100)))
            [vxp, vyp] = get_params(b, v0, Rmax)
            r = get_reduction_factor(b, v0, Rmax, t_init*(1-i/100))
            s = np.random.random()
            if s < r:
                sim.add(m = 0, x = -Rmax, y = 0, vx = vxp, vy = vyp)
                N = N + 1
        if N > 0:
            sim.integrate(wait_time - years_to_seconds(t_init*(i+1)/100))
        else:
            sim.t = wait_time - years_to_seconds(t_init*(i+1)/100)
    Evec = np.zeros(N)
    for i in range(0, N):
        Evec[i] = getEnergy(sim.particles[i], -0.0000001)
    for i in range(0, len(Evec)):
        if Evec[i] < 0:
            Ncaptured = Ncaptured + 1
    f = v0/c_s
    Nanalytic = math.pi*n0*v0*Rmax**3/(3*B*c_s)*(1+4*B/(3*f**2)) ##m^3. 
    return Ncaptured/Nanalytic


s = get_captured_number(v0, N_exp, t_init)

np.savetxt(filePath, [s])
