import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from scipy.optimize import brentq
import h5py
import math
import rebound
import sys
import time
script_R = 3.36*10**3
T = 10 ##Kelvin
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.989*10.0**30 ## kg 
G = 6.674*10.0**(-11)  ## m^3*kg^-1*s^-2
a = 0.975471932310942752606143078377

v0 = float(sys.argv[1])
#N_exp = int(sys.argv[2])
t_f = float(sys.argv[2])
#filepath = sys.argv[4]

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

def loaddata(filename):
    with open(filename, "rb") as file:
        Nsteps = np.fromfile(file, dtype=np.int32, count=1)[0]

        xvec = np.fromfile(file, dtype=np.float64, count=Nsteps)
        alphavec = np.fromfile(file, dtype=np.float64, count=Nsteps)
        vvec = np.fromfile(file, dtype=np.float64, count=Nsteps)

    return xvec, alphavec, vvec

x4, alpha4, v4 = loaddata("shusolution.bin")

alphaHelp = scipy.interpolate.interp1d(x4, alpha4)
vHelp = scipy.interpolate.interp1d(x4, v4)
 
def getAlpha(x): 
    if x>1: 
        return 2*x**(-2) 
    else:  
        return float(alphaHelp(x)) 
 
def getV(x): 
    return float(vHelp(x)) 
 
def au_to_meters(au):
    au_in_meters = 1.49e11 # meters
    return au * au_in_meters
def meters_to_au(au):
    au_in_meters = 1.49e11 # meters
    return au / au_in_meters
def years_to_seconds(years): 
    seconds = years * 31557600 
    return seconds 
def seconds_to_years(seconds): 
    years = seconds / 31557600 
    return years 
 
v_x = np.logspace(-6, 4, 10000) 

v_integral = np.loadtxt('shuInt.txt') 
 
IntHelper = scipy.interpolate.interp1d(v_x, v_integral, kind = 'cubic')
 
def get_x_integral(x):
    if x < 10.0**(-6):
        return 0
    if x > 2:
        return IntHelper(2) + 2*(x-2)
    return IntHelper(x)

def Mtot(r, t):
    r0 = c_s*t
    centralMass = .975502*c_s**2*r0/G
    Mcalc =  centralMass + r0**3/(G*t**2)*get_x_integral(r/r0)
    return min(Mcalc, M)

def getr_max(t): 
    return scipy.optimize.brentq(lambda r: (Mp(t)+getMenc(r, t)-2*10**30), .4, 20000)
 
def getPhi(r, t):  
    rmax = (G*M)/(2*c_s**2) 
    phi_max = -G*M/rmax 
    r_meters = r 
    if r >= rmax: 
        return -G*M/r_meters 
    else: 
        return phi_max - G*trapIntegrateLog(lambda rp: (Mtot(rp, t))/rp**2, r_meters, rmax, 10000)
 
 
def get_bmax(R, v0, t): 
    phi = getPhi(R, t)
    return R*math.sqrt(1-2*phi/v0**2) 
 
def get_params(b, R, v0, t):
    potential = getPhi
    v = math.sqrt(v0**2 -2*potential(R, t) ) ##m/s
    vazimuthal = v0*b/R
    vr = math.sqrt(v**2 - vazimuthal**2)
    return [vr, vazimuthal] 
 
def getEnergy(particle, t): 
    x = particle.x 
    y = particle.y 
    vx = particle.vx 
    vy = particle.vy 
    r = math.sqrt(x**2 + y**2) 
    v = math.sqrt(vx**2 + vy**2) 
    return 0.5*v**2 + getPhi(meters_to_au(r), t) ##(m/s)^2
 
def calculateFlux(n0, v0, t): #This gives the expected flux of particles onto the sphere in terms of their number density and velocity at infinity.  Make sure you understand where this comes from
    b = get_bmax(v0, t)##meters 
    return math.pi*b**2*n0*v0 ###particles/s
 
 
def get_n0(expectedNumber, v0, Rmax): #gives the number density of planetesimals given the number expected to be captured, via inversion of equation 12
    if v0 > 5*c_s: 
        return get_n0(expectedNumber, 5*c_s, Rmax)
    f = v0/c_s 
    return (9*expectedNumber*f)/(np.pi*Rmax**3*(3*f**2+16)) ###m^-3
 
def calculate_acceleration(r, t): 
    h = 10**(-4) 
    a = -((getPhi(r+h, t)-getPhi(r-h, t))/(2*h))
    return meters_to_au(a) 
 
def simulate_n_particles(v0, t_f, n): 
    Ncaptured = 0 
    Rmax = (G*M)/(2*c_s**2) 
    sim = rebound.Simulation()  
    sim.integrator = "IAS15"  
    #sim.ri_ias15.min_dt = 1e-2  
    #sim.save_to_file("simdata.bin", step = 1000) 
    sim.G = G  
  
    def extraForce(reb_sim):  
        t = sim.t
        for particle in sim.particles: 
            for i in range(0, len(sim.particles)):
                r = math.sqrt(sim.particles[i].x**2 + sim.particles[i].y**2) # in meters
                M_enc = Mtot(r, t) 
                a = G*M_enc/r**2
                dt = sim.dt

            particle.ax = -particle.x/r*a
            particle.ay = -particle.y/r*a

    sim.additional_forces = extraForce  
    N = 0   
    Nsteps = 250   
    t_final = years_to_seconds(t_f)   
    tvec = np.logspace(np.log10(t_final/10.0**5), np.log10(t_final), Nsteps+1)
    sim.t = tvec[0]   
    bmax = get_bmax(Rmax, 1.025*v0, sim.t)   
    b = np.sqrt(bmax**2 * np.random.random())
    [vr, vazimuthal] = get_params(b, Rmax, v0,sim.t)
    sim.add(m=0, x= -Rmax, y = 0, vx = vr, vy = vazimuthal)
    x_vals = []
    y_vals = []
    energy_vals = []
    
    for i in range(0, Nsteps):
        sim.integrate(tvec[i+1])
        particle = sim.particles[0]
        x_vals.append(particle.x)
        y_vals.append(particle.y)

        energy = getEnergy(particle, tvec[i])
        energy_vals.append(energy)

    Evec = np.array(energy_vals)

    for i in range(0, len(Evec)):
        if Evec[i] < 0:
            Ncaptured = 1

    plt.figure(figsize=(8, 6))
    plt.plot(x_vals, y_vals, label="Particle Motion")
    plt.xlabel("x (AU)")
    plt.ylabel("y (AU)")
    circle = plt.Circle((0, 0), Rmax, color='blue', fill=False, linestyle='--', label=f'Rmax = {Rmax:.2e} m')
    plt.gca().add_patch(circle) 
    plt.legend()
    plt.show()

    plt.figure(figsize=(8, 6))
    plt.plot(tvec[1:], Evec, label="Energy vs Time", color='red')
    plt.xlabel("Time (s)")
    plt.ylabel("Energy $\\left[kg\\cdot \\dfrac{{m^2}}{{s^2}}\\right]$")
    plt.yscale('symlog')
    plt.legend()
    plt.show()

    return Ncaptured

out = simulate_n_particles(v0, t_f, 1)   
print(out)
              
