import numpy as np
import matplotlib.pyplot as plt
import scipy
from scipy import interpolate
from scipy.optimize import brentq
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
N_exp = int(sys.argv[2])
t_f = float(sys.argv[3])
filepath = sys.argv[4]

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

def takeStep(alpha, v, x, dx):
    vp = (alpha*(x-v)-2.0/x)*(x-v)/((x-v)**2-1)
    alphap = (alpha-2.0/x*(x-v))*(x-v)*alpha/((x-v)**2-1)
    vt = v + vp*dx
    alphat = alpha + alphap*dx
    xt = x + dx
    vp1 = (alphat*(xt-vt)-2.0/xt)*(xt-vt)/((xt-vt)**2-1)
    alphap1 = (alphat-2.0/xt*(xt-vt))*(xt-vt)*alphat/((xt-vt)**2-1)
    vp = 0.5*(vp + vp1)
    alphap = 0.5*(alphap + alphap1)
    return v + vp*dx, alpha + alphap*dx

def getSolution(A, Nsteps):
    xvec = np.logspace(np.log10(2), -12.0, Nsteps)
    vvec = np.zeros(Nsteps)
    alphavec = np.zeros(Nsteps)
    vvec[0] = (-(A-2)/xvec[0])-((1-A/6)*(A-2))/(xvec[0]**3)
    alphavec[0] = A/xvec[0]**2-A*(A-2)/(2*xvec[0]**4)
    for i in range(0, Nsteps-1):
        alpha = alphavec[i]
        x = xvec[i]
        v = vvec[i]
        dx = xvec[i+1]-xvec[i]
        vvec[i+1], alphavec[i+1] = takeStep(alpha, v, x, dx)
    return xvec, alphavec, vvec

x4, alpha4, v4 = getSolution(2.000001, 100000000)

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

rmax = meters_to_au((G*M)/(2*c_s**2))


v_x = np.logspace(-6, 4, 10000)
'''
int01 = trapIntegrateLinear(lambda x: x**2*getAlpha(x), 1*10**(-12), 1, 10000)
v_integral = np.zeros(10000)
for i in range(0, 10000):
    if v_x[i] < 1:
        v_integral[i] = trapIntegrateLinear(lambda x: x**2*getAlpha(x), 1*10**(-12), v_x[i], 10000)
    else:
        v_integral[i] = int01 + trapIntegrateLog(lambda x: x**2*getAlpha(x), 1, v_x[i], 10000)

np.savetxt('vIntShu.txt', v_integral)
'''

v_integral = np.loadtxt('vIntShu.txt')

IntHelper = scipy.interpolate.interp1d(v_x, v_integral, kind = 'cubic')

def get_x_integral(x):
    if x < 10.0**(-6):
        return .467*x**(3/2)
    if x > 100:
        return 2*x
    return IntHelper(x)

c_s = np.sqrt(3.36 * 10.0 ** 3 * 10) * (31557600 / 149597870691)  # AU/YR
c_sms =np.sqrt(3.36 * 10.0 ** 3 * 10) #m/s

def Mp(t):
    return (a * c_sms**3 * years_to_seconds(t))/G #kg

def getMenc(r, t):
    x = r/(c_s*t)
    integral = get_x_integral(x)
    return (years_to_seconds(t)*c_sms**3/G)*integral

def Mtot(r, t):
    return getMenc(r, t)+Mp(t)

def getr_max(t):
    return scipy.optimize.brentq(lambda r: (Mp(t)+getMenc(r, t)-2*10**30), .4, 20000)

def getPhi(r, t): #r in AU, t in years
    rmax_meters = au_to_meters(rmax)
    phi_max = -G*M/rmax_meters
    r_meters = au_to_meters(r)
    if r >= rmax:
        return -G*M/r_meters
    else:
        return phi_max - G*trapIntegrateLog(lambda rp: (Mtot(meters_to_au(rp), t))/rp**2, r_meters, rmax_meters, 10000)


def get_rdot(t):
    h = 10**(-4) ##years
    rdot = (getr_max(t + h) - getr_max(t - h)) / (2 * h)
    return rdot ##AU/years




def get_bmax(v0): #returns the maximum impact parameter for particles impacting the edge of the collapsing cloud at time t
    R = au_to_meters(rmax) ## - meters
    return R*np.sqrt(1+(4*c_sms**2)/(v0**2)) #meters

def get_t(r):
    return scipy.optimize.brentq(lambda t: getr_max(t) - r, 0, 10000000) ##years

#R in AU, b in meters, v0 in m/s
def get_params(b, v0, R, t): # returns the components of the velocity vector as the particle enters the mass distribution
    v = math.sqrt(v0**2+4*c_sms**2) ##m/s
    psi = G*M/(au_to_meters(R)*v0**2) ##dimensionless
    theta = math.acos(meters_to_au(b)/(R*math.sqrt(1+(4*c_sms**2/v0**2)))) ##dimensionless
    vr = v*math.sin(theta) ## m/s
    vazimuthal = v*math.cos(theta) ##m/s
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

def simulate_one_particle(v0, t_f):
    Rmax = au_to_meters(rmax)
    sim = rebound.Simulation()
    sim.integrator = "IAS15"
    #sim.ri_ias15.min_dt = 1e-6
    sim.G = G

    def extraForce(reb_sim):
        t = seconds_to_years(sim.t)
        for particle in sim.particles:
            r = math.sqrt(particle.x**2 + particle.y**2)
            if r >= Rmax:
                a = G*M/r**2
            else:
                a = G*(Mp(t)+getMenc(meters_to_au(r), t))/r**2
            particle.ax += -particle.x/r*a
            particle.ay += -particle.y/r*a
    sim.additional_forces = extraForce
    bmax = get_bmax(v0)
    b = np.sqrt(bmax**2 * np.random.random())
    [vr, vazimuthal] = get_params(b, v0, meters_to_au(Rmax), 0)
    sim.add(m=0, x=-Rmax, y=0, vx=vr, vy=vazimuthal)
    t_fs = years_to_seconds(t_f)
    for i in range(0, 100):
        #[sim.particles[0].vx, sim.particles[0].vy] = get_params(b, v0, np.sqrt(sim.particles[0].x**2+sim.particles[0].y**2), 0)
        x_arr = [sim.particles[0].x]
        y_arr = [sim.particles[0].y]
        vx_arr = [sim.particles[0].vx]
        vy_arr = [sim.particles[0].vy]
        ax_arr = [sim.particles[0].ax]
        ay_arr = [sim.particles[0].ay]
        np.savetxt(filepath+"/testacceleration/ax" + str(i)+'.txt', ax_arr)
        np.savetxt(filepath+"/testacceleration/ay" + str(i)+'.txt', ay_arr)
        np.savetxt(filepath+"/testvelocity/vx" + str(i)+'.txt', vx_arr)
        np.savetxt(filepath+"/testvelocity/vy" + str(i)+'.txt', vy_arr)
        np.savetxt(filepath+"/testposition/x" + str(i)+'.txt', x_arr)
        np.savetxt(filepath+"/testposition/y" + str(i)+'.txt', y_arr)
        energy = getEnergy(sim.particles[0], seconds_to_years(sim.t+0.001))
        np.savetxt(filepath+"/testtime/t"+str(i)+'.txt', [sim.t])
        np.savetxt(filepath+"/testenergy/E" + str(i) + '.txt', [energy])
        KE = .5*(sim.particles[0].vx**2+sim.particles[0].vy**2)
        PE = getPhi(meters_to_au(math.sqrt(sim.particles[0].x**2+sim.particles[0].y**2)), seconds_to_years(sim.t)+0.001)
        np.savetxt(filepath+"/testke/KE" + str(i) + '.txt', [KE])
        np.savetxt(filepath+"/testpe/PE" + str(i) + '.txt', [PE]) 
        sim.integrate(t_fs*(i+1)/100)
        if energy<0:
            print ("capture")
        else:
            print("no capture")

simulate_one_particle(v0, t_f)

