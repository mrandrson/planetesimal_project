import math
import numpy as np
import matplotlib.pyplot as plt
import scipy

script_R = 3.36*10**3
T = 10 ##Kelvin
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.989*10.0**30 ## kg 
G = 6.674*10.0**(-11)  ## m^3*kg^-1*s^-2
nsamples = 249
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

##Positions in meters##
x = [np.loadtxt(f'./testposition/x{i}.txt') for i in range(nsamples)]
x = np.array(x)
y = [np.loadtxt(f'./testposition/y{i}.txt') for i in range(nsamples)]
y = np.array(y)
r = [np.sqrt(x[i]**2+y[i]**2) for i in range(nsamples)]
r = np.array(r)

##Velocities in meters per second##
vx = [np.loadtxt(f'./testvelocity/vx{i}.txt') for i in range(nsamples)]
vx = np.array(vx)
vy = [np.loadtxt(f'./testvelocity/vy{i}.txt') for i in range(nsamples)]
vy = np.array(vy)
v = [np.sqrt(vx[i]**2+vy[i]**2) for i in range(nsamples)]
v = np.array(v)

##Simulation time in seconds##
t = [np.loadtxt(f'./testtime/t{i}.txt') for i in range(nsamples)]
t = np.array(t)

##Simulation Energy##
Ecluster = [np.loadtxt(f'./testenergy/E{i}.txt') for i in range(nsamples)]
Ecluster = np.array(Ecluster)

##rmax in meters, then converted to AU##
rmax = meters_to_au((G*M)/(2*c_s**2))

intvec = np.logspace(-6, 4, 10000)
v_integral = np.loadtxt('./vIntShu.txt')

IntHelper = scipy.interpolate.interp1d(intvec, v_integral, kind = 'cubic')

def get_x_integral(x):
    if x < 10.0**(-6):
        return .467*x**(3/2)
    if x > 100:
        return 2*x
    return IntHelper(x)

c_s = np.sqrt(3.36 * 10.0 ** 3 * 10) * (31557600 / 149597870691)  # AU/YR
c_sms =np.sqrt(3.36 * 10.0 ** 3 * 10) #m/s

def Mp(t):
    return (.98 * c_sms**3 * years_to_seconds(t))/G #kg

def getMenc(r, t):
    x = r/(c_s*t)
    integral = get_x_integral(x)
    return (years_to_seconds(t)*c_sms**3/G)*integral

def Mtot(r, t):
    return getMenc(r, t)+Mp(t)

def getr_max(t):
    return scipy.optimize.brentq(lambda r: (Mp(t)+getMenc(r, t)-M), .4, 20000)

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

def getPhi(r, t): #r in AU, t in years
    rmax_meters = au_to_meters(rmax)
    phi_max = -G*M/rmax_meters
    r_meters = au_to_meters(r)
    if r >= rmax:
        return -G*M/r_meters
    else:
        return phi_max - G*trapIntegrateLog(lambda rp: (Mtot(meters_to_au(rp), t))/rp**2, r_meters, rmax_meters, 10000)

##Energy Arrays##
PE = [getPhi(meters_to_au(r[i]), seconds_to_years(t[i])+0.001) for i in range(nsamples)]
PE = np.array(PE)
KE = [(1/2)*(vx[i]**2+vy[i]**2) for i in range(nsamples)]
KE = np.array(KE)
Elocal = PE+KE

min_index = np.argmin(Ecluster)
min_value = Ecluster[min_index]
#plt.yscale("symlog")
plt.plot(t, KE, label='$E_k$')
plt.plot(t, Ecluster, label='$E_{Cluster}$')
plt.plot(t, PE, label = '$E_p$')
#plt.plot(t, Elocal, label = '$E_{Local}$')
#plt.plot(t[min_index], min_value, 'ro', markersize=6)
#plt.plot(t[min_index+1], Ecluster[min_index+1], 'ro', markersize=6)
plt.text(8e12, -45000, f'$\\Delta E =\\ ${Ecluster[248]-Ecluster[0]:.2e} $\\frac{{m^2}}{{s^2}}$', fontsize=12, verticalalignment='top', horizontalalignment='right', bbox=dict(facecolor='white', alpha=0.5))

plt.xlabel('$t \\ (s)$')
plt.ylabel('Energy $\\left(\\frac{m^2}{s^2}\\right)$')
plt.text(4e9, 2e4, '$v_0 = 10\\ \\frac{{m}}{{s}} \\quad N_{exp} = 20 \\quad t_f = 500,000\\ yr$')
plt.legend()

plt.show()

