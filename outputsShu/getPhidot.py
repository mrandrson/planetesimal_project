import matplotlib.pyplot as plt
import numpy as np
import math 
import scipy

script_R = 3.36*10**3
T = 10 ##Kelvin
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.989*10.0**30 ## kg 
G = 6.674*10.0**(-11)  ## m^3*kg^-1*s^-2
a = 0.975471932310942752606143078377

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

x4, alpha4, v4 = loaddata("/Users/richardanderson/workdir/planetesimal_project/ShuScripts/shusolution.bin")

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
    seconds = years * 3.15576e7
    return seconds
def seconds_to_years(seconds):
    years = seconds / 3.15576e7
    return years
def meters3toau3(cubic_meters):
    cubic_au_in_meters = (1.49e11)**3
    return cubic_meters / cubic_au_in_meters


#rmax = meters_to_au((G*M)/(2*c_s**2)) 


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

v_integral = np.loadtxt('/Users/richardanderson/workdir/planetesimal_project/ShuScripts/shuInt.txt')

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

def phidot(r, t):
    h = t*1e-4
    phi1 = getPhi(r, t+h)
    phi2 = getPhi(r, t-h)
    phidotresult = (phi1-phi2)/(2*h)
    return phidotresult

def plotenergies(index):
    file_path = f"./particle_energies/particle_dat{index}.csv"
    try:
        energytime = pd.read_csv(file_path)
    except EmptyDataError:
        return
        
    for particle in energytime['particle']:
        plt.plot(energytime['time'][energytime['particle'] == particle], energytime['energy'][energytime['particle'] == particle])





