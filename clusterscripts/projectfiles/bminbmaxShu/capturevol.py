import pandas as pd
import scipy
import os
import numpy as np
import re
import math
script_R = 3.36*10**3
T = 10
B = 8.86
M = 1.9891*10.0**30 ## kg 
G = 6.67*10.0**(-11)  ## m^3*kg^-1*s^-2
cs = math.sqrt(3.36)*100
bvec = np.logspace(9, 17, 80)
Nv = 30
vvec = np.logspace(np.log10(cs/10), np.log10(cs*5), Nv)
tf = 1.0829*10.0**13 #s

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

def get_t(r):
    return scipy.optimize.brentq(lambda t: get_r_max(t) - r, 0, 10000000) ##years


def get_params(potential, b, R, v0, t): # returns the components of the velocity vector given a potential function, and the radius R
    v = math.sqrt(v0**2 -2*potential(R, t) ) ##m/s
    vazimuthal = v0*b/R
    vr = math.sqrt(v**2 - vazimuthal**2)
    return [vr, vazimuthal]


def get_augmentation_factor_Shu(b, R, v_0, t):
    v = math.sqrt(v_0**2 -2*get_phi_Shu(R, t) ) ##m/s
    v_azimuthal = v_0*b/R
    v_r = math.sqrt(v**2 - v_azimuthal**2)
    return 1 + cs/v_r


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

shu_integral = np.loadtxt('/home/rbanderson/projectfiles/bminbmaxShu/shuInt.txt')

shu_helper = scipy.interpolate.interp1d(shu_x, shu_integral, kind = 'cubic')
def get_shu_integral(x):
    if x < 10.0**(-12):
        return 0
    if x > 2:
        return shu_helper(2) + 2*(x-2)
    return shu_helper(x)

r_out = 1.974342*10.0**15 # radius in the Shu model with initial enclosed mass equal to one solar mass

def get_Shu_enclosed_mass(r, t):
    r0 = cs*t
    centralMass = .9754516*cs**2*r0/G
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

V = {}

bmax = 0

for v in vvec:
    V[f'{v}'] = [0]

for f in os.listdir('./data/'):
    if f.endswith('.csv'):
        dat = pd.read_csv('./data/'+f)
        pattern = r"v(\d+)b(\d+)"

        match = re.search(pattern, f)

        if match:
            v_index = int(match.group(1))
            b_index = int(match.group(2))
        else:
            print("No match found")

        Ncapture = 0

        for p in dat['particle_id'].unique():
            pdat = dat[dat['particle_id'] == p]
            
            if any(dat['energy']<0):
                if bvec[b_index+1] > bmax:
                    bmax = bvec[b_index+1]
                    print('b_max index:', b_index)

                Ncapture += 1

        #Vcapture = np.pi*Ncapture*(bvec[b_index+1]**2-bvec[b_index]**2)*vvec[v_index]*tf*get_augmentation_factor_Shu(math.sqrt(bvec[b_index+1]*bvec[b_index]), r_out, vvec[v_index], r_out/cs)
        Vcapture = np.pi*Ncapture*(bvec[b_index+1]**2-bvec[b_index]**2)*vvec[v_index]*tf
        V[f'{vvec[v_index]}'].append(Vcapture)

Vcapture = []

print('bmax:', bmax/1.5e11)

for v in V.keys():
    Vtot = np.sum(np.array(V[f'{v}']))
    Vcapture.append([f'{v}', Vtot])

v0 = np.array([float(Vcapture[i][0]) for i in range(len(Vcapture))])
Vvol = np.array([2.987e-37*Vcapture[i][1] for i in range(len(Vcapture))])

data = np.array(list(zip(Vvol, v0)))
np.savetxt("output.txt", data, delimiter=" ", fmt="%.6e", header="Vvol v0")
