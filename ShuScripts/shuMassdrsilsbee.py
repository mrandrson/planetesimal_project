'''
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math

script_R = 3.36*10**3
T = 10
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.9891*10.0**30 ## kg 
G = 6.67*10.0**(-11)  ## m^3*kg^-1*s^-2

def years_to_seconds(years):
    seconds_per_year = 365.25 * 24 * 3600
    return years * seconds_per_year

def au_to_meters(au):
    au_in_meters = 1.49e11  # meters
    return au * au_in_meters


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

def get_Shu_enclosed_mass(r, t):
    r0 = c_s*t
    centralMass = .975502*c_s**2*r0/G
    Mcalc =  centralMass + r0**3/(G*t**2)*get_shu_integral(r/r0)
    return min(Mcalc, M)

r_values = [10, 100, 1000]

tsec = np.logspace(-1, 6, 1000)

for i, r in enumerate([au_to_meters(r_values[i]) for i in [0, 1, 2]]):
    t = [years_to_seconds(tsec[i]) for i in range(len(tsec))]
    M_values = np.zeros(len(t))
    for j in range(0, len(t)):
        M_values[j] = get_Shu_enclosed_mass(r, t[j]) #kg
    plt.plot(t, M_values, label=f'r = {r_values[i]} AU')


plt.xlabel('t(yr)', fontsize = 20)
plt.ylabel('$M_{enc}(kg)$', fontsize = 20)
plt.legend(fontsize = 20)
plt.xscale('log')
plt.yscale('log')
plt.show()
'''
import numpy as np
import matplotlib.pyplot as plt
import scipy
import math

script_R = 3.36e3
T = 10
B = 8.86
c_s = math.sqrt(script_R * T)
M = 1.9891e31
G = 6.67e-11

def years_to_seconds(years):
    return years * 365.25 * 24 * 3600

def au_to_meters(au):
    return au * 1.49e11

shu_x = np.logspace(-12, np.log10(2), 10000)
shu_integral = np.loadtxt('shuInt.txt')
shu_helper = scipy.interpolate.interp1d(shu_x, shu_integral, kind='cubic')

def get_shu_integral(x):
    if x < 1e-12:
        return 0
    if x > 2:
        return shu_helper(2) + 2 * (x - 2)
    return shu_helper(x)

r_out = 1.974342e15

def get_Shu_enclosed_mass(r, t):
    r0 = c_s * t
    centralMass = 0.975502 * c_s**2 * r0 / G
    Mcalc = centralMass + r0**3 / (G * t**2) * get_shu_integral(r / r0)
    return min(Mcalc, M)

r_values_au = np.logspace(np.log10(10), np.log10(r_out/1.49e11), 100)
r_values = au_to_meters(r_values_au)
tsec = np.logspace(-1, 6, 1000)
t_values = np.array([years_to_seconds(t) for t in tsec])
M_matrix = np.zeros((len(r_values), len(t_values)))

for i, r in enumerate(r_values):
    for j, t in enumerate(t_values):
        M_matrix[i, j] = get_Shu_enclosed_mass(r, t)

T, R = np.meshgrid(tsec, r_values_au)
M_flat = M_matrix.flatten()
T_flat = T.flatten()
R_flat = R.flatten()

plt.figure(figsize=(10, 6))
sc = plt.scatter(T_flat, M_flat, c=R_flat, cmap='viridis', s=10, edgecolor='none')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('t (yr)', fontsize=14)
plt.ylabel(r'$M_{\rm enc}$ (kg)', fontsize=14)
cbar = plt.colorbar(sc, label='Radius (AU)')
plt.tight_layout()
plt.show()

