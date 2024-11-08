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
    au_in_meters = 149597870.7  # meters
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
##Is this Rmax?

def get_Shu_enclosed_mass(r, t):
    r0 = c_s*t
    centralMass = .975502*c_s**2*r0/G #Why doesn't the central mass change with time?
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
