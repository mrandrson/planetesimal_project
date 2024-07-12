import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import math
x = np.zeros(99)
y = np.zeros(99)
r = np.zeros(99)
r0 = np.zeros(99)
rmaxvec = np.zeros(99)
rinit = np.zeros(99)
t = np.zeros(99)
script_R = 3.36*10**3
T = 10 ##Kelvin
B = 8.86
c_s = math.sqrt(script_R*T)
M = 1.989*10.0**30 ## kg
G = 6.674*10.0**(-11)  ## m^3*kg^-1*s^-2

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
for i in range(99):
    x[i] = np.loadtxt('./testposition/' + f'x{i}.txt')
    y[i] = np.loadtxt('./testposition/' + f'y{i}.txt')
    t[i] = np.loadtxt('./testtime/'+f't{i}.txt')
    r[i] = np.sqrt(x[i]**2+y[i]**2)
    r0[i] = c_s*t[i]
    #rinit[i] = r[1]
    rmaxvec[i] = au_to_meters(rmax)

plt.plot(t, r, label = '$r$')
#plt.plot(t, r0, label = '$r_0$')
plt.plot(t, rmaxvec, label = '$r_{max}$')
#plt.plot(t[62], r[62], 'ro', markersize=6)
#plt.plot(t[63], r[63], 'ro', markersize=6)
#plt.plot(t, rinit, label = '$r_{init}$')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('t(seconds)')
plt.ylabel('Position(meters)')
plt.legend()
plt.show()
