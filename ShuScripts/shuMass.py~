import shuSolution as sS
import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as sc
from LarsonSolution import trapIntegrateLinear, trapIntegrateLog
import scipy
c_s = np.sqrt(3.36 * 10.0 ** 3 * 10) * (31557600 / 149597870691)  # AU/YR
c_sms =np.sqrt(3.36 * 10.0 ** 3 * 10) 
r_values = [10, 100, 1000]  # AU
t_values = [10000, 100000, 350000]
def Mp(t):
    return (.98 * c_sms**3 * t)*(31557600) / (sc.G) #kg



def M_enc(r, t):
	x = r/(c_s*t)
	if x < 1:
		integral = trapIntegrateLinear(lambda x: x**2 * sS.alpha(x), 1*10**(-12), r / (c_s*t), 100)
	else:
		integral = trapIntegrateLinear(lambda x: x**2 * sS.alpha(x), 1*10**(-12), 1, 100) + trapIntegrateLog(lambda x: x**2 * sS.alpha(x), 1, r / (c_s*t), 100)
	return ((c_sms**3 * t) / (sc.G)) * integral * (31557600) #kg

plt.figure()

'''
for i, r in enumerate(r_values):
    t = np.logspace(2, 6, 100)
    M_values = np.zeros(len(t))
    for j in range(0, len(t)):
        M_values[j] = Mp(t[j])+M_enc(r, t[j]) #kg
    plt.plot(t, M_values, label=f'r = {r} AU')
'''

for i, t in enumerate(r_values):
    r = np.logspace(0, 4, 100)
    M_values = np.zeros(len(r))
    for j in range(0, len(r)):
        M_values[j] = Mp(t)+M_enc(r[j], t) #kg
    plt.plot(r, M_values, label=f't = {t} yr')

#print(scipy.optimize.brentq(lambda r: (Mp(10000)+M_enc(r, 10000)-2*10**30), .4, 20000))
plt.xscale("symlog")
plt.yscale("symlog")
plt.xlabel('r(AU)')
plt.ylabel('$M_{tot}(kg)$')
plt.legend()

plt.show()

