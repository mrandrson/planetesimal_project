import numpy as np
import os
maxEnergy = 0
os.chdir('ShuEnergies')
s = os.listdir('.')
for i in s:
    m = 0
    energies = np.loadtxt(i)
    if energies.size == 1:
        m = float(energies)
    if energies.size > 1:
        m = max(energies)
    if m > maxEnergy:
        maxEnergy = m

print(maxEnergy)
