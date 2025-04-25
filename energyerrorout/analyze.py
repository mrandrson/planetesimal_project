import pandas as pd
import os
import numpy as np
import fast_histogram as fh
import matplotlib.pyplot as plt
import re
import math

cs = math.sqrt(3.36)*100
bvec = np.logspace(9, 17, 80)
Nv = 30
vvec = np.logspace(np.log10(cs/10), np.log10(cs*5), Nv)
tf = 1.0829*10.0**11 #s

V = {}

for v in vvec:
    V[f'{v}'] = [0]

for f in os.listdir('.'):
    if f.endswith('.csv'):
        dat = pd.read_csv(f)
        pattern = r"v(\d+)b(\d+)"

        match = re.search(pattern, f)

        if match:
            v_index = int(match.group(1))
            b_index = int(match.group(2))
            print(f"v_index: {v_index}, b_index: {b_index}")
        else:
            print("No match found")

        Ncapture = 0

        for p in dat['particle_id'].unique():
            pdat = dat[dat['particle_id'] == p]
            
            if any(dat['energy']<0):
                Ncapture += 1
                #print(f'b: {pdat[impact paramter]}')

        if Ncapture > 0:
            print("Captured:", Ncapture)
            print(f'bmax: {bvec[b_index+1]}, bmin: {bvec[b_index]}')

        Vcapture = np.pi*Ncapture*(bvec[b_index+1]**2-bvec[b_index]**2)*vvec[v_index]*tf
        V[f'{vvec[v_index]}'].append(Vcapture)

Vcapture = []

#print(V)

for v in V.keys():
    Vtot = np.sum(np.array(V[f'{v}']))
    Vcapture.append([f'{v}', Vtot])


v0 = np.array([float(Vcapture[i][0]) for i in range(len(Vcapture))])
Vvol = np.array([Vcapture[i][1] for i in range(len(Vcapture))])

#plt.yscale('log')
plt.plot(v0, Vvol)
plt.show()
