import numpy as np
import math
import warnings
import os
import re
import matplotlib.pyplot as plt
cs = math.sqrt(3.36)*100

bvec = np.logspace(9, 18, 60)
Nv = 30
vvec = np.logspace(np.log10(cs / 100), np.log10(cs*5), Nv)

outdir = './data/ShuEnergies'
warnings.filterwarnings(
    "ignore",
    message="loadtxt: input contained no data*",
    category=UserWarning
)

energies = sorted(os.listdir(outdir))
deltaE = []
Efinal = []

for e in energies:
    path = os.path.join(outdir, e)
    evec = np.loadtxt(path)

    pattern = r"v(\d+)b(\d+)"

    match = re.search(pattern, e)

    if match:
        v_index = int(match.group(1))
        b_index = int(match.group(2))
    else:
        print("No match found")

    v0 = vvec[v_index]

    print(f"v_index: {v_index}, b_index: {b_index}")

    if evec.size == 0:
        continue

    if evec.size == 1:
        deltaE.append([evec-0.5*v0**2])
        Efinal.append([evec])
        continue

    deltaE.append(evec-0.5*v0**2)
    Efinal.append(0.5*v0**2*(evec-1))

deltaE = np.concatenate(deltaE)
Efinal = np.concatenate(Efinal)
def log_histogram(data, bins):
    data = np.array(data)
    positive_data = data[data > 0]
    negative_data = -data[data < 0]

    min_val = min(positive_data.min(), negative_data.min())
    max_val = max(positive_data.max(), negative_data.max())

    log_bins = np.logspace(np.log10(min_val), np.log10(max_val), num=bins+1)

    neg_bins = -log_bins[::-1]
    all_bins = np.concatenate([neg_bins, [0], log_bins])

    counts, edges = np.histogram(data, bins=all_bins)

    plt.figure(figsize=(8, 5))
    plt.hist(data, bins=all_bins, histtype='step')


vmin, vmax = round(min(vvec)/cs, 2), round(max(vvec)/cs, 1)
bmin, bmax = int(np.log10(min(bvec))), int(np.log10(max(bvec)))


np.savetxt(f"v{round(vmin)}_{round(vmax)}b{int(bmin)}_{int(bmax)}deltae.txt", deltaE)
log_histogram(Efinal, 30)
#plt.hist(deltaE, bins=30)
plt.xscale('symlog', linthresh=1e-1)
plt.yscale('log')
plt.xlabel(r"$\Delta E \left(\frac{m^2}{s^2}\right)$")
plt.xlim(3*min(Efinal), 2*max(Efinal))
plt.ylabel("Counts")
plt.title(rf"Energy Change for ${vmin}c_s \leq v_0 \leq {vmax}c_s$ and $10^{{{bmin}}} \leq b \leq 10^{{{bmax}}}$")
#plt.show()
plt.savefig(f"/Users/richardanderson/Downloads/v{vmin}_{vmax}b{bmin}_{bmax}deltae.png")
