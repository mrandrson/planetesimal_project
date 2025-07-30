import os
import glob
import numpy as np
from fast_histogram import histogram1d
import matplotlib.pyplot as plt

dirs1 = glob.glob('./*/ShuEnergies*')
dirs2 = glob.glob('./*/dataShu*')

singlepe1 = []
singlepe2 = []

for d in dirs1:
    for f in os.listdir(d):
            if f.endswith('.txt'):
                    full_path = os.path.join(d, f)
                    dat = np.loadtxt(full_path)
                    dat = np.array(dat)
                    if dat.size == 0:
                            continue
                    if dat.size == 1:
                            singlepe1.append(dat)

for d in dirs2:
    for f in os.listdir(d):
            if f.endswith('.txt'):
                    full_path = os.path.join(d, f)
                    dat = np.loadtxt(full_path)
                    dat = np.array(dat)
                    if dat.size == 0:
                            continue
                    if dat.size == 1:
                            singlepe2.append(dat)

singlepe1 = np.array(singlepe1)
singlepe2 = np.array(singlepe2)

data1 = singlepe1.flatten()
data2 = singlepe2.flatten()

def log_histogram(data, bins, label):
    data = np.array(data)
    if data.size == 0:
        print("log_histogram: no data to plot")
        return
    positive_data = data[data > 0]
    negative_data = -data[data < 0]

    min_val = min(positive_data.min(), negative_data.min())
    max_val = max(positive_data.max(), negative_data.max())

    log_bins = np.logspace(np.log10(min_val), np.log10(max_val), num=bins+1)

    neg_bins = -log_bins[::-1]
    all_bins = np.concatenate([neg_bins, [0], log_bins])

    counts, edges = np.histogram(data, bins=all_bins)

    plt.hist(data, bins=all_bins, histtype='step', label = label)

log_histogram(data1, bins=15, label=None)
log_histogram(data2, bins=15, label=None)
plt.legend()
plt.xscale('symlog', linthresh = 1e-4)
plt.yscale('log')
plt.ylabel('Counts')
plt.xlabel(r'Energy $\left[\frac{m^2}{s^2}\right]$')
plt.show()

