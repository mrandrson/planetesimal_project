import numpy as np
import os
import fast_histogram
import matplotlib.pyplot as plt

energy = []

for f in os.listdir('.'):
    if f.endswith('.txt'):
        try:
            dat = np.loadtxt(f)
            if np.isscalar(dat):
                energy.append(dat)
            elif isinstance(dat, np.ndarray):
                energy.extend(dat.flatten())
        except Exception as e:
            print(f"Error reading {f}: {e}")

energy = np.array(energy)

'''
#Linear Plotting
num_bins = 100
range_min = min(energy)
range_max = max(energy)

counts = fast_histogram.histogram1d(energy, bins=num_bins, range=(range_min, range_max))

bin_edges = np.linspace(range_min, range_max, num_bins + 1)
bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

plt.clf()
plt.yscale('log')
plt.plot(bin_centers, counts, drawstyle='steps-mid')
plt.xlabel('Energy')
plt.ylabel('Count')
plt.show()
'''

#Log Plotting

lin_threshold = 0.01
log_min = np.log10(lin_threshold)
log_max = np.log10(np.max(np.abs(energy)) + 1)

num_bins_log = 40
log_bins_pos = np.logspace(log_min, log_max, num_bins_log + 1)
log_bins_neg = -np.logspace(log_min, log_max, num_bins_log + 1)[::-1]

lin_bins = np.linspace(-lin_threshold, lin_threshold, 11)
bins = np.concatenate([log_bins_neg, lin_bins, log_bins_pos])

counts, bin_edges = np.histogram(energy, bins=bins)

bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

plt.figure()
plt.stairs(counts, bin_edges, fill=False, color='k')
plt.xscale('symlog', linthresh=lin_threshold, linscale=1)
plt.yscale('log')
plt.xlabel('Energy')
plt.ylabel('Count')
plt.show()
