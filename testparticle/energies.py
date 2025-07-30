import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

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

def hist2d(x, y, bins_x, bins_y, logx=False, logy=False, cmap='viridis', figsize=(8,6)):
    def symmetric_log_edges(data, bins):
        data = np.array(data)
        pos = data[data>0]
        neg = -data[data<0]
        if len(pos)==0 or len(neg)==0:
            raise ValueError("Data must contain both positive and negative for symmetric log")
        minval = min(pos.min(), neg.min())
        maxval = max(pos.max(), neg.max())
        log_half = np.logspace(np.log10(minval), np.log10(maxval), bins+1)
        neg_edges = -log_half[::-1]
        pos_edges =  log_half
        return np.concatenate([neg_edges, [0.], pos_edges])

    if isinstance(bins_x, int):
        if logx:
            edges_x = symmetric_log_edges(x, bins_x)
        else:
            edges_x = np.linspace(np.min(x), np.max(x), bins_x+1)
    else:
        edges_x = np.array(bins_x)

    if isinstance(bins_y, int):
        if logy:
            edges_y = symmetric_log_edges(y, bins_y)
        else:
            edges_y = np.linspace(np.min(y), np.max(y), bins_y+1)
    else:
        edges_y = np.array(bins_y)

    H, xe, ye = np.histogram2d(x, y, bins=[edges_x, edges_y])

    fig, ax = plt.subplots(figsize=figsize)
    pcm = ax.pcolormesh(
        edges_x, edges_y, H.T,
        norm=LogNorm(vmin=H[H>0].min(), vmax=H.max()),
        cmap=cmap
    )
    cbar = fig.colorbar(pcm, ax=ax)
    cbar.set_label('Counts')

    ax.set_xlabel('x')
    ax.set_ylabel('y')

    if logx:
        ax.set_xscale('log', linthresh=min(abs(edges_x[edges_x>0].min()), abs(edges_x[edges_x<0].max())))
    if logy:
        ax.set_yscale('symlog', linthresh=min(abs(edges_y[edges_y>0].min()), abs(edges_y[edges_y<0].max())))


energies = np.loadtxt('energies.txt')
v0,e =energies[:,0], energies[:,1]
ke = (v0*v0)/2

''' 
hist2d(v0, e, bins_x=50, bins_y=50, logx=False, logy=True, cmap='viridis', figsize=(8, 6))
plt.yscale('symlog', linthresh=10)
plt.ylabel(r'Energy $\left(\frac{m^2}{s^2}\right)$')
plt.xlabel(r'$v_0\ \left(\frac{m}{s}\right)$')
'''  
  
  
log_histogram(e/ke, bins=50)
plt.title("C Simulation Energy Ratio")
plt.xscale('symlog', linthresh=1e-3)
plt.yscale('log')
plt.xlabel(r'$E_f/E_0$')
plt.ylabel('Counts')

plt.show()
