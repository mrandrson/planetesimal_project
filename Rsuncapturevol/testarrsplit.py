class SplitArray:
    def __init__(self, array, m):
        self.array = array
        self.m = m
        self.n = len(array)
        self.chunk_size = self.n // m
        self.remainder = self.n % m

    def __getitem__(self, index):
        if index >= self.m or index < -self.m:
            raise IndexError("Split index out of range")

        if index < 0:
            index += self.m

        start = index * self.chunk_size + min(index, self.remainder)
        end = start + self.chunk_size
        if index < self.remainder:
            end += 1
        return self.array[start:end]

    def __len__(self):
        return self.m

import numpy as np
import matplotlib.pyplot as plt

def plotmean(x, y, nbins, vlines=False):
    xbins = SplitArray(x, nbins)
    
    def yinbin(bin):
        yi = np.linspace((bin)*(len(x)/len(xbins)),(bin+1)*(len(x)/len(xbins)), int(len(x)/len(xbins)+1))
        yvec = np.array([y[int(yii)] for yii in yi])
        return yvec

    if vlines == True:
        for bin in xbins: 
            plt.axvline(min(bin), linestyle = '--', color = 'red', alpha = 0.2)
            plt.axvline(max(bin), linestyle = '--', color = 'red', alpha = 0.2)

    ymean = [np.mean(yinbin(bindex)) for bindex in range(len(xbins)-1)]
    xmid = [(max(xbins[bindex])+min(xbins[bindex]))/2 for bindex in range(len(xbins)-1)]
    plt.plot(xmid, ymean)

def plotmin(x, y, nbins):
    xbins = SplitArray(x, nbins)
    
    def yinbin(bin):
        yi = np.linspace((bin)*(len(x)/len(xbins)),(bin+1)*(len(x)/len(xbins)), int(len(x)/len(xbins)+1))
        yvec = np.array([y[int(yii)] for yii in yi])
        return yvec

    ymin = [min(yinbin(bindex)) for bindex in range(len(xbins)-1)]
    xmid = [(max(xbins[bindex])+min(xbins[bindex]))/2 for bindex in range(len(xbins)-1)]
    plt.plot(xmid, ymin)

x = np.linspace(0, 100, 1000)

y = [np.sqrt(np.random.random())*(2*xi-xi**2+np.exp(0.09*xi)) for xi in x]
f = lambda x: 2*x-x**2+np.exp(0.09*x)
plt.clf()
#plt.plot(x, f(x))
plt.scatter(x, y, alpha =0.1)
plotmean(x, y, 20, vlines=True)
plotmin(x, y, 20)
plt.show()
