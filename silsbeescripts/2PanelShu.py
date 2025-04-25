import numpy as np
import math
import matplotlib.pyplot as plt
import os

c_s = 100*math.sqrt(3.36)
AU = 1.496*10**11
#Rmax = 1.602*10.0**15 # t_init = -100000
#Rmax = 6.19668*10.0**14 # t_init = -10000
#Rmax = 4.70168*10.0**14 # t_init = -1000



    
v_vec = np.logspace(np.log10(c_s/10), np.log10(c_s*5), 20)

def getVandB(f):
    v, b = f[1:-4].split('b')
    return float(v), float(b)



def get_volume_vec(t_init_mag, folder, b_bounds_vec, Nsims):
    Vol_vec = np.zeros(len(v_vec))
    for i in range(0, len(b_bounds_vec)-1):
        for j in range(0, len(v_vec)):
            s = float(np.loadtxt(folder + '/v' + str(j) + 'b' + str(i) + '.txt'))
            v = v_vec[j]
            Vol_vec[j] = Vol_vec[j] + s*math.pi*(b_bounds_vec[i+1]**2 - b_bounds_vec[i]**2)*v_vec[j]*t_init_mag/Nsims
    return np.divide(Vol_vec, (1.496*10**11)**3)

def getNumbers(v_index):
    f = np.zeros(len(b_bounds_vec)-1)
    for i in range(0, len(b_bounds_vec)-1):
        f[i] = float(np.loadtxt(folder + '/v' + str(v_index) + 'b' + str(i) + '.txt'))
    return f

b_bounds_vec_1 = np.logspace(9, 17, 40)
b_bounds_vec_2 = np.logspace(9, 17, 40)
b_bounds_vec_3 = np.logspace(9, 17, 40)
Volume_vec_1 = get_volume_vec(1.0829*10**11, './Shutran/Shu1', b_bounds_vec_1, 30)
Volume_vec_2 = get_volume_vec(1.0829*10**12, './Shutran/Shu2', b_bounds_vec_2, 30)
Volume_vec_3 = get_volume_vec(1.0829*10**13, './Shutran/Shu3', b_bounds_vec_3, 30)


def makeFigure():
    fig = plt.figure(figsize = (11, 6))
    ax1 = fig.add_subplot(121)
    l1, = ax1.plot(np.divide(v_vec, c_s), Volume_vec_1, color = 'red')
    l2, = ax1.plot(np.divide(v_vec, c_s), Volume_vec_2, color = 'green')
    l3, = ax1.plot(np.divide(v_vec, c_s), Volume_vec_3, color = 'blue')
    #plt.legend([l1, l2, l3], ['$R_{\\rm fin} = R_{\\rm M_\odot}/100$', '$R_{\\rm fin} = R_{\\rm M_\odot}/10$', '$R_{\\rm fin} = R_{\\rm M_\odot}$'], loc = "upper right", prop={'size':17.0}, ncol = 1, numpoints = 5, handlelength = 3.5)
    ax1.set_xlabel('$v/c_s$', fontsize = 20)
    ax1.set_ylabel('capture volume, AU$^3$', fontsize = 20)
    ax1.tick_params(axis='both', labelsize=20)
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    
    ax2 = fig.add_subplot(122)
    l1, = ax2.plot(np.divide(v_vec, c_s), np.divide(Volume_vec_1, 4*math.pi/3*(1.985*10**13/AU)**3), color = 'red')
    l2, = ax2.plot(np.divide(v_vec, c_s), np.divide(Volume_vec_2, 4*math.pi/3*(1.985*10**14/AU)**3), color = 'green')
    l3, = ax2.plot(np.divide(v_vec, c_s), np.divide(Volume_vec_3, 4*math.pi/3*(1.985*10**15/AU)**3), color = 'blue')
    plt.legend([l1, l2, l3], ['$R_{\\rm fin} = R_{\\rm M_\odot}/100$', '$R_{\\rm fin} = R_{\\rm M_\odot}/10$', '$R_{\\rm fin} = R_{\\rm M_\odot}$'], loc = "lower left", prop={'size':16.0}, ncol = 1, numpoints = 5, handlelength = 3.5)
    ax2.set_xlabel('$v/c_s$', fontsize = 20)
    ax2.set_ylabel('normalized capture volume', fontsize = 20)
    ax2.tick_params(axis='both', labelsize=20)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    fig.subplots_adjust(wspace = .4)
    fig.subplots_adjust(bottom = .12)
    fig.subplots_adjust(top = .98)
    fig.subplots_adjust(left = .12)
    fig.subplots_adjust(right = .98)
    plt.show()
    #plt.savefig('ShuCaptureRateRevised.png')


makeFigure()
      


    
            
        
    
