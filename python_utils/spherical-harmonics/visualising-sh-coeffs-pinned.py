import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import sys

# this is a simple script for quick visualization of the power spectrum 
# as a function of time and for the outputs generated with oscillations-array.sh


dir = sys.argv[1]
if(not dir[-1] == '/'):
    dir = dir +'/'
coeffs = genfromtxt(dir+'sh-coeffs.csv',delimiter=';')[:,:]
times = genfromtxt(dir+'times.csv',delimiter=';')[:,2]
print(times.shape,coeffs.shape)

fig = plt.figure()
ax = fig.add_subplot(2,1,1)

for i in range(32):
    ax.plot(times,coeffs[:,i])
#ax.set_yscale('log') # sometimes the log view is clearer
ax.legend(list(map(str,range(16))))
ax.set_ylim([0.0,0.2])
ax.set_title(f"zonal")



ax = fig.add_subplot(2,1,2)

for i in range(32,64):
    ax.plot(times,coeffs[:,i])
#ax.set_yscale('log') # sometimes the log view is clearer
ax.legend(list(map(str,range(16))))
ax.set_ylim([0.0,0.2])
ax.set_title(f"sectoral")

"""
ax = fig.add_subplot(1,2,2)

for i in range(nmodes):
    ax.plot(times,coeffs[:,32+i])
#ax.set_yscale('log') # sometimes the log view is clearer
ax.legend(list(map(str,labelmodes)))

ax.set_title(f"sectoral")
"""
plt.show()
