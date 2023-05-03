import numpy as np
from numpy import genfromtxt
import matplotlib.pyplot as plt
import matplotlib.colors as colors

# this is a simple script for quick visualization of the power spectrum 
# as a function of time and for the outputs generated with oscillations-array.sh

nmodes = 10
labelmodes = range(nmodes)

pressures = np.arange(13e4,14.1e4,1e4)
radii = np.arange(60e-6,65.1e-6,5e-6)

for p in pressures:
    for r in radii:
        print(f"results-{int(p)}-Pa-{r:.6f}-m")


        dir = f"../../build/oscillation-results/results-{int(p)}-Pa-{r:.6f}-m/"
        coeffs = genfromtxt(dir+'sh-coeffs.csv',delimiter=';')[:,:nmodes]
        times = genfromtxt(dir+'times.csv',delimiter=';')[:,2]
        print(times.shape,coeffs.shape)

        fig = plt.figure()
        ax = fig.add_subplot()

        for i in range(nmodes):
            ax.plot(times,coeffs[:,i])
        #ax.set_yscale('log') # sometimes the log view is clearer
        ax.legend(list(map(str,labelmodes)))

        ax.set_title(f"{int(p)} Pa | {r:.6f} m")

        plt.show()
