import numpy as np
import pyshtools as pysh  # pip install pyshtools makes this possible!
from numpy import genfromtxt

import sys

from multiprocessing import Pool

S = [None]*(len(sys.argv)-2)

def generate_sh_coeffs(i):
    """ this function generates spherical harmonics power spectrum of
        sample data given in the files mesh-#.csv """
    
    filename = sys.argv[i]
    print(filename)
    posa = filename.rfind('mesh-')
    posb = filename.rfind('.csv')

    ind = int(filename[posa+5:posb])

    x = genfromtxt(filename,delimiter=';').T

    y = pysh.expand.SHExpandDH(x,sampling=2)

    modulus = np.sqrt(y[0]**2 + y[1]**2)

    
    return ind,np.concatenate((modulus[:,0],np.diagonal(modulus))),i-2

# code that runs power spectrum generation in parallel for all 
# files given in the arguments (except the first two). The second
# argument gives the filepath to which the results are written.
output = 0
with Pool(48) as p:
    output = p.map(generate_sh_coeffs, range(2,len(sys.argv)))

for i,out in enumerate(output):
    S[out[2]] = out[1]

# write output
S = np.asarray(S)
np.savetxt(sys.argv[1],S,delimiter=';')
print('spherical-harmonics.py has finished.')
