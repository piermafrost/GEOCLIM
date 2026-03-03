'''
Computes the average of the 16 orbital configurations in the files
"GCM_oceanic_temp.dat" and "exchange_2.dat" in given COMBINE input directories.
Creates output files "Orb-Mean_GCM_oce_temp.dat", "Orb-Mean_exchange.dat" and
"Orb-Mean_exchange_human-friendly.dat" in all given input directories.

USAGE:
    `python compute_orbital_average.py INPUT/COMBINE/DIR1 INPUT/COMBINE/DIR2`
    (as many "DIR" as you wish)
'''

import numpy as np
import sys

list_rep = sys.argv[1:]

for rep in list_rep:

    # Temperature
    # -----------

    T = np.loadtxt(rep+'/GCM_oceanic_temp.dat')
    n,k = T.shape # (configs, boxes)
    if n == 32:
        Tave = np.zeros((2,k), dtype=T.dtype)
        Tave[0,:] = T[:n//2, :].mean(0)
        Tave[1,:] = T[n//2:, :].mean(0)
    else:
        Tave = T.mean(0)
        if n != 16:
            print('WARNING: found {0:} configs. This script was meant for 1x16 or 2x16 configs.'.format(n))
            print(' => average all {0:} configurations'.format(n))

    np.savetxt(rep+'/Orb-Mean_GCM_oce_temp.dat', Tave, fmt='%10.5f')

    # === #

    # Oceanic circulation
    # -------------------

    X = np.loadtxt(rep+'/exchange_2.dat')
    nk,k = X.shape # (configs*boxes, boxes)
    n = nk//k
    if n*k != nk:
        raise ValueError('Cannot understand a series of exchanges matrices of shape ({0:},{1:}).\nExpect (nconfigs*nboxes, nboxes).'.format(nk,k))

    # (configs, rows, colums)
    X = X.reshape((n,k,k))

    if n == 32:
        Xave = np.zeros((2*k,k), dtype=X.dtype)
        Xave[:k,:] = X[:n//2, :, :].mean(0)
        Xave[k:,:] = X[n//2:, :, :].mean(0)
    else:
        Xave = X.mean(0)
        if n != 16:
            print('WARNING: found {0:} configs. This script was meant for 1x16 or 2x16 configs.'.format(n))
            print(' => average all {0:} configurations'.format(n))

    np.savetxt(rep+'/Orb-Mean_exchange.dat', Xave, fmt='%.7e', delimiter='\t')
    np.savetxt(rep+'/Orb-Mean_exchange_human-friendly.dat', Xave, fmt='%6.2f')

