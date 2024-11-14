'''
This script is meant as a test for GEOCLIM (COMBINE) boundary conditions
before running the full model.
It computes the age model of GEOCLIM basins based on their dimensions and
on the exchange matrix, both taken from COMBINE input files.
'''

import numpy as np
from matplotlib import pyplot as plt


# ++++++++++++++++++++++++++++++++++++++++++++++ #
# Which COMBINE input repertory + name of figure #
# ++++++++++++++++++++++++++++++++++++++++++++++ #

# historical BC
#DIR = '../../INPUT/COMBINE/historical/'
#XFILE = 'exchange_2.dat'
#NAME = 'GEOCLIM-hist'
#
# revised idealized BC
#DIR = '../../INPUT/COMBINE/idealized/'
#XFILE = 'exchange_2.dat'
#NAME = 'GEOCLIM-idl'
#
# IPSL-CM5A2 Pre-ind BC
#DIR = '../../INPUT/COMBINE/IPSL-PI-ref/'
#XFILE = 'exchange_2.dat'
#NAME = 'IPSL-PI-ref'
#
# IPSL-CM5A2 90Ma split-epic BC
#DIR='../../INPUT/COMBINE/IPSL-90Ma-splitepic-2X/'
#XFILE = 'Orb-Mean_exchange.dat'
#NAME = 'IPSL_90Ma_single-midlat'
#
# IPSL-CM5A2 90Ma 3bas-Arct3-epol
#DIR='../../INPUT/COMBINE/IPSL-90Ma-3bas-Arct3-2X/'
#XFILE = 'Orb-Mean_exchange.dat'
#NAME = 'IPSL_90Ma_3bas-Arct3'
#
# IPSL-CM5A2 90Ma 3bas-Arct3-epol-APx0.25
#DIR='../../INPUT/COMBINE/IPSL-90Ma-3bas-Arct3-2X/'
#XFILE = 'Orb-Mean_exchange_APx0.25.dat'
#NAME = 'IPSL_90Ma_3bas-Arct3-AP025'
#
# IPSL-CM5A2 90Ma 3bas-Arct3-epol-APx0
DIR='../../INPUT/COMBINE/IPSL-90Ma-3bas-Arct3-2X/'
XFILE = 'Orb-Mean_exchange_APx0.dat'
NAME = 'IPSL_90Ma_3bas-Arct3-AP0'


# To add a name to boxes
def basin_name(i):
    '''
    Function to customize: return a name associated to the basin number "i"
    '''
    return ''


#############
# Load data #
#############

srf = np.loadtxt(DIR+'indice_surface.dat')[:-1].astype(bool)
vol = 1e15*np.loadtxt(DIR+'oce_vol.dat')[:-1]
Xch = 365.2422*24*60*60*1e6*np.loadtxt(DIR+XFILE, delimiter='\t', max_rows=vol.size+1, dtype=float)[:-1,:-1].transpose()


############################
# Advection-aging equation #
############################

inv_vol = 1/vol

dt = 1e-1 # year

t = np.arange(0, 1e4, dt)
age = np.zeros((vol.size,t.size), vol.dtype)

# resolution
for i in range(1, t.size):
    age[:,i] = age[:,i-1] + dt*(inv_vol*(np.matmul(Xch.transpose(), age[:,i-1]) - Xch.sum(1)*age[:,i-1]) + 1)
    age[srf,i] = 0


############
# Plotting #
############

fig, ax = plt.subplots()
for i in range(vol.size):
    if not srf[i]:
        ax.plot(t, age[i,:], label='box #'+str(i+1)+basin_name(i))

ax.set_xlabel('model time (yr)')
ax.set_ylabel('box age (yr)')
fig.suptitle('age model -- '+NAME)
ax.legend()
fig.savefig('age_model_'+NAME+'.pdf', format='pdf', transparent=True)
plt.show()


##########
# Saving #
##########

np.savetxt('age_model_'+NAME+'.dat', age[:,-1])

