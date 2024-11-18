import numpy as np

DT = 1e3 # (yrs)
t_stop = 1e6 # (yrs)
t = np.arange(0, t_stop, DT)


# Common parameter ranges #
# ======================= #

ecc_period = 405e3 # (yrs)
ecc_range = (0.015, 0.06)
ecc_start_phase = 0. # (radians)
#
prec_period = 23e3 # (yrs)
prec_start  = 0. # (radians)
#
obl_period = 41e3 # (yrs)
obl_range = (22.1, 24.5) # (degrees)
obl_start_phase = 0. # (radians)


# high eccentricity, precession-only:
# -----------------------------------

ecc = 0.06
obl = 22.1 # (°) #=> low obliquity
prec = prec_start + 360*t/prec_period
params = np.array([ecc*np.ones(prec.shape), prec, obl*np.ones(prec.shape)]).transpose()
np.savetxt('ideal_ecc-0.06_prec-23kyr_obl-22.1_1Myrs.dat', params,
           fmt='%.15e', delimiter='\t', comments='#',
           header=' Idealized 23-kyr precession cycle, with eccentricity at 0.06 and fixed obliquity at 22.1°\n time-step: 1 kyr; length: 1 Myr\n eccentricity    \tprecession (°)   \tobliquity (°)')


# eccentricity and precession:
# ----------------------------

obl = 22.1 # (°) #=> low obliquity
ecc  = ecc_range[0] + (ecc_range[1] - ecc_range[0])*(1 - np.cos(2*np.pi*t/ecc_period + ecc_start_phase))/2
prec = prec_start + 360*t/prec_period
params = np.array([ecc, prec, obl*np.ones(prec.shape)]).transpose()
np.savetxt('ideal_ecc-405kyr_prec-23kyr_obl-22.1_1Myrs.dat', params,
           fmt='%.15e', delimiter='\t', comments='#',
           header=' Idealized 23-kyr precession cycle, 405-kyrs eccentricity cycle (sinusoidal in 0.015–0.06), with fixed obliquity at 22.1°\n time-step: 1 kyr; length: 1 Myr\n eccentricity    \tprecession (°)   \tobliquity (°)')


# low eccentricity, Sep Equinox perihel (ORB7); obliquity-only 
# ------------------------------------------------------------

ecc = 0.015
prec = 0.
obl  = obl_range[0] + (obl_range[1] - obl_range[0])*(1 - np.cos(2*np.pi*t/obl_period + obl_start_phase))/2
params = np.array([ecc*np.ones(obl.shape), prec*np.ones(obl.shape), obl]).transpose()
np.savetxt('ideal_ecc-0.015_prec-0_obl-41kyr_1Myrs.dat', params,
           fmt='%.15e', delimiter='\t', comments='#',
           header=' Idealized 41-kyrs obliquity cycle (sinusoidal in 22.1°–24.5°)\n time-step: 1 kyr; length: 1 Myr\n eccentricity    \tprecession (°)   \tobliquity (°)')


# All three parameters
# --------------------

ecc  = ecc_range[0] + (ecc_range[1] - ecc_range[0])*(1 - np.cos(2*np.pi*t/ecc_period + ecc_start_phase))/2
prec = prec_start + 360*t/prec_period
obl  = obl_range[0] + (obl_range[1] - obl_range[0])*(1 - np.cos(2*np.pi*t/obl_period + obl_start_phase))/2
params = np.array([ecc, prec, obl]).transpose()
np.savetxt('ideal_ecc-405kyr_prec-23kyr_obl-41kyr_1Myrs.dat', params,
           fmt='%.15e', delimiter='\t', comments='#',
           header=' Idealized 23-kyr precession cycle, 405-kyrs eccentricity cycle (sinusoidal in 0.015–0.06), 41-kyrs obliquity cycle (sinusoidal in 22.1°–24.5°)\n time-step: 1 kyr; length: 1 Myr\n eccentricity    \tprecession (°)   \tobliquity (°)')

