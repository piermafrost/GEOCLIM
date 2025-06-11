'''
Generate a "paleo" slope field, with the following inputs:
* paleo topography (meter)
* paleo tectonics division (masks of all categories, same grid)
* modern topography (meter)
* modern tectonics division (similar to paleo)
* modern slope

The algorithm use the tectonics division, and add to it an elevation division
(typically: lower than 500m, 500-1000m, 1000-2000m, 2000m and higher).
Then, each paleo grid cell is assigned a slope value randomly selected among the
grid cells of same tectonics and elevation category on the modern slope field.

All the modern files (topography, categories masks and slope) must have the same
geographic grid, and similarly, all the paleo files must have the same grid.
Ideally, the modern and paleo grid should be the same, but it is not necessary.

All the files must be netCDF format.

The user is expected to indicate all input data in the block
"USER CUSTOMIZABLE DATA AND PARAMETERS"
'''

# Required packages
import netCDF4 as nc
import numpy as np

# Optional package (for plotting)
try:
    from matplotlib import pyplot as plt
    PLOT_PACKAGE = True
except ModuleNotFoundError:
    PLOT_PACKAGE = False



## <><><><><><><><><><><><><><><><><><><><> ##
## <><><><><><><><><><><><><><><><><><><><> ##
##                                          ##
##  USER CUSTOMIZABLE DATA AND PARAMETERS   ##
##                                          ##
## <><><><><><><><><><><><><><><><><><><><> ##
## <><><><><><><><><><><><><><><><><><><><> ##


# Wether or not plotting "check" figures at each step
# ---------------------------------------------------

CHECK_PLOT = True


# Name of output slope file
# =========================

OUT_SLOPE_FILE = 'slope_Mid-Cretaceous_90Ma.nc'
FILLVAL = -9e33


# Elevation categories
# ====================

ELEV_SPLIT = (9.8*500, 9.8*1000, 9.8*2000) # splitting values => 4 categories
#             ^^^^     ^^^^      ^^^^      -->  because input data is PHIS (geopotential height, m2/s2), not elevation


# Paleo files
# ===========

# Elevation => reference file, only close it at the end!
fref = nc.Dataset('../INPUT/IPSL/90Ma_Laugie/2X/Orb1/CPL-90Ma-ORB1-2X_SE_6950_7049_1Y_atm.nc')
ZPALEO_NCVAR = fref['phis']
z_paleo = ZPALEO_NCVAR[0,:,:] # first dimension is "TIME_COUNTER" (size=1)

# Land-sea mask (expect >1 for lands and 0 for oceans)
f = nc.Dataset('../INPUT/IPSL/90Ma_Laugie/2X/Orb1/CPL-90Ma-ORB1-2X_SE_6950_7049_1Y_sechiba.nc')
lsm_paleo = f['Contfrac'][:,:]
f.close()


# Masks of all tectonics categories
masks_paleo = []
for file_name, var_name in zip(['example_tecto_split/active_orogen_90Ma_LMDz-96x96.nc',
                                'example_tecto_split/ancient_relief_90Ma_LMDz-96x96.nc',
                                'example_tecto_split/continent_90Ma_LMDz-96x96.nc'],
                               3*['polygon_fraction']):
    f = nc.Dataset(file_name)
    masks_paleo.append(f[var_name][:,:].data >= 0.5) 
    f.close()


# Modern files
# ============

# Elevation
f = nc.Dataset('../INPUT/IPSL/CTRL_CM5A2/piControl_SE_2750_2849_1Y_atm.nc')
z_modern = f['phis'][0,:,:] # first dimension is "TIME_COUNTER" (size=1)
f.close()

# Slope
f = nc.Dataset('../INPUT/slope/slope_PI_IPSL-96x96.nc')
SLPMODERN_NCVAR = f['slope']
slp_modern = SLPMODERN_NCVAR[:,:]
f.close()

# Land-sea mask (expect >1 for lands and 0 for oceans)
f = nc.Dataset('../INPUT/IPSL/CTRL_CM5A2/piControl_SE_2750_2849_1Y_sechiba.nc')
lsm_modern = f['Contfrac'][:,:]
f.close()


# Masks of all tectonics categories (must be the same categories as the paleo ones!)
masks_modern = []
for file_name, var_name in zip(['example_tecto_split/active_orogen_CTRL_LMDz-96x96.nc',
                                'example_tecto_split/ancient_relief_CTRL_LMDz-96x96.nc',
                                'example_tecto_split/continent_CTRL_LMDz-96x96.nc'],
                               3*['polygon_fraction']):
    f = nc.Dataset(file_name)
    masks_modern.append(f[var_name][:,:].data >= 0.5) 
    f.close()


## <><><><><><><><><><><><><><><><><><><><> ##
## <><><><><><><><><><><><><><><><><><><><> ##





#####################################
# Construction of paleo slope field #
#####################################


# Check-plots
# -----------

PLOT = (PLOT_PACKAGE and CHECK_PLOT)

try:
    ydim, xdim = ZPALEO_NCVAR.dimensions
except ValueError: 
    _, ydim, xdim = ZPALEO_NCVAR.dimensions # in case there is a degenerate time dimension

try:
    x, y = (fref[var][:] for var in (xdim, ydim))
except IndexError:
    x, y = None, None

def plotmap(var, ax=None):
    if ax is None:
        _, ax = plt.subplots()

    if x is None or y is None:
        pid = ax.pcolormesh(var)
    else:
        pid = ax.pcolormesh(x, y, var)

    plt.colorbar(pid, ax=ax)


# Mask non-lands points
# ---------------------

if lsm_modern is not None:
    lsm_modern[slp_modern.mask] = 0
    z_modern = np.ma.masked_where(lsm_modern==0, z_modern)

if lsm_paleo is not None:
    z_paleo = np.ma.masked_where(lsm_paleo==0, z_paleo)

if PLOT:
    fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(10,5))
    plotmap(z_modern, ax=ax[0])
    ax[0].set_title('elevation - MODERN')
    plotmap(z_paleo, ax=ax[1])
    ax[1].set_title('elevation - PALEO')
    fig.tight_layout()
    plt.show()


# build masks of elevation categories
# -----------------------------------

elevmasks_modern = []
elevmasks_paleo  = []

elevmasks_modern.append(z_modern < ELEV_SPLIT[0])
elevmasks_paleo.append(z_paleo < ELEV_SPLIT[0])

z0 = ELEV_SPLIT[0]
for z1 in ELEV_SPLIT[1:]:
    elevmasks_modern.append(np.logical_and(z_modern >= z0, z_modern < z1))
    elevmasks_paleo.append(np.logical_and(z_paleo >= z0, z_paleo < z1))
    z0 = 1*z1

elevmasks_modern.append(z_modern >= ELEV_SPLIT[-1])
elevmasks_paleo.append(z_paleo >= ELEV_SPLIT[-1])


# Initialization of paleo slope field
# -----------------------------------

slp_paleo = np.ma.zeros(z_paleo.shape, dtype=slp_modern.dtype)
slp_paleo.mask = z_paleo.mask


# Assign slope value for each category
# ------------------------------------

# double loop on all tectonics/elevation categories
for i, (msk_mod, msk_pal) in enumerate(zip(masks_modern, masks_paleo)):
    for j, (elmsk_mod, elmsk_pal) in enumerate(zip(elevmasks_modern, elevmasks_paleo)):

        if PLOT:
            fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(10,5))
            plotmap(np.logical_and(msk_mod, elmsk_mod), ax=ax[0])
            ax[0].set_title('tect. categ. #{:0} & elev. categ. #{:1} - modern'.format(i, j))
            plotmap(np.logical_and(msk_pal, elmsk_pal), ax=ax[1])
            ax[1].set_title('tect. categ. #{:0} & elev. categ. #{:1} - paleo'.format(i, j))
            fig.tight_layout()
            plt.show()

        # get indices of current tecto/elev masks
        idx_modern = np.argwhere(np.logical_and(msk_mod, elmsk_mod))
        idx_paleo  = np.argwhere(np.logical_and(msk_pal, elmsk_pal))

        # put random selection of modern slope in paleo slope field
        if idx_paleo.shape[0] > 0:
            if idx_modern.shape[0] == 0:
                print(idx_paleo.shape, idx_modern.shape)
                raise ValueError('No modern land points found in tectonics category #{:0} and elevation category #{:1}'.format(i, j))
            else:
                kselect = np.random.randint(low=0, high=idx_modern.shape[0], size=idx_paleo.shape[0])
                slp_paleo[idx_paleo[:,0], idx_paleo[:,1]]  =  slp_modern[idx_modern[kselect,0], idx_modern[kselect,1]]


# Record paleo-slope field in netCDF file
# ---------------------------------------

# plot of modern and built paleo slope
if PLOT:
    fig, ax = plt.subplots(ncols=2, sharex=True, sharey=True, figsize=(10,5))
    plotmap(slp_modern, ax=ax[0])
    ax[0].set_title('slope - modern')
    plotmap(slp_paleo, ax=ax[1])
    ax[1].set_title('slope - paleo')
    fig.tight_layout()
    plt.show()


# create file
fout = nc.Dataset(OUT_SLOPE_FILE, mode='w', data_model='NETCDF3_CLASSIC')

# create dimensions and dimension variables
for dim in (ydim, xdim):
    fout.createDimension(dim, size=fref.dimensions[dim].size)
    try:
        fout.createVariable(dim, datatype=fref[dim].datatype, dimensions=(dim,))
        fout[dim][:] = fref[dim][:]
        for att in fref[dim].ncattrs():
            if att != '_FillValue':
                fout[dim].setncattr(att, fref[dim].getncattr(att))

    except IndexError:
        print('Warning: cannot find dimension variable "{:}"'.format(dim))

# set fill-value
slp_paleo.data[slp_paleo.mask] = FILLVAL

# create slope variable
fout.createVariable('slope', datatype=slp_paleo.dtype, dimensions=(ydim, xdim), fill_value=FILLVAL)
fout['slope'][:] = slp_paleo[:]
for att in SLPMODERN_NCVAR.ncattrs():
    if att != '_FillValue':
        fout['slope'].setncattr(att, SLPMODERN_NCVAR.getncattr(att))

# close output in reference files
fref.close
fout.close()


