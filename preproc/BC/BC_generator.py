'''
Functions to generate GEOCLIM oceanic boundary condition files from netCDF
ouputs of ocean module of any climate model.

The central function of this file is "oce_output_2_GEOCLIM_BC"

See file "examples_oceanic_input.py" for examples of use.
'''

import netCDF4 as nc
import numpy as np
from units import temperature_units, length_units, area_units, volume_units, latitude_units, velocity_units, flux_units, UnknownUnitError
from local_functions import remove_orca_north_fold, find_bathy, get_var, get_axis_weighting, average, minimum
from conservative_flux_correction import PositiveLeastSquare_correction, RelativeLeastSquare_correction
from itertools import cycle
import os




###############################
# SPECIFIC GEOCLIM PARAMETERS #
###############################

# GEOCLIM's default oceanic basin definition
# ------------------------------------------
NBASIN = 9 # not counting atmosphere box
H_EPICONT = 200  # (m) maximum ocean floor depth (bathymetry) of the epicontinental box
H_SURF    = 100  # (m) depth of surface boxes (ie, starting depth of thermocline or epicont deep box)
H_THERMO  = 1000 # (m) depth of thermocline boxes (ie, starting depth of deep box)
MIDLAT_RANGE = (-60, 60) # (degrees N): latitudinal range of inner ocean (ie, excluding high-latitude oceans)

# Modern total ocean volume and areas
# -----------------------------------
MODERN_VOLUME = 1.335e18 # m3
MODERN_AREA   = 3.619e14 # m2

# Physical parameters
# -------------------
RHO_x_G = 1027*9.806

# Parameterization of vertical mixing
# -----------------------------------
W_PARAM_ADVECTIVE_MIXING = 8e-8 #3e-8 # m/s
# => Tuned to get a bidirectional vertical mixing flux of ~30 Sv world-wide
#    This value fit the oceanic age model

# Arbitrary and and technical GEOCLIM parameters (must not be modified)
# ---------------------------------------------------------------------
ATM_VOL = 1.
ATM_AREA = 0.363e15
ATM_SEDI_AREA = 0.357e15
DEFAULT_VOLUME = 1. # pseudo-volume for non-concentration variables
AREA_CONVERSION_FACTOR      = 1e-15 # areas expressed in 1e9 km2
VOLUME_CONVERSION_FACTOR    = 1e-15 # volume expressed in 1e6 km3
FLUX_CONVERSION_FACTOR      = 1e-6  # water fluxes expressed in Sv (1e6 m3/s)
DEPTH_CONVERSION_FACTOR     = 1e-3  # box depths expressed in km
PRESSURE_CONVERSION_FACTOR  = 1e-5  # pressure expressed in bar



##################################################
##----------------------------------------------##
##-- FUNCTION DEFINING MASK OF GEOCLIM BASINS --##
##----------------------------------------------##
##################################################

def geoclim_default_basin_mask(nav_lat, nav_z, nav_bathy, lat_split=MIDLAT_RANGE):
    '''
    Create a mask telling if points are in (T) or outside (F) each of the 9
    GEOCLIM basins, pre_defined by default. The mask is built on the grid
    "nav_z", and seafloor defined by the latitude array "nav_lat", the vertical
    column depth array depth (bathymetry) array "nav_bathy" (longitude array is
    not needed since there is no longitudinal split).
    The returned mask will rank-4, the 1st dimension being GEOCLIM basins, the
    last 3 being the dimensions of the 3-D input fields (i.e., latitude, z and
    bathymetry).

    nav_lat, nav_z and nav_bathy MUST BE RANK-3, even if they are not defined on
    all dimensions (e.g., latitude and bathymetry should be defined on
    horizontal, dimensions, whereas z should be defined on the vertical
    dimension). Use degenerated (size-1) dimensions for those extra dimensions.
    The 3 dimensions must obviously correspond (ie, be in the same order) to
    the dimensions of 3D oceanic variables (e.g., temperature).

    The GEOCLIM "historical" basin definition must be the following:
      0. N high-latitutde, surface (incl. thermo)
      1. N high-latitutde, deep
      2. open ocean, surface
      3. open ocean, thermocline
      4. open ocean, deep
      5. epicontinental, surface
      6. epicontinental, deep
      7. S high-latitude, surface (incl. thermo)
      8. S high-latitude, deep

    RETURN:

        mask: 4-D boolean numpy-array
            mask of all defined GEOCLIM boxes. First dimension is "box",
            remaining ones are 3D ocean grid dimensions.

        i_surf: 1-D integer numpy-array
            Indices of all surface boxes (used for intern computations in
            'oce_output_2_GEOCLIM_BC')

        deepbox_mask: 1-D integer numpy-array
            Whether or not each box is "deep" (1|0)

        surfbox_mask: 1-D integer numpy-array
            Whether or not each box is "surface" (1|0)

        thrmbox_mask: 1-D integer numpy-array
            Whether or not each box is "thermocline" ("intermediate") (1|0)

        polebox_mask: 1-D integer numpy-array
            Whether or not each box is "polar" (1|0)

        epicbox_mask: 1-D integer numpy-array
            Whether or not each box is "epicontinental" (1|0)

        appcbox_mask: 1-D integer numpy-array
            Whether or not each box "receives continental inputs" (1|0)

        box_depth: 1-D float numpy-array
            (bottom) depth of each box 
    '''

    # "lat_split" can be a scalar (=> symmetric splitting)
    if np.size(lat_split) == 1:
        lat_split = (-np.abs(lat_split), np.abs(lat_split))


    ##############################
    i_surf = np.array([0, 2, 5, 7])
    i_deep = np.array([1, 4, 8])
    i_N_HL = np.array([0, 1])
    i_S_HL = np.array([7, 8])
    i_epic = np.array([5, 6])
    ##############################

    deepbox_mask, surfbox_mask, thrmbox_mask, polebox_mask, epicbox_mask, appcbox_mask = np.zeros((6,NBASIN+1), dtype='int8')
    #                           vv----------------- one more basin for "atm box" --------------------------^^
    box_depth = np.zeros((NBASIN+1), dtype=float)

    deepbox_mask[i_deep] = 1
    surfbox_mask[i_surf] = 1
    polebox_mask[i_S_HL] = 1
    polebox_mask[i_N_HL] = 1
    epicbox_mask[i_epic] = 1
    thrmbox_mask[:] = 1
    thrmbox_mask[i_surf] = 0
    thrmbox_mask[i_deep] = 0
    thrmbox_mask[i_epic] = 0
    thrmbox_mask[-1] = 0 # atm box
    appcbox_mask[np.logical_and(surfbox_mask==1, epicbox_mask==1)] = 1


    shp = np.maximum(np.maximum(nav_lat.shape, nav_z.shape), nav_bathy.shape)

    mask = np.ones(np.concatenate(([NBASIN], shp)), dtype=bool)


    epicont_mask    = (nav_bathy <= H_EPICONT)
    surf_mask       = (nav_z <= H_SURF)
    surfthermo_mask = (nav_z <= H_THERMO)
    thermo_mask     = (nav_z > H_SURF)
    deep_mask       = (nav_z > H_THERMO)
    bottom_mask     = (nav_z <= nav_bathy)
    north_mask      = (nav_lat > lat_split[1])
    south_mask      = (nav_lat < lat_split[0])

    # Northern high-lat basins
    for i in i_N_HL:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], north_mask)

    # Inner ocean basins
    for i in range(NBASIN):
        if i not in i_N_HL and i not in i_S_HL and i not in i_epic:
            mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], np.logical_and(~north_mask, ~south_mask))

    # Southern high-lat basins
    for i in i_S_HL:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], south_mask)

    # Epicontinental basins
    for i in i_epic:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], epicont_mask)

    # Non epicontinental basins
    for i in range(NBASIN):
        if i not in i_epic:
            mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], ~epicont_mask)

    # Surface basins
    for i in i_surf:
        if i in i_N_HL or i in i_S_HL: # (high latitude)
            mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], surfthermo_mask)
            box_depth[i] = H_THERMO
        else:                          # (non high latitude)
            mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], surf_mask)
            box_depth[i] = H_SURF

    # Epicontinental deep basins and thermocline (intermediate) basins
    for i in range(NBASIN):
        if i not in i_surf and i in i_epic: # epic deep
            mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], thermo_mask)
            box_depth[i] = -1

        if i not in i_surf and i not in i_deep and i not in i_epic: # intermediate
            mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], thermo_mask)
            mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], surfthermo_mask)
            box_depth[i] = H_THERMO

    # Other ocean deep basin
    for i in i_deep:
        mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], deep_mask)
        box_depth[i] = -1

    ## + epicontinental (remove points below seafloor)
    #for i in range(NBASIN):
    #    mask[i,:,:,:] = np.logical_and(mask[i,:,:,:], bottom_mask)

    return mask, i_surf, deepbox_mask, surfbox_mask, thrmbox_mask, polebox_mask, epicbox_mask, appcbox_mask, box_depth



# ----- #


def generate_basin_mask(nav_lat, nav_z, nav_bathy,
                        mask_2D_field=None, epicmask_2D_field=None, mask_values=None, epicmask_values=None,
                        depth_split=(H_SURF, H_THERMO), lat_split=MIDLAT_RANGE, epicont_depth=H_EPICONT,
                        split_epicont=False,
                        polar_auto_identification=False, polar_mask_values=None, polar_epicmask_values=None,
                        merge_polar_upper_boxes=True):
    '''
    Create a mask telling if points are inside (T) or outside (F) each basin
    (GEOCLIM box). Basins are defined accordingly to the depth cutting points,
    the latitude cutting points, the optional given basin masks, and the
    epicontinental boxes.
    The mask is built on the grid defined by the latitude array "nav_lat", the
    vertical column depth array "nav_z", and seafloor depth (bathymetry) array
    "nav_bathy" (longitude array is not needed since there is no longitudinal
    split).

    INPUTS ARGUMENTS:

        nav_lat: numpy-array (rank-3)
            Array indicating the latitude of each grid point

        nav_z: numpy-array (rank-3)
            Array indicating the depth of each grid point

        nav_bathy: numpy-array (rank-3)
            Array indicating the seafloor depth of each grid point.
            A degenerated (size-1) dimension on the vertical axis is expected.

        mask_2D_field: numpy-array (rank-3), or None (default)
            Optional horizontal "mask" array, with 1 value per ocean basin to
            define as GEOCLIM box.
            A degenerated (size-1) dimension on the vertical axis is expected.

        epicmask_2D_field: numpy-array (rank-3), or None (default)
            Same as 'mask_2D_field' for epicontinental boxes.

        mask_values: list of values, or None (default)
            The list of values to consider in "mask_2D_fields".

        epicmask_values: list of values, or None (default)
            Same as 'mask_values' for epicontinental boxes.
        
        depth_split: list of floats, default: (H_SURF, H_THERMO)
            List of cutting depths for defining GEOCLIM boxes.

        lat_split: scalar or list of floats, default: MIDLAT_RANGE
            List of cutting latitudes for defining GEOCLIM boxes.
            If scalar, assume symmetric splitting (same latitude in North and
            South).
            If empty list, False or None, no division will be performed (except
            the ones potentially provided with "mask_2D_fields").
            Note that if "mask_2D_field" is provided, the latitude division
            will be performed FOR EACH "mask_2D_field" values, while excluding
            empty basin.

        epicont_depth: float, default: H_EPICONT
            Seafloor depth defining the limit of epicontinental basins
            (bathymetry < epicont_depth).

        split_epicont: bool, default: False
            True: epicontinental basins will be separated according to
            the horizontal division of open-ocean basins ("lat_split" and
            "mask_2D_field"). False: single epicontinental basin.

        polar_auto_identification: bool, default: False,
            Whether or not automatically identifying polar basins from their
            latitude range.

        polar_mask_values: list of values, or None (default)
            The list of values in "mask_2D_fields" that will be considered as
            polar basins.

        polar_epicmask_values: list of values, or None (default)
            Same as "polar_mask_values" for "epicmask_2D_field".

        merge_polar_upper_boxes: bool, default: True
            Whether or not merging the first boxes of polar water columns.

    RETURN:

        mask: 4-D boolean numpy-array
            mask of all defined GEOCLIM boxes. First dimension is "box",
            remaining ones are 3D ocean grid dimensions.

        i_surf: 1-D integer numpy-array
            Indices of all surface boxes (used for intern computations in
            'oce_output_2_GEOCLIM_BC')

        deep_boxmask: 1-D integer numpy-array
            Whether or not each box is "deep" (1|0)

        surf_boxmask: 1-D integer numpy-array
            Whether or not each box is "surface" (1|0)

        intm_boxmask: 1-D integer numpy-array
            Whether or not each box is "intermediate" ("thermocline") (1|0)

        pole_boxmask: 1-D integer numpy-array
            Whether or not each box is "polar" (1|0)

        epic_boxmask: 1-D integer numpy-array
            Whether or not each box is "epicontinental" (1|0)

        appc_boxmask: 1-D integer numpy-array
            Whether or not each box "receives continental inputs" (1|0)

        box_depth: 1-D float numpy-array
            (bottom) depth of each box 

    NOTES:

    nav_lat, nav_z and nav_bathy, mask_2D_field and epicmask_2D_field MUST ALL
    BE RANK-3, even if they are not defined on all dimensions (e.g., latitude
    and bathymetry should be defined on horizontal, dimensions, whereas z should
    be defined on the vertical dimension).

    mask_2D_field and epicmask_2D_field are optional, but if used, both must be
    given.
    '''


    # ===================== #
    # Check given arguments #
    # ===================== #

    if (mask_2D_field is None) or (epicmask_2D_field is None):

        for var,name in zip([polar_mask_values,   mask_values,   polar_epicmask_values,   epicmask_values],
                           ['polar_mask_values', 'mask_values', 'polar_epicmask_values', 'epicmask_values']):
            if var is not None:
                print('argument "{:}" passed to the function will be ignored: no mask field provided.'.format(name))

        # Trick to cover the whole ocean
        mask_2D_field = np.ones((1,1,1), dtype='int8')
        mask_values = (1,)
        epicmask_2D_field = mask_2D_field
        epicmask_values = mask_values

    elif (mask_2D_field is not None) and (epicmask_2D_field is not None):

        if mask_values is None:
            mask_values = np.unique(mask_2D_field)

        if epicmask_values is None:
            epicmask_values = np.unique(epicmask_2D_field)

        if polar_mask_values is not None and polar_auto_identification:
            print('Automatic identification of polar boxes ("polar_auto_identification=True") desactivated \
                   since the list of values of polar boxes in the mask field ("polar_mask_values") was provided')

    else:
        raise ValueError('Both "mask_2D_field" and "epicmask_2D_field" must be provided, or none.')

    if lat_split is None or lat_split is False or lat_split==() or lat_split==[]:
        lat_split = None
    else:
        # "lat_split" can be a scalar
        if np.size(lat_split) == 1:
            try:
                _, = lat_split
            except TypeError:
                lat_split = (lat_split,)


    # ============== #
    # Initialization #
    # ============== #

    # output variables:
    #  - indices of surface and app_cont boxes
    i_surf = []
    #  - 1/0 indices for each box
    deep_boxmask = []
    surf_boxmask = []
    intm_boxmask = []
    pole_boxmask = []
    epic_boxmask = []
    appc_boxmask = []
    #  - box depth
    box_depth = []
    #  - 3D mask boolean array for all boxes
    mask = []

    # data shape
    shp = np.maximum(np.maximum(nav_lat.shape, nav_z.shape), nav_bathy.shape)

    # mask fields for latitudinal discretization
    if lat_split is None:
        latit_ispolar = cycle([False])
        latit_fieldmask = [True]
    else:
        latit_fieldmask = []
        y0 = -1e12
        for y1 in np.concatenate((lat_split, [1e12])):
            latit_fieldmask.append(np.logical_and(nav_lat>=y0, nav_lat<y1))
            y0 = y1

        if not polar_auto_identification:
            latit_ispolar = [True] + (len(latit_fieldmask)-2)*[False] + [True]
        else:
            latit_ispolar = cycle([None])

    # box ordering from North to South
    latit_fieldmask.reverse()

    # mask fields for vertical discretization
    if depth_split is None:
        depth_fieldmask = [True]
    else:
        depth_fieldmask = []
        z0 = -1e12
        for z1 in np.concatenate((depth_split, [1e12])):
            depth_fieldmask.append(np.logical_and(nav_z>=z0, nav_z<z1))
            z0 = z1

    # Coastal (epicontinental) mask
    epic_maskfield = (nav_bathy <= epicont_depth)


    # ============================================== #
    # Nested loops to identify each individual boxes #
    # ============================================== #

    ibasin = 0
    mask.append(np.zeros(shp, dtype=bool))


    # * non-epicontinental nested loop
    #   ==============================

    # loop on 2D basin mask
    for imsk in mask_values:

        # Loop on latitude
        for latmsk,ispol in zip(latit_fieldmask, latit_ispolar):

            # Loop on depth
            merging = False
            for iz,zmsk in enumerate(depth_fieldmask):

                locmsk = np.logical_and(~epic_maskfield,
                                        np.logical_and(np.logical_and(zmsk, latmsk),
                                                       (mask_2D_field==imsk)))
                if locmsk.any():

                    # Record GEOCLIM basin
                    mask[ibasin] = np.logical_or(mask[ibasin], locmsk)

                    if merging:
                        # Box characteristics are already recorded
                        merging = False # => end of box merging process

                    else:
                        if polar_auto_identification:
                            #
                            # Consider box is polar if >50% of the box is poleward of 60°N/S 
                            #
                            hl_locmsk = locmsk.any(0) # Any point on the vertical column
                            try:
                                hl_locmsk[np.abs(nav_lat[0,:,:])<60] = False
                            except IndexError:
                                hl_locmsk[np.abs(nav_lat[0,:,0])<60, :] = False

                            ispolar = (np.count_nonzero(hl_locmsk) / np.count_nonzero(locmsk.any(0)) > 0.5)

                        elif polar_mask_values is not None:
                            ispolar = (imsk in polar_mask_values)

                        else:
                            ispolar = ispol

                        # record box characteristics:
                        # ---------------------------

                        # - surface flag
                        if iz == 0:
                            surf_boxmask.append(1)
                            i_surf.append(ibasin)
                            if ispolar and merge_polar_upper_boxes and len(depth_fieldmask) >= 2:
                                merging = True # => initiate box merging process
                        else:
                            surf_boxmask.append(0)

                        # - deep and intermediate (thermo) flags
                        if iz == len(depth_fieldmask)-1 or (merging and iz==len(depth_fieldmask)-2):
                            deep_boxmask.append(1)
                            intm_boxmask.append(0)
                        else:
                            deep_boxmask.append(0)
                            if iz > 0:
                                intm_boxmask.append(1)
                            else:
                                intm_boxmask.append(0)

                        # - polar flag
                        pole_boxmask.append(1 if ispolar else 0)

                        # - epicontinental and contiental input flags
                        epic_boxmask.append(0)
                        appc_boxmask.append(0)

                        # - box depth
                        if deep_boxmask[ibasin] == 1:
                            box_depth.append(-1.)
                        else:
                            if merging:
                                box_depth.append(depth_split[iz+1])
                            else:
                                box_depth.append(depth_split[iz])


                    # iterate basin number
                    if not merging:
                        ibasin += 1
                        mask.append(np.zeros(shp, dtype=bool))


    # * epicontinental nested loop
    #   ==========================

    if not split_epicont:
        # erase latitude and 2D mask data
        latit_fieldmask = [True]
        epicmask_2D_field = 1
        epicmask_values = [1]

    if depth_split is None:
        izmax = 1
    else:
        izmax = 1 + np.searchsorted(depth_split, epicont_depth)

    # loop on 2D basin mask
    for imsk in epicmask_values:

        # Loop on latitude
        for latmsk in latit_fieldmask:

            # Loop on depth
            for iz,zmsk in enumerate(depth_fieldmask[:izmax]):

                locmsk = np.logical_and(epic_maskfield,
                                        np.logical_and(np.logical_and(zmsk, latmsk),
                                                       (epicmask_2D_field==imsk)))
                if locmsk.any():

                    if polar_epicmask_values is not None:
                        ispolar = (imsk in polar_epicmask_values)
                    else:
                        ispolar = False
                        # -> no epicontinental polar boxes, unless an epicontinental mask field is provided.

                    # Record GEOCLIM basin
                    mask[ibasin] = np.logical_or(mask[ibasin], locmsk)

                    # record box characteristics:

                    # - surface flag
                    if iz == 0:
                        surf_boxmask.append(1)
                        i_surf.append(ibasin)
                    else:
                        surf_boxmask.append(0)

                    # - deep and intermediate (thermo) flags
                    deep_boxmask.append(0)
                    intm_boxmask.append(0)

                    # - polar flag
                    pole_boxmask.append(1 if ispolar else 0)

                    # - epicontinental flag
                    epic_boxmask.append(1)

                    # - continental input flags
                    appc_boxmask.append(1 if surf_boxmask[ibasin]==1 else 0)

                    # - box depth
                    box_depth.append(-1.)


                    # iterate basin number
                    ibasin += 1
                    mask.append(np.zeros(shp, dtype=bool))


    # Remove extra mask field
    del(mask[-1])


    # Add atmospheric box
    deep_boxmask.append(0)
    surf_boxmask.append(0)
    intm_boxmask.append(0)
    pole_boxmask.append(0)
    epic_boxmask.append(0)
    appc_boxmask.append(0)
    box_depth.append(0)


    return np.array(mask), np.array(i_surf, dtype='int8'), \
           np.array(deep_boxmask, dtype='int8'), np.array(surf_boxmask, dtype='int8'), np.array(intm_boxmask, dtype='int8'), \
           np.array(pole_boxmask, dtype='int8'), np.array(epic_boxmask, dtype='int8'), np.array(appc_boxmask, dtype='int8'), np.array(box_depth)



#######################################################################
##-------------------------------------------------------------------##
##-- MAIN FUNCTION: READ NETCDF FILES AND WRITE GEOCLIM INPUT FILE --##
##-------------------------------------------------------------------##
#######################################################################


def oce_output_2_GEOCLIM_BC(input_files, latitude, z, temperature,
                            bathy=None, time_dim=None, root='', co2=None,
                            x_weight=None, y_weight=None, horiz_weight=None, z_weight=None, z_bounds=None, weight=None,
                            grid_file=None, input_u_files=None, input_v_files=None, input_w_files=None,
                            u=None, v=None, w=None, w_method='divergence', inverted_w=False,
                            flux_correction='relative least square', h_mix_reduct_factor=1.,
                            periodic_x=True, x_overlap=0, special_wrap=None,
                            basins_kw='historical',
                            mask_nc_file=None, mask_nc_var=None, epicmask_nc_var=None, mask_values=None, epicmask_values=None,
                            outdir='./',
                            temp_outfile='GCM_oceanic_temp.dat', surf_outfile='oce_surf.dat', sedsurf_outfile='surf_sedi.dat',
                            vol_outfile='oce_vol.dat', watflux_outfile='exchange_2.dat', bndlen_outfile='box_sf_bnd_len.dat',
                            deepbox_outfile='indice_deep.dat', surfbox_outfile='indice_surface.dat',
                            thermobox_outfile='thermocline.dat', polarbox_outfile='indice_polar.dat',
                            epicontbox_outfile='indice_epicont.dat', appcontbox_outfile='apport_ct.dat',
                            depth_press_file='press_box.dat',  mask_outfile='basins_masks.nc',
                            export_mask_fields=False, check_plot=False):
    '''
    This function read the oceanic output and grid definition of a GCM, compute
    the volumes, surfaces and basin-average temperature of GEOCLIM basins, and
    save them in 4 GEOCLIM-readable ascii files, that can be used directly as
    GEOCLIM input: temperature file, box surface file, box sedimentary surface
    file, and box volume file.

    MANDATORY INPUT ARGUMENTS:

        input_files: list of string
            list of names of GCM output netCDF files, each corresponding to one
            "climate state" (combination of CO2 and climatic parameters).
                    
        latitude: string
            Name of the variable giving the latitude of each point. That
            variable should be either in the main input files, or in the
            (optional) grid file (see argument "grid_file"). It can be 1D, or
            2D (if the horizontal grid is not longitude-latitude, for instance).

        z: string
            Name the variable giving the depth of each point of the 3D grid.
            That variable should be either in the main input files, or in the
            (optional) grid file (see argument "grid_file"). It must be 1D.

        temperature: string
            Name the variable giving the temperature of each point of the 3D
            grid. That variable should be either in all the main input climate
            files Must be 3D or 4D if time-dependent. In that case, it will be
            averaged along the time dimension.

    OPTIONAL INPUT ARGUMENTS:

        bathy: string
            Name the variable giving, at each "horizontal" grid point, the depth
            of ocean floor. That variable should be either in the main input
            files, or in the (optional) grid file (see argument "grid_file").
            It must be 2D (horizontal). If not provided, the seafloor depth will
            be computed according to the valid points of temperature field.

        time_dim: string
            Name of the time dimension if temperature needs to be averaged on
            time. Note that the program will by default try to identify it using
            common time dimension names.

        root: string
            Root directory where to the input files are (or more generally,
            common part to all file paths)

        co2: list of floats
            Archaism. CO2 levels corresponding to each item of the list(s) of
            climate files

        x_weight: string
            Name of the variable giving the weighting of the grid in the x
            direction (ie, width). Must be 1D or 2D.

        y_weight: string
            Name of the variable giving the weighting of the grid in the y
            direction (ie, width). Must be 1D or 2D.

        horiz_weight: string
            Name of the variable giving the weighting of the grid in the
            horizontal directions (x and y. ie, area). Must be 2D.

        z_weight: string
            Name of the variable giving the weighting of the grid in the z
            direction (ie, thickness). Must be 1D.

        z_bounds: string
            Alternative to "z_weight", name of the variable giving the upper
            and lower bounds of the vertical levels. z-weighting will then be
            computed as diff(z_bounds). Must be 2D, the 2nd dimension being the
            bounds (size-2).

        weight: string
            Name of the variable giving the 3-D weighting of each point of the
            grid. Must be 3D.

        grid_file: string
            Path of an alternative file to get the grid variables (latitude, z,
            bathy and weightings)

        input_u_files, input_v_files, input_w_files : lists of string
            List of names (paths) of GCM output netCDF files for U, V and W
            variables (respectively), IN SAME ORDER THAN "input_files" ARGUMENT.

        u, v, w: string
            Name of the netCDF variable for oceanic U/V/W velocity (or flux)
            field (respectively)

        w_method: string,
           'data' or 'divergence' (default)
           How to compute the vertical fluxes: from loaded data W ("data") or
           as the divergence of U and V ("divergence")

        inverted_w: bool, default: False
           Whether the orientation of positive W values is inverted with respect
           to Z axis

        flux_correction: string or None
            How to compute the correction of the matrix of water fluxes between
            boxes X to ensure mass conservation (=> X.sum(0) == X.sum(1)).
                None/'none':
                    No flux correction
                'pls'|'positive least square'|'positive least-square':
                    Additive correction dX (X => X+dX), computed by minimizing
                    the sum of dF^2, and modified afterward to avoid negative
                    fluxes in X+dX.
                    This method is unconditionally stable, but can generate
                    fluxes between unconnected boxes.
                'rls'|'relative least square'/relative least-square':
                    Default method.
                    Multiplicative correction F (X => X*F, element-wise
                    multiplication), computed by minimizing the sum of (1-F)^2.
                    This method has better properties: correction computed
                    relatively to absolute fluxes, keep 0 fluxes between
                    unconnected boxes. However, it is not always stable (it
                    requires a numerical matrix inversion), and is known not
                    to work if they are cluster of boxes totally disconnected.

        h_mix_reduct_factor: float, default: 1.
            Reduction factor applied to bidirectional horizontal water fluxes.
            Note: fluxes are DIVIDED by 'h_mix_reduct_factor', not multiplied by

        periodic_x: bool, default: True
            Whether the "x" axis is handle as periodic (for water fluxes)

        x_overlap: integer, default: 0
            How much "columns" (in the X direction) are duplicated, and will NOT

        special_wrap: string or None (default)
            Specificity of the North pole fold of the IPSL "ORCA" grid.
            Indicates the code how to connect to northn cells, and which ones
            not to consider. Legal values are:
                None: no North Pole fold
                'ORCA-T'|'ORCA-Tpivot'|'ORCA T-pivot': NP fold on a "T" point
                'ORCA-F'|'ORCA-Fpivot'|'ORCA F-pivot': NP fold on a "F" point
            NOTE: whatever the size "ny" of the grid, the pivot point is
            expected to be at 'ny//2'.

        basins_kw: dictionary, or string "historical"
            If "historical", basins are defined as the historical GEOCLIM
            9 ocean + 1 atmosphere boxes.
            Otherwise, expect a dictionay of key-word arguments that will be
            passed to the function 'generate_basin_mask'. Keys and items are:
                'depth_split':
                    tuple of cutting depths values.
                'lat_split':
                    scalar cutting latitude value.
                'epicont_depth':
                    seafloor depth defining epicontinental basins.
                'split_epicont':
                    True/False, single epicontinental basin, or split
                    accordingly to the provided basin mask.
                    This argument is ignored if an epicontinental mask field is
                    provided (see "epicmask_nc_var" argument).
                'polar_auto_identification':
                    True/False, to automatically identify polar boxes from their
                    latitude range.
                'polar_mask_values':
                    None (if automatic polar identification) or tuple of "mask"
                    values corresponding to polar boxes (if a "mask" is given
                    to the main function, see main arguments "mask_nc_file",
                    "mask_nc_var" and "mask_values").
                'polar_epicmask_values':
                    Same as 'polar_mask_values' for epicontinental boxes.
                    This input is only used if an epicontinental mask field is
                    given. Note thatn in the default settings of GEOCLIM, none
                    of the epicontinental boxes have a "polar" flag.
                'merge_polar_upper_boxes':
                    True/False, to merge the upper two boxes of polar water
                    columns (in open ocean only).
            See "generate_basin_mask" function for more details, and default
            values.

        mask_nc_file, mask_nc_var, epicmask_nc_var: string
            Those arguments are only used if basins_kw!='historical'.
            mask_nc_file: path of the netCDF file containing the variables:
            mask_nc_var:     name of the netCDF variable of the basins mask.
                             This variable must be 2D, and have a small number
                             of values, each value corresponding to the
                             horizontal division of a basin.
            epicmask_nc_var: name of the (optional) netCDF variable of the
                             epicontinental basins mask. Similar to mask_nc_var.
            Open ocean basins are always the areas from "mask_nc_var" deeper
            than 'epicont_depth' (see "basins_kw" argument).
            Epicontinental basins are the areas from "mask_nc_var" shallower
            than 'epicont_depth' if "epicmask_nc_var" is None, otherwise, the
            areas from "epicmask_nc_var" shallower than 'epicont_depth'.

        mask_values: list of numbers
            List of values of the basin "mask" variable ("mask_nc_var" argument)
            that will be considered for the definition of GEOCLIM boxes.

        epicmask_values: list of numbers
            List of values of the epicontinental basin "mask" variable
            ("epicmask_nc_var" argument) that will be considered for the
            definition of GEOCLIM epicontinental boxes.

        outdir: string
            Path of the directory where the output files will be written.

        temp_outfile, surf_outfile, sedsurf_outfile, vol_outfile,
        watflux_outfile, bndlen_outfile, deepbox_outfile, surfbox_outfile,
        thermobox_outfile, polarbox_outfile, epicontbox_outfile,
        appcontbox_outfile, depth_press_file,  mask_outfile:
            string
            Name (ie, path) of the output files for (respectively):
            box mean temperature, box top area, box sedimenting (seafloor)
            area, box volume, water exchange between boxes, length of horizontal
            boundaries between boxes, "is deep" box indices, "is surface" box
            indices, "is intermediate" box indices, "is polar" box indices, "is
            epicontinental" box indices, "receives continental flux" box
            indices, box depth and pressure, netCDF box mask file.
            See function definition for default names.

        export_mask_fields: bool
            Whether or not exporting basins definition on the native ocean grid
            as a mask in a netCDF file ("mask_outfile")

        check_plot: bool
            Whether or not plotting basin definition and characteristic during
            their computation.

        NOTES:

        The priority order of the variables for the grid cell weighting (in
        case of conflict) is the following:
          weight > z_weight,horiz_weight > x_weight,y_weight
        In case of missing information for the weighting in one or several
        dimensions, A UNIFORM WEIGHTING WILL BE APPLIED.

        Whether U/V/W variables are water velocity or fluxes will be determined
        by identifying the physical units of the variables.

        Output pressure and depth units are "bar" and "km"
        Output volume unit is "Mkm3" ("1e15 m3")
        Output area unit is "Gkm2" ("1e15 m2").
        Output temperature units is "°C".
        Output water flux units is "Sv".
    '''


    # Chek-plotting module
    if check_plot:
        try:
            from matplotlib import pyplot as plt
        except ModuleNotFoundError:
            print('')
            print('Warning: Matplotlib module not found, cannot draw check-plots')
            check_plot = False


    # Output directory
    if outdir[-1]!='/':
        outdir = outdir+'/'

    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        print('Created directory "'+outdir+'"')



    # Check optional arguments
    # ------------------------

    if basins_kw != 'historical' and not isinstance(basins_kw, dict):
        raise TypeError('"basins_kw" argument must either be "historical" or be a dictionary of keyword-arguments passed to the function "generate_basin_mask"')



    # ====================== #
    # Get inpput information #
    # ====================== #


    # CO2 array (for GEOCLIM without climatic parameters)
    # ---------------------------------------------------

    if co2 is None:
        nclim = len(input_files)
        co2 = cycle([-0.12345])
    else:
        nclim = len(co2)


    # Number of input files
    # ---------------------

    if len(input_files) != nclim:
        raise ValueError('Number of input files inconsistent with CO2 axis')

    if input_u_files is None:
        input_u_files = cycle([None])
    elif len(input_u_files) != nclim:
        raise ValueError('Number of input "u" files inconsistent with number of input "Temp" files')

    if input_v_files is None:
        input_v_files = cycle([None])
    elif len(input_v_files) != nclim:
        raise ValueError('Number of input "v" files inconsistent with number of input "Temp" files')

    if input_w_files is None:
        input_w_files = cycle([None])
    elif len(input_w_files) != nclim:
        raise ValueError('Number of input "w" files inconsistent with number of input "Temp" files')

    all_inputs = (input_files, input_u_files, input_v_files, input_w_files)


    # Grid information:
    # -----------------

    # Grid file (if not provided, get variables in the 1st main files)
    if grid_file is None:
        fgrid = None
    else:
        try:
            fgrid  = nc.Dataset(root+grid_file)
        except FileNotFoundError as err:
            fgrid  = nc.Dataset(grid_file)

    fgrid2 = nc.Dataset(root+input_files[0])

    xwg = get_var(fgrid, fgrid2, varname=x_weight)
    ywg = get_var(fgrid, fgrid2, varname=y_weight)
    zwg = get_var(fgrid, fgrid2, varname=z_weight)
    zbd = get_var(fgrid, fgrid2, varname=z_bounds)
    hwg = get_var(fgrid, fgrid2, varname=horiz_weight)
    wgh = get_var(fgrid, fgrid2, varname=weight)

    # Main variables (try the main file before grid file)
    lati = get_var(fgrid2, fgrid, varname=latitude, raise_error=True)
    zvar = get_var(fgrid2, fgrid, varname=z, raise_error=True)
    baty = get_var(fgrid2, fgrid, varname=bathy)

    # temperature of the 1st main file (to check existence, shape and dimensions)
    temp = fgrid2[temperature]

    # Check units (perform dummy conversion)
    _ = temperature_units.convert(0, temp.units)

    # Dimension checks:
    if temp.ndim == 3:
        idx_t = () # trick: empty tuple, so that temp.mean(idx_t) does nothing.
        temp_shape = temp.shape
        temp_dims = temp.dimensions
    elif temp.ndim == 4:
        if time_dim is None:
            idx_t = None
            for tname in ['time', 'time_counter', 'month']:
                if tname in temp.dimensions:
                    time_dim = tname
                    break

            if time_dim is None:
                raise ValueError('Cannot identify time dimension for temperature variable.')

        idx_t = temp.dimensions.index(time_dim)
        temp_shape = list(temp.shape)
        del(temp_shape[idx_t])
        temp_shape = tuple(temp_shape)
        temp_dims = list(temp.dimensions)
        del(temp_dims[idx_t])
        temp_dims = tuple(temp_dims)
                
    else:
        raise ValueError('temperature variable (argument "temperature") must be rank-3, or rank-4 if time-dependent.')


    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#
    #@ Expected dimensions (time excluded): @#
    #@  0: z                                @#
    #@  1: y                                @#
    #@  2: x                                @#
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@#


    # Special grid boundary conditions (x-overlapping, ...)
    # -----------------------------------------------------

    wrap_pivot_point = None

    if special_wrap in ('ORCA-T', 'ORCA-Tpivot', 'ORCA T-pivot', 'ORCA-F', 'ORCA-Fpivot', 'ORCA F-pivot'):
        ix0 = 1
        ix1 = temp_shape[2]-1
        x_overlap = 2
        #
        wrap_pivot_point = special_wrap[5]
        special_wrap = 'ORCA'

    else:
        ix0 = 0
        ix1 = temp_shape[2] - x_overlap
        if special_wrap is not None:
            print('WARNING: unknown special wrapping case "'+str(special_wrap)+'".')
            print('         Will be ignored.')
            special_wrap = None


    horiz_shp = temp_shape[1:]


    # create "nav_z"
    # --------------
    shp = list(temp_shape)
    if zvar.ndim == 3:
        nav_z = zvar[:,:,:]

    elif zvar.ndim == 2:
        if zvar.dimensions[1] == temp_dims[2]: # "z" dimensions are {x,z}
            shp[1] = 1 # remove y dimension
        elif zvar.dimensions[1] == temp_dims[1]: # "z" dimensions are {y,z}
            shp[2] = 1 # remove x dimension
        else:
            raise ValueError('Cannot identify dimensions of "z" variable (do not match "temperature" variable)')

        nav_z = zvar[:,:].reshape(shp)

    elif zvar.ndim == 1:
        shp[1] = 1 # remove y dimension
        shp[2] = 1 # remove x dimension
        nav_z = zvar[:].reshape(shp)

    else:
        raise ValueError('Illegal number of dimensions for variable "z"')

    # z must be positive
    if (nav_z<=0).all():
        nav_z = -nav_z

    # check z ordering:
    reverse_z = (nav_z[0,:,:] > nav_z[1,:,:]).any()
    isurf = -1 if reverse_z else 0

    # Units conversion
    nav_z = length_units.convert(nav_z, zvar.units)

    # create "nav_lat"
    # ----------------
    shp = list(temp_shape)
    if lati.ndim == 3:
        raise ValueError('Cannot handle 3-dimensionnal latitude field. Must be 1D or 2D')
    elif lati.ndim == 2:
        shp[0] = 1 # no z dimension
    elif lati.ndim == 1:
        shp[0] = 1 # no z
        shp[2] = 1 # nor x dimension
    else:
        raise ValueError('Illegal number of dimensions for variable "latitude"')

    nav_lat = latitude_units.convert(lati[:], lati.units).reshape(shp)



    # ======================================================== #
    # Volumetric (3-D) and/or area (horizontal) grid weighting #
    # ======================================================== #


    # Volumetric
    # ----------
    if wgh is not None:
        if wgh.ndim != 3:
            raise ValueError('weight variable (argument "weight") must be rank-3')
        else:
            if wgh.shape != temp_shape:
                raise ValueError('temperature and weight variables (arguments "temperature" and "weight") must have the same shape')
            else:
                if temp_dims != wgh.dimensions: 
                    print('WARNING: temperature and weight variables dimensions does not have the same names')

                w_vol = wgh[:,:,:]
                try:
                    w_vol = volume_units.convert(w_vol, wgh.units)
                    volunits = True
                except (AttributeError, UnknownUnitError):
                    volunits = False

    else:
        w_vol = None

    # Area
    # ----
    if hwg is None:
        w_horiz = None
    else:
        if hwg.ndim != 2:
            raise ValueError('horizontal weight variable (argument "horiz_weight") must be rank-2')
        else:
            for d in lati.dimensions:
                if d not in hwg.dimensions:
                    print('WARNING: latitude and horizontal-weight variables dimensions does not have the same names')

            w_horiz = hwg[:,:]
            try:
                w_horiz = area_units.convert(w_horiz, hwg.units)
                xyunits = True
            except (AttributeError, UnknownUnitError):
                xyunits = False


    # x, y and z axis weighting
    # -------------------------

    if xwg is None and ywg is None:
        if w_horiz is None:
            w_horiz = np.ones(horiz_shp)
            xyunits = False
            if w_vol is None:
                print('WARNING: uniform weighting assumed in horizontal dimensions (x and y)')

    # x weighting
    w_x, xunits = get_axis_weighting('x', xwg, temp_dims, temp_shape, t_dim=idx_t)

    # y weighting
    w_y, yunits = get_axis_weighting('y', ywg, temp_dims, temp_shape, t_dim=idx_t)

    # x weighting
    w_z, zunits = get_axis_weighting('z', zwg, temp_dims, temp_shape, axis_bnds_var=zbd, axis_var=zvar, t_dim=idx_t)


    # Compute areal (horizontal) weighting if it wasn't loaded
    if w_horiz is None:
        if xwg is None and ywg is None:
            xyunits = False
            if w_vol is None:
                print('WARNING: uniform weighting assumed in horizontal dimensions (x and y)')
            else:
                # consider first "surface" volumetric weighting
                w_horiz = w_vol[isurf,:,:]
                
        else:
            w_horiz = w_x[isurf,:,:] * w_y[isurf,:,:]
            xyunits = (xunits and yunits)

        if not xyunits:
            print('WARNING: x and y units not understood. Areas will be computed as a fraction of modern oceanic areas')


    # Compute volumetric weighting if it wasn't loaded
    if w_vol is None:
        w_vol = w_horiz.reshape((1,)+horiz_shp) * w_z
        volunits = (xyunits and zunits)
        if not volunits:
            print('WARNING: horizontal or vertical units not understood. Volume will be computed as a fraction of modern oceanic areas')


    # Special domain definition
    # + + + + + + + + + + + + +
    if special_wrap == 'ORCA':
        # => irregular inner domain mask at North fold boundary in ORCA grid
        remove_orca_north_fold(w_horiz, pivot=wrap_pivot_point)
        remove_orca_north_fold(w_vol,   pivot=wrap_pivot_point)
        remove_orca_north_fold(w_x,     pivot=wrap_pivot_point)
        remove_orca_north_fold(w_y,     pivot=wrap_pivot_point)
        remove_orca_north_fold(w_z,     pivot=wrap_pivot_point)



    # Get depth and create "nav_bathy" (that may be computed from z weighting)
    # ========================================================================

    if bathy is None:
        # bathymetry is computed from cell thickness, if possible, rather than determined from cell depth
        if zwg is not None and zunits:
            baty, bottom_idx = find_bathy(var_mask=temp[:].mask.all(idx_t), cell_thick=w_z, invert_zaxis=reverse_z)
            print('NOTE: "bathy" field (bathymetry) computed from cell thickness, and determined')
            print('      from the last valid point of temperature field on each vertical column')
        else:
            baty, bottom_idx = find_bathy(var_mask=temp[:].mask.all(idx_t), z=nav_z, invert_zaxis=reverse_z)
            print('NOTE: "bathy" field (bathymetry) automatically computed as the z-value of')
            print('      the last valid point of temperature field on each vertical column')

    else:
        bottom_idx = find_bathy(var_mask=temp[:].mask.all(idx_t), invert_zaxis=reverse_z)

        if baty.ndim != 2:
            raise ValueError('seafloor depth variable (argument "bathy") must be rank-2')
        else:
            for d in lati.dimensions:
                if d not in baty.dimensions:
                    print('WARNING: latitude and seafloor depth variables dimensions does not have the same names')
            
            baty = length_units.convert(baty[:], baty.units)

    shp = list(temp_shape)
    shp[0] = 1 # no z dimension
    nav_bathy = baty.reshape(shp)



    # ================================= #
    # Get mask of each basin of GEOCLIM #
    # ================================= #


    if basins_kw == 'historical':
        Gmask, i_surfbox, deep_mask, surf_mask, intm_mask, pole_mask, epic_mask, appc_mask, box_depth = \
        geoclim_default_basin_mask(nav_lat, nav_z, nav_bathy)
        nbasin = NBASIN

    else:

        if mask_nc_file is None or mask_nc_var is None:
            basins_kw['mask_2D_field'] = None
            basins_kw['mask_values']   = None
        else:
            try:
                fmsk = nc.Dataset(root+mask_nc_file)
            except FileNotFoundError as err:
                fmsk = nc.Dataset(mask_nc_file)
            
            basins_kw['mask_2D_field'] = fmsk[mask_nc_var][:,:].reshape((1,)+horiz_shp)
            basins_kw['mask_values']   = mask_values

            if epicmask_nc_var is None:
                basins_kw['epicmask_2D_field'] = basins_kw['mask_2D_field']
                basins_kw['epicmask_values']   = basins_kw['mask_values']
            else:
                basins_kw['epicmask_2D_field'] = fmsk[epicmask_nc_var][:,:].reshape((1,)+horiz_shp)
                basins_kw['epicmask_values']   = epicmask_values
                basins_kw['split_epicont'] = True

            fmsk.close()

        Gmask, i_surfbox, deep_mask, surf_mask, intm_mask, pole_mask, epic_mask, appc_mask, box_depth = \
          generate_basin_mask(nav_lat, nav_z, nav_bathy, **basins_kw)
        nbasin = deep_mask.size - 1 # do not count atmospheric box 

    # Remove missing-points from basin masks:
    Gmask[:, w_vol==0] = False



    # ======================================================================= #
    # Compute and export temperature average per basin -- loop on input files #
    # ======================================================================= #


    # Initialization
    geoclim_input = np.zeros((nclim,nbasin+1), dtype=float)
    # Note: For retro-compatibility, the first column of each row in oceanic temperature
    # GEOCLIM input file is expected to be CO2 levels. the others columns are basin temperatures

    # Check: draw maps of masks
    # <+++++++++++++++++++++++++++++++++++++++> #
    if check_plot:
        fin = nc.Dataset(root+input_files[0])
        tempmask = fin[temperature][:].mask.all(idx_t)
        fin.close()
        _, ax = plt.subplots()
        pid = ax.pcolormesh(np.ma.masked_where(tempmask[isurf,:,:], nav_bathy[0,:,:]))
        plt.colorbar(pid, label='m')
        ax.set_title('bathymetry')
        for k in range(nbasin):
            locmask = np.logical_and(Gmask[k,:,:,:], ~tempmask)
            print('Basin #{0:}, volume fraction: {1:f}'.format(k, w_vol[locmask].sum()/w_vol.sum()) \
                  + (' POLAR' if pole_mask[k] else ''))
            fig, ax = plt.subplots(nrows=3, figsize=(6,8))
            cax = [fig.add_axes([0.95, 0.9-0.8*(i+1)/3, 0.01, 0.2]) for i in range(3)]
            pid = ax[0].pcolormesh(np.count_nonzero(locmask, axis=0))
            ax[0].set_xlabel('x')
            ax[0].set_ylabel('y')
            plt.colorbar(pid, cax=cax[0])
            pid = ax[1].pcolormesh(np.count_nonzero(locmask, axis=2))
            ax[1].set_xlabel('y')
            ax[1].set_ylabel('z')
            plt.colorbar(pid, cax=cax[1])
            pid = ax[2].pcolormesh(np.count_nonzero(locmask, axis=1))
            ax[2].set_xlabel('x')
            ax[2].set_ylabel('z')
            plt.colorbar(pid, cax=cax[2])
            fig.suptitle('Basin #{0:}'.format(k))
            if not reverse_z:
                ax[1].invert_yaxis()
                ax[2].invert_yaxis()
        plt.show()
    # <+++++++++++++++++++++++++++++++++++++++> #

    for k,(inputfile,c) in enumerate(zip(input_files, co2)):

        fname = root+inputfile

        fin = nc.Dataset(fname)

        # Main variable
        temp = fin[temperature]
        tempvar = temp[:].mean(idx_t)

        # CO2 level
        geoclim_input[k,0] = c

        # basin-average temperature
        for i in range(nbasin):
            geoclim_input[k,i+1] = (tempvar[:,:,ix0:ix1][Gmask[i,:,:,ix0:ix1]] * w_vol[:,:,ix0:ix1][Gmask[i,:,:,ix0:ix1]]).sum() / \
                                   w_vol[:,:,ix0:ix1][np.logical_and(Gmask[i,:,:,ix0:ix1], ~tempvar.mask[:,:,ix0:ix1])].sum()

        geoclim_input[k,1:] = temperature_units.convert(geoclim_input[k,1:], temp.units)

        fin.close()


    # Write in file
    np.savetxt(outdir+temp_outfile, geoclim_input, fmt='%10.5f')



    # =============================== #
    # Compute and export water fluxes #
    # =============================== #

    w_method = w_method.lower()
    if w_method == 'divergence':
        need_w = False
        inverted_w = False
    elif w_method == 'data':
        need_w = True
    else:
        print('WARNING: unkown w computation method "{:}".'.format(w_method))
        print('CANNOT COMPUTE EXCHANGE FLUXES')
        print('')
        w_method = '~#!ERROR!#~'


    if (u is None) or (v is None) or (need_w and w is None):
        print('Missing "u", "v", or "w" fields. Exchange fluxes not computed.')


    elif w_method != '~#!ERROR!#~':


        # Flux correction method
        # ----------------------

        flux_correction = flux_correction.lower()

        if flux_correction is None or flux_correction == 'none':
            flx_message = ''
            def flxcorr(*args, **kwargs):
                pass

        elif flux_correction in ('pls', 'positive least square', 'positive least-square'):
            flx_message = '\n\t--> positive least-square flux correction'
            flxcorr = PositiveLeastSquare_correction

        elif flux_correction in ('rls', 'relative least square', 'relative least-square'):
            flx_message = '\n\t--> relative least-square flux correction'
            flxcorr = RelativeLeastSquare_correction

        else:
            print('')
            print('Warning: unknown flux correction method "{:}".'.format(flux_correction))
            print('Use "relaive least square" method instead')

            # default flux correction method
            flx_message = '\n\t--> relative least-square (default) flux correction'
            flxcorr = RelativeLeastSquare_correction


        fout = open(outdir+watflux_outfile, mode='w')

        # Finite volume approach:
        # compute areas of the surfaces between cells (in the perpendiular direction)
        #    -> "y-z" area between cell x=i and cell x=i+1 is ("y-z" area[cell i]  +  "y-z" area[cell i+1])/2
        # except that for vertical dimension, the minimum "thickness" between cell i and cell i+1 is considered:
        #    vertical cell thickness is expected to be uniform along x and y dimension, except for bottom cells, that intercept
        #    sea-floor with different depths. -> take the minimum thickness between 2 adjacent cell as the exchange surface

        if yunits and zunits: # average/minimum in the x direction
            w_yz = np.zeros(w_z.shape, w_z.dtype)
            w_yz[:,:,ix0:ix1-1] = average(w_y[:,:,ix0:ix1-1], w_y[:,:,ix0+1:ix1]) * minimum(w_z[:,:,ix0:ix1-1], w_z[:,:,ix0+1:ix1])
            if periodic_x:
                w_yz[:,:,ix1-1] = average(w_y[:,:,ix1-1], w_y[:,:,ix0]) * minimum(w_z[:,:,ix1-1], w_z[:,:,ix0])

        if xunits and zunits: # average/minimum in the y direction
            w_xz = np.zeros(w_z.shape, w_z.dtype)
            w_xz[:,:-1,ix0:ix1] = average(w_x[:,:-1,ix0:ix1], w_x[:,1:,ix0:ix1]) * minimum(w_z[:,:-1,ix0:ix1], w_z[:,1:,ix0:ix1])
            if special_wrap == 'ORCA':
                ihalf = w_x.shape[-1]//2
                if wrap_pivot_point == 'T':
                    w_xz[:,-2,ix0+1:ihalf] = average(w_x[:,-2,ix0+1:ihalf], w_x[:,-3,ix1-1:ihalf:-1]) * minimum(w_z[:,-2,ix0+1:ihalf], w_z[:,-3,ix1-1:ihalf:-1])
                else: # => wrap_pivot_point == 'F':
                    w_xz[:,-2,ix0:ihalf] = average(w_x[:,-2,ix0:ihalf], w_x[:,-2,ix1-1:ihalf-1:-1]) * minimum(w_z[:,-2,ix0:ihalf], w_z[:,-2,ix1-1:ihalf-1:-1])

        if xunits and yunits: # average in the z direction
            w_xy = np.zeros(w_z.shape, w_z.dtype)
            if inverted_w:
                w_xy[1:,:,ix0:ix1] = average(w_x[:-1,:,ix0:ix1], w_x[1:,:,ix0:ix1]) * average(w_y[:-1,:,ix0:ix1], w_y[1:,:,ix0:ix1])
            else:
                w_xy[:-1,:,ix0:ix1] = average(w_x[:-1,:,ix0:ix1], w_x[1:,:,ix0:ix1]) * average(w_y[:-1,:,ix0:ix1], w_y[1:,:,ix0:ix1])

        units = lambda flx: velocity_units if flx is None else flux_units


        print('')


        # <------------------------------------------------------------------------------------->
        # <-- Loop on all u,v,w input files => compute exchange matrix for all climate states -->
        # <------------------------------------------------------------------------------------->

        missing_axis_units = False

        for allfiles in zip(*all_inputs):

            allf = [None if fname is None else nc.Dataset(root+fname) for fname in allfiles]

            u_var = get_var(*allf, varname=u, raise_error=True)
            v_var = get_var(*allf, varname=v, raise_error=True)
            if need_w:
                w_var = get_var(*allf, varname=w, raise_error=True)

            try:
                u_field = velocity_units.convert(u_var[:].mean(idx_t), u_var.units)
                uisveloc = True
            except UnknownUnitError:
                u_field = flux_units.convert(u_var[:].mean(idx_t), u_var.units)
                uisveloc = False

            try:
                v_field = velocity_units.convert(v_var[:].mean(idx_t), v_var.units)
                visveloc = True
            except UnknownUnitError:
                v_field = flux_units.convert(v_var[:].mean(idx_t), v_var.units)
                visveloc = False

            if need_w:
                try:
                    w_field = velocity_units.convert(w_var[:].mean(idx_t), w_var.units)
                    wisveloc = True
                except UnknownUnitError:
                    w_field = flux_units.convert(w_var[:].mean(idx_t), w_var.units)
                    wisveloc = False

            elif w_method == 'divergence':
                wisveloc = True


            # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
            # Matrix of fluxes: X_mat[i,j] = flux FROM box i TO box j #
            # + + + + + + + + + + + + + + + + + + + + + + + + + + + + #
            X_mat = np.zeros((nbasin+1,nbasin+1), dtype=float)

            msk_uij, msk_vij, msk_wij = np.zeros((3,)+Gmask.shape[1:], dtype=bool)

            for i in range(nbasin):
                for j in list(range(i))+list(range(i+1,nbasin)):

                    uij, vij, wij = (), (), ()
                    msk_uij[:], msk_vij[:], msk_wij[:] = False, False, False

                    # u-fluxes
                    msk_uij[:,:,ix0:ix1-1] = Gmask[i,:,:,ix0:ix1-1] & Gmask[j,:,:,ix0+1:ix1] & ~u_field[:,:,ix0:ix1-1].mask
                    if periodic_x:
                        msk_uij[:,:,ix1-1] = Gmask[i,:,:,ix1-1] & Gmask[j,:,:,ix0] & ~u_field[:,:,ix1-1].mask

                    if msk_uij.any():
                        uij = u_field[msk_uij].data
                        if uisveloc:
                            if yunits and zunits:
                                uij *= w_yz[msk_uij]
                            else:
                                print('Warning: cannot compute water fluxes from U velocity field: missing units for y or z axis weigthing')
                                missing_axis_units = True

                        pos = (uij >= 0)
                        X_mat[i,j] += uij[pos].sum()
                        X_mat[j,i] -= uij[~pos].sum()

                    # v-fluxes
                    msk_vij[:,:-1,ix0:ix1] = Gmask[i,:,:-1,ix0:ix1] & Gmask[j,:,1:,ix0:ix1] & ~v_field[:,:-1,ix0:ix1].mask
                    if special_wrap == 'ORCA':
                        ihalf = w_vol.shape[2]//2
                        if wrap_pivot_point == 'T':
                            msk_vij[:,-2,ix0+1:ihalf] = Gmask[i,:,-2,ix0+1:ihalf] & Gmask[j,:,-3,ix1-1:ihalf:-1] & ~v_field[:,-2,ix0+1:ihalf].mask
                        else: # => wrap_pivot_point == 'F'
                            msk_vij[:,-2,ix0:ihalf] = Gmask[i,:,-2,ix0:ihalf] & Gmask[j,:,-2,ix1-1:ihalf-1:-1] & ~v_field[:,-2,ix0:ihalf].mask

                    if msk_vij.any():
                        vij = v_field[msk_vij].data
                        if visveloc:
                            if xunits and zunits:
                                vij *= w_xz[msk_vij]
                            else:
                                print('Warning: cannot compute water fluxes from U velocity field: missing units for x or z axis weigthing')
                                missing_axis_units = True

                        pos = (vij > 0)
                        X_mat[i,j] += vij[pos].sum()
                        X_mat[j,i] -= vij[~pos].sum()

                    # w-fluxes
                    if need_w:
                        if inverted_w:
                            msk_wij[1:,:,ix0:ix1] = Gmask[i,1:,:,ix0:ix1] & Gmask[j,:-1,:,ix0:ix1] & ~w_field[1:,:,ix0:ix1].mask
                        else:
                            msk_wij[:-1,:,ix0:ix1] = Gmask[i,:-1,:,ix0:ix1] & Gmask[j,1:,:,ix0:ix1] & ~w_field[:-1,:,ix0:ix1].mask

                    else:
                        msk_wij[:-1,:,ix0:ix1] = Gmask[i,:-1,:,ix0:ix1] & Gmask[j,1:,:,ix0:ix1] & \
                                                 ~tempvar[:-1,:,ix0:ix1].mask & ~tempvar[1:,:,ix0:ix1].mask

                    if msk_wij.any():

                        if need_w:
                            wij = w_field[msk_wij].data
                        elif w_method == 'divergence':
                            wij = W_PARAM_ADVECTIVE_MIXING

                        if wisveloc:
                            if xunits and yunits:
                                wij *= w_xy[msk_wij]
                                # note: "w_xy" was designed according to the "inverted_w" case (that is "False" is need_w is False)
                            else:
                                print('Warning: cannot compute water fluxes from W velocity field: missing units for x or y axis weigthing')
                                missing_axis_units = True

                        if w_method == 'data':
                            pos = (wij >= 0)
                            X_mat[i,j] += wij[pos].sum()
                            X_mat[j,i] -= wij[~pos].sum()
                        elif w_method == 'divergence':
                            X_mat[i,j] += wij.sum()
                            X_mat[j,i] += wij.sum()


                    # <++++++++++++++++++++++++++++++++++++++++++++++++++++++> #
                    if check_plot:
                        if msk_uij.any() or msk_vij.any() or msk_wij.any():
                            fig, ax = plt.subplots(3, 3, figsize=(15,8))
                            for k,(uvw,msk,txt) in enumerate(zip([uij, vij, wij], [msk_uij, msk_vij, msk_wij], 'UVW')):
                                x = np.zeros(msk.shape)
                                x[msk] = np.abs(uvw)
                                #pid = ax[0,k].pcolormesh(np.count_nonzero(msk, axis=0))
                                pid = ax[0,k].pcolormesh(np.ma.masked_values(x.max(axis=0), 0))
                                ax[0,k].set_xlabel('x')
                                ax[0,k].set_ylabel('y')
                                plt.colorbar(pid)
                                #pid = ax[1,k].pcolormesh(np.count_nonzero(msk, axis=2))
                                pid = ax[1,k].pcolormesh(np.ma.masked_values(x.max(axis=2), 0))
                                ax[1,k].set_xlabel('y')
                                ax[1,k].set_ylabel('z')
                                plt.colorbar(pid)
                                #pid = ax[2,k].pcolormesh(np.count_nonzero(msk, axis=1))
                                pid = ax[2,k].pcolormesh(np.ma.masked_values(x.max(axis=1), 0))
                                ax[2,k].set_xlabel('x')
                                ax[2,k].set_ylabel('z')
                                plt.colorbar(pid)
                                if not reverse_z:
                                    ax[1,k].invert_yaxis()
                                    ax[2,k].invert_yaxis()
                                ax[0,k].set_title(txt)
                            fig.suptitle('Basin {0:} => Basin {1:}'.format(i,j))
                            plt.show()
                    # <++++++++++++++++++++++++++++++++++++++++++++++++++++++> #

                    if missing_axis_units:
                        break

                if missing_axis_units:
                    break

            if missing_axis_units:
                break



            if w_method == 'divergence':
                # Add to each box the exchange flux derived from box divergence, top-bottom scheme
                # + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
                idx = np.concatenate((i_surfbox, [nbasin]))
                for i in range(i_surfbox.size):
                    for k in range(idx[i], idx[i+1]-1):
                        div = X_mat[:-1,k].sum() - X_mat[k,:-1].sum() 
                        # add flux divergence of box k as a downward flux (from box k to box k+1)
                        if div > 0:
                            X_mat[k,k+1] += div
                        else:
                            X_mat[k+1,k] -= div

            # Reduction of horizontal bidirectional fluxes
            # --------------------------------------------

            # get birectional (mixing) fluxes --> min(X[i,j], X[j,i])
            Xsub = np.minimum(X_mat[:-1,:-1], X_mat[:-1,:-1].transpose())
            # compute values to subtract, to get a reduction of 1/h_mix_reduct_factor
            Xsub *= (h_mix_reduct_factor - 1.)/h_mix_reduct_factor
            # do not modify vertical fluxes (=> set Xsub to 0)
            i = list(range(nbasin-1))
            for k in i_surfbox[1:]: i.remove(k-1)
            i = np.array(i)
            Xsub[i,i+1] = 0
            Xsub[i+1,i] = 0
            # REDUCE HORIZONTAL MIXING FLUXES:
            X_mat[:-1,:-1] -= Xsub


            # Check that sum of incoming minus outgoping fluxes is zero
            # ---------------------------------------------------------

            dx = X_mat[:-1,:-1].sum(1) - X_mat[:-1,:-1].sum(0)
            if (dx != 0).any():
                print('mass conservation issue in computed water fluxes. Max error = {:.1e} Sv'.format( \
                                                            FLUX_CONVERSION_FACTOR*(np.abs(dx)).max()) + flx_message)
                flxcorr(X_mat[:-1,:-1], dx)


            # Save in file (append file for each climate state)
            # -------------------------------------------------
            np.savetxt(fout, FLUX_CONVERSION_FACTOR*X_mat.transpose(), fmt='%.7e', delimiter='\t') # => for nicer file visu: fmt='%.2f'
            # note: tranpose "X_mat" because of Fortran opposite {i,j} reading order
            fout.write('\n')


        # <-- end of loop on each input files (ie, climate states) -->


        fout.close()

        if missing_axis_units:
            os.remove(fout.name)

        print('')


    if fgrid is not None:
        fgrid.close()



    # > ============================================================================== < #
    # > Following input files depend only on oceanic basins definition, not on climate < #
    # > ============================================================================== < #


    # Update basin mask to remove missing points of last loaded temperature
    for k in range(nbasin):
        Gmask[k,:,:,:][tempvar.mask] = False


    # Compute and export box areas
    # ----------------------------

    area_top = np.zeros((nbasin+1,), dtype=w_horiz.dtype)
    area_sed = np.zeros((nbasin+1,), dtype=w_horiz.dtype)

    # "area_top" is the horizontal area at the top of boxes
    for k in range(nbasin):
        area_top[k] = w_horiz[:,ix0:ix1][Gmask[k,:,:,ix0:ix1].any(0)].sum()

    # "area_sed" is the area where boxes intercep seafloor
    jj, ii = np.indices(np.array(horiz_shp) - np.array([0,x_overlap]))
    for k in range(nbasin):
        area_sed[k] = w_horiz[:,ix0:ix1][Gmask[k,:,:,ix0:ix1][bottom_idx[:,ix0:ix1], jj, ii]].sum()

    if not xyunits:
        # normalize areas and scale them to modern total ocean volume
        print('areas scaled to modern ocean area')
        area_top *= MODERN_AREA/area_top[i_surfbox].sum()
        area_sed *= MODERN_AREA/area_sed.sum()

    # atmosphere box
    area_top[-1] = ATM_AREA
    area_sed[-1] = ATM_SEDI_AREA

    # Conversion to GEOCLIM units
    area_top *= AREA_CONVERSION_FACTOR
    area_sed *= AREA_CONVERSION_FACTOR

    # Write in file
    np.savetxt(outdir+surf_outfile,    area_top, fmt='%.12e', delimiter='\n')
    np.savetxt(outdir+sedsurf_outfile, area_sed, fmt='%.12e', delimiter='\n')


    # Compute and export box volumes
    # ------------------------------

    box_vol = np.zeros((nbasin+1,), dtype=w_vol.dtype)

    totvol = 0.
    for k in range(nbasin):
        vol = w_vol[:,:,ix0:ix1][Gmask[k,:,:,ix0:ix1]].sum()
        totvol += vol
        box_vol[k] = vol

    if not volunits:
        # normalize volume and scale it to modern total ocean volume
        print('volumes scaled to modern ocean volume')
        box_vol *= MODERN_VOLUME/totvol

    # atmosphere box
    box_vol[-1] = ATM_VOL

    # Conversion to GEOCLIM units
    box_vol *= VOLUME_CONVERSION_FACTOR

    # Write in file
    np.savetxt(outdir+vol_outfile, box_vol, fmt='%.12e', delimiter='\n')


    # Compute and export box depths and average pressure
    # --------------------------------------------------

    box_press = np.zeros(box_depth.shape, box_depth.dtype)

    for k in range(nbasin):
        #
        if k in i_surfbox:
            # impose null pressure for surface boxes 
            box_press[k] = 0
        else:
            # compute average hydrostatic pressure
            volbis = np.zeros(w_vol.shape, w_vol.dtype)
            volbis[:,:,ix0:ix1][Gmask[k,:,:,ix0:ix1]] = w_vol[:,:,ix0:ix1][Gmask[k,:,:,ix0:ix1]]
            zweight = volbis.sum((1,2))
            box_press[k] =  RHO_x_G * (nav_z[:,0,0]*zweight).sum() / zweight.sum()
        #
        if box_depth[k] == -1:  # => bottom boxes
            # compute depth (average bathymetry)
            hmsk = Gmask[k,:,:,ix0:ix1].any(0)
            box_depth[k] = (nav_bathy[0,:,ix0:ix1][hmsk]*w_horiz[:,ix0:ix1][hmsk]).sum() / w_horiz[:,ix0:ix1][hmsk].sum()

    # conversion to GEOCLIM units
    box_depth *= DEPTH_CONVERSION_FACTOR
    box_press *= PRESSURE_CONVERSION_FACTOR

    # Write in file
    np.savetxt(outdir+depth_press_file, np.transpose([box_depth, box_press]), fmt='%.2f', delimiter='\t')


    # Compute and export box-2-box seafloor boundary lengths
    # ------------------------------------------------------

    boundary_length = np.zeros((nbasin+1, nbasin+1), dtype='float32')

    Bmask = Gmask[:,:,:,ix0:ix1][:, bottom_idx[:,ix0:ix1], jj, ii]
    # note: ii, jj are np.indices(...)
    Bmask = np.logical_and(Bmask, ~tempvar[:,:,ix0:ix1][bottom_idx[:,ix0:ix1], jj, ii].mask.reshape((1,)+ii.shape))
    Bmap = np.ma.masked_values(Bmask.astype('int8')[::-1,:,:].cumsum(0).sum(0) - 1,  -1)

    for i in range(nbasin):
        for j in list(range(i))+list(range(i+1,nbasin)):

            # "x" connections
            if periodic_x:
                msk_xij = np.logical_and(Bmask[i,:,-1], Bmask[j,:,0])
                if msk_xij.any():
                    dyij = w_y[isurf,:,ix1-1][msk_xij]
                    boundary_length[i,j] += dyij.sum()
            #
            msk_xij = np.logical_and(Bmask[i,:,:-1], Bmask[j,:,1:])
            if msk_xij.any():
                dyij = w_y[isurf,:,ix0:ix1-1][msk_xij]
                boundary_length[i,j] += dyij.sum()

            # "y" connections
            msk_yij = np.logical_and(Bmask[i,:-1,:], Bmask[j,1:,:])
            if msk_yij.any():
                dxij = w_x[isurf,:-1,ix0:ix1][msk_yij]
                boundary_length[i,j] += dxij.sum()

            # <++++++++++++++++++++++++++++++++++++++++++++++++++++++> #
            if check_plot:
                if msk_xij.any() or msk_yij.any():
                    fig, ax = plt.subplots(3, figsize=(15,18), sharex=True, sharey=True)
                    for k,(wxy,msk,txt) in enumerate(zip([w_x[isurf,:,ix0:ix1-1], w_y[isurf,:-1,ix0:ix1]],
                                                         [msk_xij, msk_yij],
                                                         ['X connection', 'Y connection'])):
                        x = np.zeros(msk.shape)
                        x[msk] = wxy[msk]
                        pid = ax[k].pcolormesh(x)
                        plt.colorbar(pid, ax=ax[k])
                        ax[k].set_title(txt)
                    pid = ax[2].pcolormesh(Bmap, cmap='jet')
                    plt.colorbar(pid, ax=ax[2])
                    ax[2].set_title('Map of basin # on seafloor')
                    for axk in ax:
                        axk.set_ylabel('y')
                    ax[2].set_xlabel('x')
                    fig.suptitle('Basin {0:} => Basin {1:}'.format(i,j))
                    plt.show()
            # <++++++++++++++++++++++++++++++++++++++++++++++++++++++> #

    # "Symmetrize" the matrix (box boundary lengths are unoriented)
    boundary_length += boundary_length.transpose()

    #Write in file
    np.savetxt(outdir+bndlen_outfile, boundary_length, fmt='%.7e', delimiter='\t')


    # Write box definition files
    # --------------------------

    np.savetxt(outdir+deepbox_outfile,    deep_mask, fmt='%i', delimiter='\n')
    np.savetxt(outdir+surfbox_outfile,    surf_mask, fmt='%i', delimiter='\n')
    np.savetxt(outdir+thermobox_outfile,  intm_mask, fmt='%i', delimiter='\n')
    np.savetxt(outdir+polarbox_outfile,   pole_mask, fmt='%i', delimiter='\n')
    np.savetxt(outdir+epicontbox_outfile, epic_mask, fmt='%i', delimiter='\n')
    np.savetxt(outdir+appcontbox_outfile, appc_mask, fmt='%i', delimiter='\n')


    # Save mask fields in netCDF file, if asked to
    # --------------------------------------------

    if export_mask_fields:

        # Create file
        fout = nc.Dataset(outdir+mask_outfile, mode='w', data_model='NETCDF3_CLASSIC')

        # Global attribute
        fout.setncattr('title', 'masks of GEOCLIM basins on GCM oceanic grid')
        fout.setncattr('Nbasins', nbasin)
        fout.setncattr('inputs', '; '.join(([] if mask_nc_file is None else [mask_nc_file]) + \
                                        ([] if grid_file is None else [grid_file]) + \
                                        input_files))

        # Dimensions
        for dim,siz in zip(['x', 'y', 'z'], temp_shape[::-1]):
            fout.createDimension(dim, size=siz)
            fout.createVariable(dim, datatype='int32', dimensions=(dim,))
            fout.variables[dim].setncattr('axis', dim.upper())
            fout.variables[dim].setncattr('units', '-')
            fout.variables[dim][:] = np.arange(1, siz+1, dtype='int32')

        fout.variables['z'].setncattr('positive', 'up' if reverse_z else 'down')

        # Variables
        fout.createVariable('weight', datatype='float32', dimensions=('z', 'y', 'x'),
                            fill_value=(w_vol.fill_value if hasattr(w_vol, 'fill_value') else None))
        fout.variables['weight'].setncattr('long_name', 'volumetric weight of grid cells')
        fout.variables['weight'].setncattr('standard_name', 'weight')
        fout.variables['weight'].setncattr('units', 'm3')
        fout.variables['weight'][:,:,:] = w_vol
        #
        fout.createVariable('mask', datatype='int32', dimensions=('z', 'y', 'x'), fill_value=0)
        fout.variables['mask'].setncattr('long_name', '3D field mask of GEOCLIM basins')
        fout.variables['mask'].setncattr('units', '-')
        fout.variables['mask'][:,:,:] = Gmask[::-1,:,:,:].cumsum(0).sum(0)

        # Close file
        fout.close()

