import numpy as np
from units import length_units, UnknownUnitError



##########################################
# FUNCTIONS SPECIFIC TO PARTICULAR GRIDS #
##########################################

def remove_orca_north_fold(x: np.ndarray, pivot: str):
    '''
    Put 0 in points of numpy-array "x" that are outside the regular
    'inner domain' mask at North fold boundary of ORCA grid, which depends on
    whether the North fold is on a T point (pivot='T') or an F point (pivot='F')
    NOTE: x-axis ("longitude") is expected to be the last dimension of "x", and
          must have an even number of points, y-axis ("latitude") is expected to
          be the second-to-last dimension of "x".
    '''

    if pivot not in ('T', 'T-point', 'F', 'F-point'):
        raise ValueError('"pivot" must be "T" ["T-point"] or "F" ["F-point"]')

    x[..., -1, :] = 0

    if pivot[:1] == 'T':
        ihalf = x.shape[-1]//2
        x[..., -2, ihalf+1:] = 0



#######################
# AUXILIARY FUNCTIONS #
#######################

def find_bathy(var_mask=None, z=None, cell_thick=None, invert_zaxis=False):
    '''
    Return a 2D (horizontal) array giving the seafloor depth (bathymetry) of each point,
    and the 2D array of vertical indices at which the seafloor is reached.
    The bathymetry is determined as follows:
      1. Find the index the last 'False' of 3D boolean array -- either "var_mask" (if given),
         or z.mask or cell_thick.mask -- in the vertical (1st) dimension
         => yields a 2D "horizontal" index array `idx`.
      2. Pick the value of depth for each "bottom index". Depth being eihter directly "z"
         (if given) or computed from cell thickness "cell_thick".
         ie, `seafloor_depth = z[idx[i,j], i, j]` for all i,j
         If all the value of 3D mask in one column are 'True', the depth of that point
         will be '0'
    "var_mask" must be a rank-3 boolean array, and is optional if a 3D "z" or "cell_thick"
    array is given.
    "z" and "cell_thick" must be either rank-1 or rank-3 array, in which case the last two
    dimensions can be denerated (size-1).
    If neither "z" nor "cell_thick" is given, "var_mask" is expected. Then, only the "bottom
    indices" array is computed and returned.

    The vertical dimension must be the 1st one (C-indexing), it must be positive downward,
    and it must be in increasing order (i.e., from surface to bottom), unless specified
    "invert_zaxis=True".
    '''

    # order of z axis
    # ---------------

    slc1D = slice(None, None, -1) if invert_zaxis else slice(None)
    slc3D = (slc1D, slice(None), slice(None))
     

    # Check call signature and arguments type and shape:
    # -------------------------------------------------

    # Argument type:
    for arg in (var_mask, z, cell_thick):
        if not (arg is None or isinstance(arg, np.ndarray)):
            raise TypeError('all input arguments of function "find_bathy" mut be numpy-arrays')

    # Call signature and argument shape

    compute_bathy = True

    if z is None and cell_thick is None:

        if var_mask is None:
            raise ValueError('at least one argument must be passed to function "find_bathy"')
        elif var_mask.ndim != 3:
            raise ValueError('"var_mask" argument must be rank-3')

        nz = var_mask.shape[0]
        horiz_shape = var_mask.shape[1:]
        compute_bathy = False

    elif cell_thick is None: # -> only "z" is given

        if z.ndim==1 or (z.ndim==3 and z.shape[1]==1 and z.shape[2]==1):
            if var_mask is None:
                raise ValueError('argument "var_mask" needed to compute bathymetry from 1D depth array')
            else:
                if var_mask.ndim != 3:
                    raise ValueError('"var_mask" argument must be rank-3')
                if var_mask.shape[0] != z.shape[0]:
                    raise ValueError('1st dimension of array "var_mask" must match the size of "z" array')

            nz = z.size
            horiz_shape = var_mask.shape[1:]
            z_is_1D = True

        elif z.ndim==3:
            if var_mask is not None:
                if var_mask.shape != z.shape:
                    raise ValueError('3D input argument "var_mask" and "z" must have the same shape')

            nz = z.shape[0]
            horiz_shape = z.shape[1:]
            z_is_1D = False

        else:
            raise ValueError('"z" argument must be rank 1 or 3')

    else:
        if z is not None:
            print('warning: in function "find_bathy", argument "z" is ignored, bathymetry computed with argument "cell_thick"')

        if cell_thick.ndim==1 or (cell_thick.ndim==3 and cell_thick.shape[1]==1 and cell_thick.shape[2]==1):
            if var_mask is None:
                raise ValueError('argument "var_mask" needed to compute bathymetry from 1D cell thickness array')
            else:
                if var_mask.ndim != 3:
                    raise ValueError('"var_mask" argument must be rank-3')
                if var_mask.shape[0] != cell_thick.shape[0]:
                    raise ValueError('1st dimension of array "var_mask" must match the size of "cell_thick" array')

            nz = cell_thick.size
            horiz_shape = var_mask.shape[1:]
            z_is_1D = True

        elif cell_thick.ndim==3:
            if var_mask is not None:
                if var_mask.shape != cell_thick.shape:
                    raise ValueError('3D input argument "var_mask" and "cell_thick" must have the same shape')

            nz = cell_thick.shape[0]
            horiz_shape = cell_thick.shape[1:]
            z_is_1D = False

        else:
            raise ValueError('"cell_thick" argument must be rank 1 or 3')

        # Compute z (depth at bottom of cells) from cell thickness
        #<><><><><><><><><><><><><><><><><><><><><>#
        if cell_thick.ndim == 1:
            z = cell_thick[slc1D].cumsum(0)[slc1D]
        else:
            z = cell_thick[slc3D].cumsum(0)[slc3D]
       #<><><><><><><><><><><><><><><><><><><><><>#

        cell_thick = None


    # Computation
    # -----------

    # Initialization
    if var_mask is not None:
        cellmask = var_mask
    else:
        if z is not None:
            if hasattr(z, 'mask'):
                cellmask = z.mask
            else:
                cellmask = False

    if np.size(cellmask)==1:
        cellmask = np.zeros((nz,1,1), dtype=bool)

    if invert_zaxis:
        bottom_idx = -1*np.ones(horiz_shape, dtype='int16')
    else:
        bottom_idx = nz*np.ones(horiz_shape, dtype='int16')

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
    # Loop on vertical dimension of "var_mask" => DETERMINE POSITION OF SEAFLOOR (bottom indices)
    for i in range(nz)[slc1D]:
        bottom_idx[~cellmask[i,:,:]] = i
    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #


    if compute_bathy:

        # Trick: add an extra element "0" and the end of z, so that if
        # bottom_idx=="last index" (all elements are True) => z[bottom_idx] = 0
        if invert_zaxis:
            z = np.concatenate((np.zeros((1,)+z.shape[1:], z.dtype), z), axis=0)
            bottom_idx += 1
        else:
            z = np.concatenate((z, np.zeros((1,)+z.shape[1:], z.dtype)), axis=0)

        #<><><><><><><><><><><><><><><>#
        if z_is_1D:
            if z.ndim == 1:
                key = (bottom_idx,)
            else: # z.ndim==3
                key = (bottom_idx, 0, 0)

        else:
            ii, jj = np.indices(horiz_shape)
            key = (bottom_idx, ii, jj)

        bathy = z[key]
        #<><><><><><><><><><><><><><><>#

        if invert_zaxis:
            bottom_idx -= 1
            # Replace "last index" values by "0"
            bottom_idx[bottom_idx==-1] = 0
        else:
            # Replace "last index" values by "-1"
            bottom_idx[bottom_idx==nz] = -1

        return bathy, bottom_idx


    else:

        if invert_zaxis:
            # Replace "last index" values by "0"
            bottom_idx[bottom_idx==-1] = 0
        else:
            # Replace "last index" values by "-1"
            bottom_idx[bottom_idx==nz] = -1

        return bottom_idx


# ========== #


def get_var(*args, varname=None, raise_error=False):
    '''
    Try to get a variable "varname" (string) in one of the provided netCDF4 dataset
    (*args is all the datasets, provided as separated arguments).
    A None is got in the following cases:
      - dataset is None
      - varname is None
      - varname not found in dataset
    The first dataset that doesn't yield a None will be kept. If all do, return None.
    '''
    for dataset in args:
        if varname is None or dataset is None:
            var = None
        else:
            try:
                var = dataset[varname]
            except IndexError:
                var = None

        if var is not None:
            return var

    if raise_error:
        raise IndexError('Variable "'+str(varname)+'" cannot be found in any of the input datasets')


# ========== #


def get_axis_weighting(axis: str, wghvar, ref_dims, ref_shp, axis_bnds_var=None, axis_var=None, t_dim=None):
    '''
    Get the weighting of one axis (x, y or z), check that the dimensions and size
    match the 3D reference shape.
    input argument:
      - "axis": name of the axis (string)
      - "wghvar": netCDF4 variable (in netCDF4 Dataset) of the weighting.
        > if wghvar is None: create uniform weighting 
      - "ref_dims": tuple of dimensions of the reference 3D variable
      - "ref_shp": shape (tuple) of the reference 3D variable arrays
    optional input arguments:
      - "axis_bnds_var": netCDF4 variable (in netCDF4 Dataset) giving the axis 
          bounds (must be 2D, with 2nd dimension, size-2, being "bounds").
          Get axis weighting as the difference (np.diff) of axis bounds in case
          "wghvar" is None.
      - "axis_var": netCDF4 variable (in netCDF4 Dataset) giving the axis
          coordinates. Used to get "axis_bnds_var" units in case "axis_bnds_var"
          doesn't have a "units" attribute. 
      - "t_dim": potential time dimension in case the "wghvar" is 4D, and will
          be averaged on the time dimension.
    returns:
      - 3D array (with shape "ref_shp") of the axis weighting
      - boolean, indicated whether axis weighting has units (expect "meter")
    '''

    if wghvar is not None:

        shp = list(ref_shp)

        if wghvar.ndim == 3:
            wgh = wghvar[:,:,:]

        elif wghvar.ndim == 2:
            if wghvar.dimensions[0]==ref_dims[1] and wghvar.dimensions[1]==ref_dims[2]: # weighting dimensions are {x,y}
                shp[0] = 1 # remove z dimension
            elif wghvar.dimensions[0]==ref_dims[0] and wghvar.dimensions[1]==ref_dims[2]: # weighting dimensions are {x,z}
                shp[1] = 1 # remove y dimension
            elif wghvar.dimensions[0]==ref_dims[0] and wghvar.dimensions[1]==ref_dims[1]: # weighting dimensions are {y,z}
                shp[2] = 1 # remove x dimension
            else:
                raise ValueError('Cannot identify dimensions of "'+axis+'_weight" variable (do not match "temperature" variable)')

            wgh = wghvar[:,:].reshape(shp)

        elif wghvar.ndim == 1:
            if wghvar.dimensions[0]==ref_dims[2]: # weighting dimension is {x}
                shp[0] = 1 # remove z dimension
                shp[1] = 1 # remove y dimension
            elif wghvar.dimensions[0]==ref_dims[1]: # weighting dimension is {y}
                shp[0] = 1 # remove z dimension
                shp[2] = 1 # remove x dimension
            elif wghvar.dimensions[0]==ref_dims[0]: # weighting dimension is {z}
                shp[1] = 1 # remove y dimension
                shp[2] = 1 # remove x dimension
            else:
                raise ValueError('Cannot identify dimensions of "'+axis+'_weight" variable (do not match "temperature" variable)')

            wgh = wghvar[:].reshape(shp)

        else:
            if wghvar.ndim == 4 and t_dim is not None:
                wgh = wghvar[:,:,:,:].mean(t_dim)
            else:
                raise ValueError('Illegal number of dimensions for variable "'+axis+'_weight"')

        try:
            wgh = length_units.convert(wgh, wghvar.units)
            has_units = True
        except (AttributeError, UnknownUnitError):
            has_units = False

    elif axis_bnds_var is not None:
        if axis_bnds_var.ndim != 2:
            raise ValueError('Cannot handle "'+axis+'_bounds" variable that is not 2D')
        elif axis_bnds_var.shape[1] != 2:
            raise ValueError('2nd dimension of '+axis+' bounds variable (argument "'+axis+'_bounds") must be size-2')
        else:
            wgh = np.abs(np.diff(axis_bnds_var[:], axis=1))
            try:
                wgh = length_units.convert(wgh, axis_bnds_var.units)
                has_units = True
            except AttributeError:
                if axis_var is None:
                    has_units = False
                else:
                    try:
                        wgh = length_units.convert(wgh, axis_var.units)
                        has_units = True
                        print('WARNING: units of '+axis+' bounds variable not found. ASSUME SAME UNIT AS '+axis+' VARIABLE.')
                    except UnknownUnitError:
                        has_units = False

            except UnknownUnitError:
                has_units = False

            shp = list(ref_shp)
            if axis_bnds_var.dimensions[0]==ref_dims[2]: # weighting dimension is {x}
                shp[0] = 1 # remove z dimension
                shp[1] = 1 # remove y dimension
            elif axis_bnds_var.dimensions[0]==ref_dims[1]: # weighting dimension is {y}
                shp[0] = 1 # remove z dimension
                shp[2] = 1 # remove x dimension
            elif axis_bnds_var.dimensions[0]==ref_dims[0]: # weighting dimension is {z}
                shp[1] = 1 # remove y dimension
                shp[2] = 1 # remove x dimension
            else:
                raise ValueError('Cannot identify dimensions of "'+axis+'_bounds" variable (do not match "temperature" variable)')

            wgh = wgh.reshape(shp)

    else:
        wgh = 1
        print('WARNING: uniform weighting assumed in '+axis+' dimension')
        has_units = False

    return wgh*np.ones(ref_shp), has_units


# ========== #


def average(*args):
    '''
    Return the unweighted average of all input arguments (that must be
    numpy.ma.core.MaskedArray of shame shape), removing masked value from
    the average. Where all arrays are masked, the average is 0.
    '''
    ave = np.zeros(args[0].shape, dtype=float)
    cnt = np.zeros(args[0].shape, dtype='int8')
    msk = np.zeros(args[0].shape, dtype=bool)
    for x in args:
        #x = np.ma.array(x)
        ave[~x.mask] += x.data[~x.mask]
        cnt[~x.mask] += 1
        msk[~x.mask] = True

    ave[msk] /= cnt[msk]
    return ave


# ========== #


def minimum(*args):
    '''
    Return the point-by-point minimum of all input arguments (that must be
    numpy.ma.core.MaskedArray of shame shape), removing masked value from
    the list to compute the miminum.
    Where all arrays are masked, the minimum is 0.
    '''
    mini = np.zeros(args[0].shape, dtype=float)
    msk  = np.zeros(args[0].shape, dtype=bool)
    for x in args:
        #x = np.ma.array(x)
        msk1 = np.logical_and(~x.mask, ~msk)
        msk2 = np.logical_and(~x.mask, msk)
        mini[msk1] = x.data[msk1]
        mini[msk2] = np.minimum(mini[msk2], x.data[msk2])
        msk[msk1] = True

    return mini


