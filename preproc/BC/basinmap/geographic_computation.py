import numpy as np


#RE = 6371e3 # (m) Earth Radius
#PI = np.pi


def cos(x): return np.cos(np.radians(x))
def sin(x): return np.sin(np.radians(x))


def dist(lon0, lat0, lon1, lat1):
    """
    This function computes the normalized distance between 2 points "0" and "1"
    on Earth with spherical approximation.
    The distance is normalized for a sphere radius of 1.
    Arguments: lon0, lat0, lon1, lat1 (all in degrees)
    """
    return np.arccos(sin(lat0)*sin(lat1) + cos(lat0)*cos(lat1)*cos((lon1-lon0))) 


class NoNeighbourgError(Exception):
    """
    Raise this exception when no neighbourg are found by "get_closest_B_from_A"
    function.
    """
    pass


def get_closest_B_from_A(lon_A: np.ndarray, lat_A: np.ndarray, lon_B: np.ndarray, lat_B: np.ndarray,
                         mask_A=None, mask_B=None, list_A=None, list_B=None):
    """
    For each points of a 2-D grid "A" where "mask_A" is True, determine the
    closest point of grid "B" where "mask_B" is True, under spherical Earth
    approximation.
    Alternatively to "mask_A"("mask_B"), a list of indices "list_A"("list_B")
    can be given, indicating the indices of the True points in grid "A"("B").
    The coordinates of the grid must be longitude-latitude, in DEGREES. This
    function also works for curvilinear grid.
    Returns two arrays, "i" and "j", of same shape than "mask_A", indicating
    for each "A" point the indices of the closet "B" point (on "B" grid).
    as longitude and latitude are given for each point.

    INPUT ARGUMENTS:

        lon_A: numpy-array, rank 1 or 2
            Longitude of "A" grid. Its shape must be "(njA,)" (rectilinear grid)
            or "(niA,njA)" (curvilinear grid).

        lat_A: numpy-array, rank 1 or 2
            Latitude of "A" grid. Its shape must be "(niA,)" (rectilinear grid)
            or "(niA,njA)" (curvilinear grid).

        lon_B: numpy-array, rank 1 or 2
            Longitude of "B" grid. Its shape must be "(njB,)" (rectilinear grid)
            or "(niB,njB)" (curvilinear grid).

        lat_B: numpy-array, rank 1 or 2
            Latitude of "B" grid. Its shape must be "(niB,)" (rectilinear grid)
            or "(niB,njB)" (curvilinear grid).

        mask_A: numpy-array, rank-2, type bool
            Points of "A" grid where the closest "B" points will be searched.
            Its shape must be "(niA,njA)".

        mask_B: numpy-array, rank-2, type bool
            Selected points of "B" grid that can be closest to "A" points.
            Other "B" points will be ignored.
            Its shape must be "(niB,njB)".

        list_A: numpy-array, or list of indices (integer)
            Alternative to "mask_A": list of "mask_A" indices.
            Must be an array of shape (n,2), or tuple of 2 sequences of same
            length.

        list_B: numpy-array, or list of indices (integer)
            Alternative to "mask_B": list of "mask_B" indices.
            Must be an array of shape (m,2), or tuple of 2 sequences of same
            length.

    OUTPUT ARGUMENTS:

        iclosest, jclosest: numpy-arrays, shape (niA,njA), type integer
            i and j indices of the closest "B" point, for each point of grid "A"
    """


    # Compatibility checks
    # ====================
    
    # "A" points coordinates
    # ----------------------
    if lon_A.ndim==1 and lat_A.ndim==1:
        curvilinear_A = False
        nAi,nAj = lat_A.size, lon_A.size

    elif lon_A.ndim==2 and lat_A.ndim==2:
        curvilinear_A = True
        if lon_A.shape != lat_A.shape:
            raise ValueError('curvilinear coordinates "lon_A" and "lat_A" must have same shape')
        else:
            nAi,nAj = lon_A.shape

    else:
        raise ValueError('"lon_A" and "lat_A" must be either both rank-1 or both rank-2')

    # "B" points coordinates
    # ----------------------
    if lon_B.ndim==1 and lat_B.ndim==1:
        curvilinear_B = False
        nBi,nBj = lat_B.size, lon_B.size

    elif lon_B.ndim==2 and lat_B.ndim==2:
        curvilinear_B = True
        if lon_B.shape != lat_B.shape:
            raise ValueError('curvilinear coordinates "lon_B" and "lat_B" must have same shape')
        else:
            nBi,nBj = lon_B.shape

    else:
        raise ValueError('"lon_B" and "lat_B" must be either both rank-1 or both rank-2')

    # "A" points mask
    # ---------------
    if mask_A is not None:
        if mask_A.ndim != 2:
            raise ValueError('masks arguments must be rank-2')
        elif mask_A.shape != (nAi,nAj):
            raise ValueError('shape of "mask_A" inconsistent with coordinates')

        list_A = np.argwhere(mask_A)

    elif list_A is not None:

        wrong = True

        if isinstance(list_A, np.ndarray):
            if list_A.ndim == 2:
                if list_A.shape[1] == 2:
                    wrong = False

        elif len(list_A) == 2:
            list_A = np.array(list_A).transpose()
            wrong = False

        if wrong:
            raise ValueError('Cannot interpret "list_A" as an numpy-array or list of 2D indices')

    else:
        raise ValueError('Expect either "mask_A" or "list_A" argument.')

    # "B" points mask
    # ---------------
    if mask_B is not None:
        if mask_B.ndim != 2:
            raise ValueError('masks arguments must be rank-2')
        elif mask_B.shape != (nBi,nBj):
            raise ValueError('shape of "mask_B" inconsistent with coordinates')

        list_B = np.argwhere(mask_B)

    elif list_B is not None:

        wrong = True

        if isinstance(list_B, np.ndarray):
            if list_B.ndim == 2:
                if list_B.shape[1] == 2:
                    wrong = False

        elif len(list_B) == 2:
            list_B = np.array(list_B).transpose()
            wrong = False

        if wrong:
            raise ValueError('Cannot interpret "list_B" as an numpy-array or list of 2D indices')

    else:
        raise ValueError('Expect either "mask_B" or "list_B" argument.')


    # Computation
    # ===========

    # Initialize output arguments
    iclosest, jclosest = np.zeros((2,nAi,nAj), dtype='int32')

    # pre-computation
    if list_B.size==0:
        raise NoNeighbourgError('no "B" point in the given mask or list')
    #
    iA, jA = list_A.transpose()
    n = iA.size
    if curvilinear_A:
        coords_A = lon_A[iA,jA].reshape((n,1)), lat_A[iA,jA].reshape((n,1))
    else:
        coords_A = lon_A[jA].reshape((n,1)), lat_A[iA].reshape((n,1))
    #
    iB, jB = list_B.transpose()
    n = iB.size
    if curvilinear_B:
        coords_B = lon_B[iB,jB].reshape((1,n)), lat_B[iB,jB].reshape((1,n))
    else:
        coords_B = lon_B[jB].reshape((1,n)), lat_B[iB].reshape((1,n))

    # Get closest 'B' point (axis #1) from each 'A' point (axis #0)
    dmat = dist(*coords_A, *coords_B)
    k = np.argmin(dmat, axis=1)
    iclosest[iA,jA], jclosest[iA,jA] = iB[k], jB[k]


    return iclosest, jclosest

