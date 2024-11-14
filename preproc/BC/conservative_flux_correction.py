'''
Module that contains a couple of methods for correction of flux matrices to ensure mass conservation
'''

import numpy as np


def __get_shape_and_flxdiv(M, dx):
    '''
    local function to get shape of matrix M, vector dx, and check consistency
    '''
    n = M.shape[0]
    if M.shape != (n,n):
        raise ValueError('Function expect a square matrix as first argument')

    # compute the flux divergence "dx", if not provided
    if dx is None:
        dx = M.sum(1) - M.sum(0)
    else:
        if dx.shape != (n,):
            raise ValueError('Argument "dx" must be rank-1 with same dimension as matrix "M"')

    return n, dx


def PositiveLeastSquare_correction(M: np.ndarray, dx=None):
    '''
    Compute and apply to matrix "M" an additive correction 'delta', so that the
    corrected flux divergence "(M+delta).sum(1) - (M+delta).sum(0)" is 0.
    'delta' is computed by minimizing its L2 norm (sum of the square of all its elements)
    The correction 'delta" is modified afterward to ensure that M stay positive.

    "dx" is the flux divergence "M.sum(1) - M.sum(0)", that can be passed to the
    function, if already computed.
    '''

    # Get shape and flux divergence "dx" (if not provided)
    n, dx = __get_shape_and_flxdiv(M, dx)

    # Compute and apply the least-square additive correction
    delta = (dx.reshape((1,n)) - dx.reshape((n,1))) / (2*n)
    M += delta

    # Modify the correction to keep "M" positive
    i,j = np.argwhere(M<0).transpose()
    M[j,i] -= M[i,j]
    M[i,j] = 0



def RelativeLeastSquare_correction(M: np.ndarray, dx=None):
    '''
    Compute and apply to matrix "M" a (element-by-element) multiplicative correction
    'F', so that the corrected flux divergence "(F*M).sum(1) - (F*M).sum(0)" is 0
    ('F*M' being the element-by-element multiplication, not matrix multiplication).
    'F' is computed by minimizing the sum of "(F[i,j] - 1)^2"

    "dx" is the flux divergence "M.sum(1) - M.sum(0)", that can be passed to the
    function, if already computed.

    WARNING: this method have shown to return erroneous results in the case where
    "M" contains indepedent clusters, with no exchange between them.
    '''

    # Get shape and flux divergence "dx" (if not provided)
    n, dx = __get_shape_and_flxdiv(M, dx)

    # Build matrix to invert
    N = (M**2 + M.transpose()**2)/2
    S = N[:-1,-1:] - N[:-1,:-1]
    S[np.diag_indices(n-1)] += N[:-1,:].sum(1)

    # Intermediate vector
    x = np.zeros((n,), M.dtype)
    # Inversion:
    # +++++++++++++++++++++++++++++++++++++++++++ #
    x[:-1] = np.matmul(np.linalg.inv(S), dx[:-1])
    x[-1] = -x[:-1].sum()
    # +++++++++++++++++++++++++++++++++++++++++++ #

    # Correction factor:
    F = 1 + (x.reshape((1,n)) - x.reshape((n,1)))*M/2
    if (F < 0).any():
        raise ValueError('"relative least-square" correction yields negative exchange fluxes\nThe needed corrections appear to be to big for this method. Try a more stable correction method (Positive least-square for instance)')

    ######
    M *= F
    ######

    if (np.abs((M.sum(1)-M.sum(0)) / np.maximum(M.sum(0), M.sum(1))) > 1e-7).any():
        raise ValueError('"Relative least-square" correction method have returned erroneous results')
