import numpy as np
import sys

list_file = sys.argv[1:]

#IBOX_ATL = (3, 4, 5)
#IBOX_PAC = (9, 10, 11)
IBOX_ATL = (3, 4,  5,  3,  4,  5)
IBOX_PAC = (9, 10, 11, 12, 13, 14)

for fname in list_file:

    if fname[-19:] == '_human-friendly.dat':
        ispl = -19
        fmtargs = {'fmt': '%6.2f'}
    elif fname[-4:] == '.dat':
        ispl = -4
        fmtargs = {'fmt': '%.7e', 'delimiter': '\t'}
    else:
        raise ValueError('Unrecognized name "{:}"'.format(fname))

    X = np.loadtxt(fname)
    nk,k = X.shape # (configs*boxes, boxes)

    #if k != 29:
    #    raise ValueError('{0:} boxes found in file {1:}.\nThis script was designed for, and can only be applied to GEOCLIM 90Ma 29 boxes configuration'.format(k, fname))

    n = nk//k
    if n*k != nk:
        raise ValueError('Cannot understand a series of exchanges matrices of shape ({0:},{1:}).\nExpect (nconfigs*nboxes, nboxes).'.format(nk,k))

    # (configs, rows, colums)
    X = X.reshape((n,k,k))

    list_XAP, list_XPA, list_Fbid = [], [], []
    for ia, ip in zip(IBOX_ATL, IBOX_PAC):
        list_XAP.append(X[:,ia,ip].copy())
        list_XPA.append(X[:,ip,ia].copy())
        list_Fbid.append(np.minimum(list_XAP[-1], list_XPA[-1]))

    for redfact in (0.25,):#(0.5, 0.25, 0):
        for ia, ip, XAP, XPA, Fbid in zip(IBOX_ATL, IBOX_PAC, list_XAP, list_XPA, list_Fbid):
            X[:,ia,ip] = XAP - (1-redfact)*Fbid
            X[:,ip,ia] = XPA - (1-redfact)*Fbid

        np.savetxt(fname[:ispl]+'_APTx{:}'.format(redfact)+fname[ispl:], X.reshape((n*k,k)), **fmtargs)

