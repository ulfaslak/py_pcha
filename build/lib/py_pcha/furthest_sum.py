u"""Python implemtation of Morten Morups matlab code for PCHA.

Translated to Python by Ulf Aslak Jensen
Copyright (C) Morten Morup and Technical University of Denmark, 2010
"""

from __future__ import division
import numpy as np
from numpy.matlib import repmat


def furthest_sum(K, noc, i, exclude=[]):
    """Furthest sum algorithm, to efficiently generat initial seed/archetypes.

    Parameters
    ----------
    K : numpy 2d-array
        Either a data matrix or a kernel matrix.

    noc : int
        Number of candidate archetypes to extract.

    i : int
        inital observation used for to generate the FurthestSum.

    exclude : numpy.1darray
        Entries in K that can not be used as candidates.

    Output
    ------
    i : int
        The extracted candidate archetypes
    """
    def max_ind_val(l):
        return max(zip(range(len(l)), l), key=lambda x: x[1])

    I, J = K.shape
    index = np.array(range(J))
    index[exclude] = 0
    index[i] = -1
    ind_t = i
    sum_dist = np.zeros((1, J), np.complex128)

    if J > noc * I:
        Kt = K
        Kt2 = np.sum(Kt**2, axis=0)
        for k in range(1, noc + 11):
            if k > noc - 1:
                Kq = np.dot(Kt[:, i[0]], Kt)
                sum_dist -= np.lib.scimath.sqrt(Kt2 - 2 * Kq + Kt2[i[0]])
                index[i[0]] = i[0]
                i = i[1:]
            t = np.where(index != -1)[0]
            Kq = np.dot(Kt[:, ind_t].T, Kt)
            sum_dist += np.lib.scimath.sqrt(Kt2 - 2 * Kq + Kt2[ind_t])
            ind, val = max_ind_val(sum_dist[:, t][0].real)
            ind_t = t[ind]
            i.append(ind_t)
            index[ind_t] = -1
    else:
        if I != J or np.sum(K - K.T) != 0:  # Generate kernel if K not one
            Kt = K
            K = np.dot(Kt.T, Kt)
            K = np.lib.scimath.sqrt(
                repmat(np.diag(K), J, 1) - 2 * K + \
                repmat(np.mat(np.diag(K)).T, 1, J)
            )

        Kt2 = np.diag(K)  # Horizontal
        for k in range(1, noc + 11):
            if k > noc - 1:
                sum_dist -= np.lib.scimath.sqrt(Kt2 - 2 * K[i[0], :] + Kt2[i[0]])
                index[i[0]] = i[0]
                i = i[1:]
            t = np.where(index != -1)[0]
            sum_dist += np.lib.scimath.sqrt(Kt2 - 2 * K[ind_t, :] + Kt2[ind_t])
            ind, val = max_ind_val(sum_dist[:, t][0].real)
            ind_t = t[ind]
            i.append(ind_t)
            index[ind_t] = -1

    return i

if __name__ == "__main__":
    X=np.array(
        [[0, 0, 1, 0.5, 4, 5, 2, 5, 3, 6, 2, 6, 3, 6, 8, 8, 6, 5, 3, 5, 3, 5, 2, 4, 5],
         [0, 1, 0, 0.6, 5, 2, 6, 7, 2, 6, 1, 3, 4, 5, 6, 6, 4, 7, 8, 8, 6, 4, 4, 4, 2]]
    )
    noc = 3
    N, M = X.shape
    I = range(M)
    U = range(M)
    print PCHA(X, noc, [0])