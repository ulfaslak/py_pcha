"""Principal Convex Hull Analysis (PCHA) / Archetypal Analysis."""

from __future__ import division
import numpy as np
from scipy.sparse import csr_matrix
from datetime import datetime as dt
import time

from .furthest_sum import furthest_sum


def PCHA(X, noc, I=None, U=None, delta=0, verbose=False, conv_crit=1E-6, maxiter=500):
    """Return archetypes of dataset.

    Note: Commonly data is formatted to have shape (examples, dimensions).
    This function takes input and returns output of the transposed shape,
    (dimensions, examples).

    Parameters
    ----------
    X : numpy.2darray
        Data matrix in which to find archetypes

    noc : int
        Number of archetypes to find

    I : 1d-array
        Entries of X to use for dictionary in C (optional)

    U : 1d-array
        Entries of X to model in S (optional)


    Output
    ------
    XC : numpy.2darray
        I x noc feature matrix (i.e. XC=X[:,I]*C forming the archetypes)

    S : numpy.2darray
        noc x length(U) matrix, S>=0 |S_j|_1=1

    C : numpy.2darray
        noc x length(U) matrix, S>=0 |S_j|_1=1

    SSE : float
        Sum of Squared Errors

    varexlp : float
        Percent variation explained by the model
    """
    def S_update(S, XCtX, CtXtXC, muS, SST, SSE, niter):
        """Update S for one iteration of the algorithm."""
        noc, J = S.shape
        e = np.ones((noc, 1))
        for k in range(niter):
            SSE_old = SSE
            g = (np.dot(CtXtXC, S) - XCtX) / (SST / J)
            g = g - e * np.sum(g.A * S.A, axis=0)

            S_old = S
            while True:
                S = (S_old - g * muS).clip(min=0)
                S = S / np.dot(e, np.sum(S, axis=0))
                SSt = S * S.T
                SSE = SST - 2 * np.sum(XCtX.A * S.A) + np.sum(CtXtXC.A * SSt.A)
                if SSE <= SSE_old * (1 + 1e-9):
                    muS = muS * 1.2
                    break
                else:
                    muS = muS / 2

        return S, SSE, muS, SSt

    def C_update(X, XSt, XC, SSt, C, delta, muC, mualpha, SST, SSE, niter=1):
        """Update C for one iteration of the algorithm."""
        J, nos = C.shape

        if delta != 0:
            alphaC = np.sum(C, axis=0).A[0]
            C = np.dot(C, np.diag(1 / alphaC))

        e = np.ones((J, 1))
        XtXSt = np.dot(X.T, XSt)

        for k in range(niter):

            # Update C
            SSE_old = SSE
            g = (np.dot(X.T, np.dot(XC, SSt)) - XtXSt) / SST

            if delta != 0:
                g = np.dot(g, np.diag(alphaC))
            g = g.A - e * np.sum(g.A * C.A, axis=0)

            C_old = C
            while True:
                C = (C_old - muC * g).clip(min=0)
                nC = np.sum(C, axis=0) + np.finfo(float).eps
                C = np.dot(C, np.diag(1 / nC.A[0]))

                if delta != 0:
                    Ct = C * np.diag(alphaC)
                else:
                    Ct = C

                XC = np.dot(X, Ct)
                CtXtXC = np.dot(XC.T, XC)
                SSE = SST - 2 * np.sum(XC.A * XSt.A) + np.sum(CtXtXC.A * SSt.A)

                if SSE <= SSE_old * (1 + 1e-9):
                    muC = muC * 1.2
                    break
                else:
                    muC = muC / 2

            # Update alphaC
            SSE_old = SSE
            if delta != 0:
                g = (np.diag(CtXtXC * SSt).T / alphaC - np.sum(C.A * XtXSt.A)) / (SST * J)
                alphaC_old = alphaC
                while True:
                    alphaC = alphaC_old - mualpha * g
                    alphaC[alphaC < 1 - delta] = 1 - delta
                    alphaC[alphaC > 1 + delta] = 1 + delta

                    XCt = np.dot(XC, np.diag(alphaC / alphaC_old))
                    CtXtXC = np.dot(XCt.T, XCt)
                    SSE = SST - 2 * np.sum(XCt.A * XSt.A) + np.sum(CtXtXC.A * SSt.A)

                    if SSE <= SSE_old * (1 + 1e-9):
                        mualpha = mualpha * 1.2
                        XC = XCt
                        break
                    else:
                        mualpha = mualpha / 2

        if delta != 0:
            C = C * np.diag(alphaC)

        return C, SSE, muC, mualpha, CtXtXC, XC

    N, M = X.shape
    

    if I is None:
        I = range(M)
    if U is None:
        U = range(M)

    SST = np.sum(X[:, U] * X[:, U])

    # Initialize C
    try:
        i = furthest_sum(X[:, I], noc, [int(np.ceil(len(I) * np.random.rand()))])
    except IndexError:
        class InitializationException(Exception): pass
        raise InitializationException("Initialization does not converge. Too few examples in dataset.")

    j = range(noc)
    C = csr_matrix((np.ones(len(i)), (i, j)), shape=(len(I), noc)).todense()

    XC = np.dot(X[:, I], C)

    muS, muC, mualpha = 1, 1, 1

    # Initialise S
    XCtX = np.dot(XC.T, X[:, U])
    CtXtXC = np.dot(XC.T, XC)
    S = -np.log(np.random.random((noc, len(U))))
    S = S / np.dot(np.ones((noc, 1)), np.mat(np.sum(S, axis=0)))
    SSt = np.dot(S, S.T)
    SSE = SST - 2 * np.sum(XCtX.A * S.A) + np.sum(CtXtXC.A * SSt.A)
    S, SSE, muS, SSt = S_update(S, XCtX, CtXtXC, muS, SST, SSE, 25)

    # Set PCHA parameters
    iter_ = 0
    dSSE = np.inf
    t1 = dt.now()
    varexpl = (SST - SSE) / SST

    if verbose:
        print('\nPrincipal Convex Hull Analysis / Archetypal Analysis')
        print('A ' + str(noc) + ' component model will be fitted')
        print('To stop algorithm press control C\n')

    dheader = '%10s | %10s | %10s | %10s | %10s | %10s | %10s | %10s' % ('Iteration', 'Expl. var.', 'Cost func.', 'Delta SSEf.', 'muC', 'mualpha', 'muS', ' Time(s)   ')
    dline = '-----------+------------+------------+-------------+------------+------------+------------+------------+'

    while np.abs(dSSE) >= conv_crit * np.abs(SSE) and iter_ < maxiter and varexpl < 0.9999:
        if verbose and iter_ % 100 == 0:
            print(dline)
            print(dheader)
            print(dline)
        told = t1
        iter_ += 1
        SSE_old = SSE

        # C (and alpha) update
        XSt = np.dot(X[:, U], S.T)
        C, SSE, muC, mualpha, CtXtXC, XC = C_update(
            X[:, I], XSt, XC, SSt, C, delta, muC, mualpha, SST, SSE, 10
        )

        # S update
        XCtX = np.dot(XC.T, X[:, U])
        S, SSE, muS, SSt = S_update(
            S, XCtX, CtXtXC, muS, SST, SSE, 10
        )

        # Evaluate and display iteration
        dSSE = SSE_old - SSE
        t1 = dt.now()
        if iter_ % 1 == 0:
            time.sleep(0.000001)
            varexpl = (SST - SSE) / SST
            if verbose:
                print('%10.0f | %10.4f | %10.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10.4f \n' % (iter_, varexpl, SSE, dSSE/np.abs(SSE), muC, mualpha, muS, (t1-told).seconds))

    # Display final iteration
    varexpl = (SST - SSE) / SST
    if verbose:
        print(dline)
        print(dline)
        print('%10.0f | %10.4f | %10.4e | %10.4e | %10.4e | %10.4e | %10.4e | %10.4f \n' % (iter_, varexpl, SSE, dSSE/np.abs(SSE), muC, mualpha, muS, (t1-told).seconds))

    # Sort components according to importance
    ind, vals = zip(
        *sorted(enumerate(np.sum(S, axis=1)), key=lambda x: x[0], reverse=1)
    )
    S = S[ind, :]
    C = C[:, ind]
    XC = XC[:, ind]

    return XC, S, C, SSE, varexpl