# Import Python Libraries & Modules
import numpy as np

import quotonic.fock as fock
from quotonic.utilities import memoized

# Import Cython Utilities

cimport cython
cimport numpy as np
from cython.parallel cimport parallel, prange
from libc.stdlib cimport abort, free, malloc
from libc.string cimport memset


@memoized
def factorial(int x):
    # Initialize the result of the factorial
    cdef long long result = 1

    # Multiply by each integer up to x
    for i in range(2, x + 1):
        result = result * i
    return result

def fockState_to_inds(int numPhotons, np.ndarray[np.int_t, ndim=1] state):
    # Initialize array of matrix indices to build U_T or U_{S,T}
    cdef np.ndarray inds = np.zeros([numPhotons, ], dtype=np.int)
    cdef int mode
    cdef int s = 0

    # Iterate through the modes of the Fock basis state
    for i in range(state.shape[0]):
        mode = state[i]
        if mode == 0:
            continue

        # Insert index i as many times as there are photons in mode i
        for j in range(mode):
            inds[s] = i
            s += 1
    
    return inds

@memoized
def prepNorms_and_USTinds(int numPhotons, int numModes, int fockDim):
    # Retrieve a list of all Fock basis states
    cdef np.ndarray fockBasis
    fockBasis = np.array(fock.basis(numPhotons, numModes), dtype=np.int)

    # Define required Numpy arrays
    cdef np.ndarray[np.double_t, ndim=1] factorialProducts
    cdef np.ndarray[np.double_t, ndim=2] norms
    cdef np.ndarray[np.int_t, ndim=2] inds
    
    # Initialize all Numpy arrays with zeroes
    factorialProducts = np.zeros([fockDim, ], dtype=np.double)
    norms = np.zeros([fockDim, fockDim], dtype=np.double)
    inds = np.zeros([fockDim, numPhotons], dtype=np.int)

    # For each basis state, compute factorial product and obtain indices for transformation
    cdef int i
    for i in range(fockDim):
        state = fockBasis[i]
        factorialProduct = 1

        # Compute product of factorials for each mode in the basis state
        for s in state:
            factorialProduct *= factorial(s)
        # Take the square root of the factorial product and store
        factorialProducts[i] = np.sqrt(factorialProduct)

        # Convert basis state to indices required for transformation
        inds[i] = fockState_to_inds(numPhotons, state)

    # Compute normalization constants for each element of \Phi(U)
    norms = 1.0 / np.outer(factorialProducts, factorialProducts)
    
    return (norms, inds)

@memoized
def prepGrayCode(int numPhotons):
    # Define required Numpy arrays and integers
    cdef np.ndarray[np.int_t, ndim = 1] grayInds
    cdef np.ndarray[np.int_t, ndim = 1] graySgns
    cdef int k, currGrayVal, prevGrayVal, xorGrayVals, grayInd

    # Initialize Gray code indices and signum results as empty
    grayInds = np.empty([2 ** numPhotons, ], dtype=np.int)
    graySgns = np.empty([2 ** numPhotons, ], dtype=np.int)

    # The initial Gray code index and signum result are always 0 and 1 respectively
    grayInds[0] = 0
    graySgns[0] = 1

    # Compute all Gray code indices and signum results in parallel
    with nogil, parallel():
        for k in prange(1, 2 ** numPhotons, schedule='dynamic'):
            # Must compute previous Gray code value since loop is in parallel
            prevGrayVal = (k-1) ^ ((k-1) >> 1)
            currGrayVal = k ^ (k >> 1)

            # Compute the Gray code signum result for this iteration
            if currGrayVal < prevGrayVal:
                graySgns[k] = -1
            else:
                graySgns[k] = 1

            # Compute the Gray code index for this iteration
            xorGrayVals = currGrayVal ^ prevGrayVal
            grayInd = 0
            while (xorGrayVals & 1) == 0:
                xorGrayVals = xorGrayVals >> 1
                grayInd = grayInd + 1
            grayInds[k] = grayInd

    return (grayInds, graySgns)

@cython.boundscheck(False)
@cython.wraparound(False)
def multiPhotonUnitary(int numPhotons, np.ndarray[np.complex128_t, ndim=2] U):
    # If there is only one photon, this transformation is not required
    if numPhotons == 1:
        return U

    # Define integers and get the number of optical modes, Fock basis dimension
    cdef int numModes, fockDim, ryserSgn
    numModes = U.shape[0]
    fockDim = fock.getDim(numPhotons, numModes)

    # Compute (-1)^n as necessary for Ryser's algorithm
    if numPhotons % 2 == 0:
        ryserSgn = 1
    else:
        ryserSgn = -1

    # Define and retrieve normalization constants, indices required for $U_T$, $U_{S,T}$ preparation
    cdef np.ndarray[np.double_t, ndim = 2] norms
    cdef np.ndarray[np.int_t, ndim = 2] inds
    norms, inds = prepNorms_and_USTinds(numPhotons, numModes, fockDim)

    # Define and retrieve Gray indices, Gray signums to use Gray code with Ryser's algorithm
    cdef np.ndarray[np.int_t, ndim= 1] grayInds
    cdef np.ndarray[np.int_t, ndim= 1] graySgns
    cdef int grayInd, graySgn
    grayInds, graySgns = prepGrayCode(numPhotons)

    # Define and intialize multi-photon unitary \Phi(U) with empty elements
    cdef np.ndarray[np.complex128_t, ndim = 2] PhiU
    PhiU = np.empty([fockDim, fockDim], dtype=np.complex128)

    # Define variables to store $U_T$, $U_{S,T}$ and terms of Ryser's formula
    cdef complex * U_T
    cdef complex * U_ST
    cdef complex * termSums
    cdef complex termProduct, perm

    # Define integers used for all loops
    cdef int row, col, i, j, k, I, J

    # Compute all elements of \Phi(U) in parallel using 12 threads
    with nogil, parallel(num_threads=12):
        
        # Allocate memory for matrix U_T (stored as array)
        U_T = <complex *> malloc(sizeof(complex) * numModes * numPhotons)
        if U_T == NULL:
            abort()

        # Allocate memory for matrix U_{S,T} (stored as array)
        U_ST = <complex *> malloc(sizeof(complex) * numPhotons * numPhotons)
        if U_ST == NULL:
            abort()

        # Allocate memory for list of elements to sum for a term of Ryser's formula
        termSums = <complex *> malloc(sizeof(complex) * numPhotons)
        if termSums == NULL:
            abort()

        # Separate columns of \Phi(U) to be computed in parallel
        for col in prange(fockDim, schedule='dynamic'):
            # Construct matrix U_T for the given column of \Phi(U)
            for j in range(numPhotons):
                J = inds[col, j]
                for i in range(numModes):
                    U_T[i + j * numModes] = U[i, J]
            
            # Iterate through all \Phi(U) elements in the given column
            for row in range(fockDim):
                # Construct matrix U_{S,T} for given element of \Phi(U)
                for i in range(numPhotons):
                    I = inds[row, i]
                    for j in range(numPhotons):
                        U_ST[i + j * numPhotons] = U_T[I + j * numModes]

                # Initialize the permanent and all elements of termSums
                perm = 0
                memset(termSums, 0, sizeof(complex) * numPhotons)
                for k in range(1, 2 ** numPhotons):
                    # Obtain current Gray index and Gray signum result for faster use
                    grayInd = grayInds[k]
                    graySgn = graySgns[k]

                    # Update all elements of termSums in accordance with Ryser's formula
                    for i in range(numPhotons):
                        termSums[i] = termSums[i] + graySgn * U_ST[i + grayInd * numPhotons]

                    # Compute (-1)^{|X|} as necessary for Ryser's formula
                    termProduct = 1 - 2 * (k % 2)

                    # Compute product in Ryser's formula for this iteration
                    for i in range(numPhotons):
                        termProduct = termProduct * termSums[i]

                    # Add product for this iteration of Ryser's formula
                    perm = perm + termProduct

                # Complete permanent calculation, multiply by normalization constant, and store element
                PhiU[row, col] = ryserSgn * norms[row, col] * perm

        # Free memory allocated for U_T, U_{S,T} and termSums respectively
        free(U_T)
        free(U_ST)
        free(termSums)

    return PhiU