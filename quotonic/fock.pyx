# Import Packages & Modules
import itertools as it

import numpy as np

from quotonic.utilities import memoized

# Import Cython Utilities

cimport numpy as np


@memoized
def getDim(int numPhotons, int numModes):
    """Calculate Fock basis dimension.

    Given a number of photons $n$ and a number of optical modes $m$, this function efficiently 
    computes the dimension of the corresponding Fock basis using Cython type definitions. 
    It is also memoized to improve runtime.

    Args:
        numPhotons (int): Number of photons $n$
        numModes (int): Number of optical modes $m$

    Returns:
        int: Fock basis dimension
    """

    # Define all necessary integers
    cdef int i, n, nInit, top, numerator, denominator
    cdef int dim = 0

    # Store the top of {n + m - 1 \choose n}
    top = numPhotons + numModes - 1

    # Evaluate the simplified version of {n + m - 1 \choose n}
    i = 0
    numerator = 1
    denominator = 1
    while top - i >= numModes:
        numerator *= top - i
        i += 1
        denominator *= i
    dim += numerator // denominator
    
    return dim

@memoized
def basis(numPhotons: int, numModes: int) -> list[list[int]]:
    """Generate a list of Fock basis states.

    Given a number of photons $n$ and a number of optical modes $m$, this function creates each
    state in the Fock basis as a list of length $m$ where each element is an integer corresponding
    to the number of photons occupying a given optical mode. All of the states are then placed in a
    list that contains the entire basis.

    Args:
        numPhotons: Number of photons $n$
        numModes: Number of optical modes $m$

    Returns:
        List of states in the Fock basis for $n$ photons and $m$ modes.
    """

    # Initialize list of the Fock basis states
    fockBasis = []

    # Generate a list of tuples of all combinations of the modes for a given number of photons  
    modeBasis = list(it.combinations_with_replacement(range(numModes), numPhotons))
    
    # For each combination of modes, compute the number of photons in each mode 
    # and append a list that is representative of the Fock basis state
    for modeState in modeBasis:
      fockState = []
      for i in range(numModes):
        photonsInMode = modeState.count(i)
        fockState.append(photonsInMode)
      fockBasis.append(fockState)
      
    return fockBasis