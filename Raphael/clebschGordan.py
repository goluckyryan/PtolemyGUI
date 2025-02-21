#!/usr/bin/env python3

import numpy as np
from scipy.special import gamma

# from sympy.physics.quantum.cg import CG
# from sympy import S
# def clebsch_gordan(j1, m1, j2, m2, j, m):
#     cg = CG(S(j1), S(m1), S(j2), S(m2), S(j), S(m))
#     result = cg.doit()
#     return np.complex128(result)

import numpy as np
from math import sqrt

def quantum_factorial(n):
    """
    Calculate factorial for integer or half-integer numbers using gamma function.
    For integer n: n! = n * (n-1) * ... * 1
    For half-integer n: n! = Γ(n + 1)
    """
    if n < 0:
        return 0.0
    return gamma(n + 1)

def clebsch_gordan(j1, m1,j2, m2, j, m):
    """
    Calculate Clebsch-Gordan coefficient <j1 m1 j2 m2 | j m>
    
    Parameters:
    j1, j2: angular momentum quantum numbers
    m1, m2: magnetic quantum numbers
    j: total angular momentum quantum number
    m: total magnetic quantum number
    
    Returns:
    float: Clebsch-Gordan coefficient value
    """
    # Check validity of inputs using triangular inequalities and conservation
    if not np.isclose(m, m1 + m2, atol=1e-10):
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(m) > j:
        return 0.0
    if not (abs(j1 - j2) <= j <= j1 + j2):
        return 0.0
    if j1 < 0 or j2 < 0 or j < 0:
        return 0.0
    
    # Ensure all quantum numbers are either integer or half-integer
    if not (np.mod(2*j1, 1) < 1e-10 or np.isclose(np.mod(2*j1, 1), 1, atol=1e-10)):
        return 0.0
    if not (np.mod(2*j2, 1) < 1e-10 or np.isclose(np.mod(2*j2, 1), 1, atol=1e-10)):
        return 0.0
    if not (np.mod(2*j, 1) < 1e-10 or np.isclose(np.mod(2*j, 1), 1, atol=1e-10)):
        return 0.0

    # Calculate the coefficient
    prefactor = sqrt((2*j + 1) * quantum_factorial(j1 + j2 - j) * 
                    quantum_factorial(j1 - j2 + j) * 
                    quantum_factorial(-j1 + j2 + j) / 
                    quantum_factorial(j1 + j2 + j + 1))
    
    prefactor *= sqrt(quantum_factorial(j + m) * quantum_factorial(j - m) * 
                     quantum_factorial(j1 - m1) * quantum_factorial(j1 + m1) * 
                     quantum_factorial(j2 - m2) * quantum_factorial(j2 + m2))
    
    # Sum over k
    sum_result = 0.0
    k_min = max(0, max(j2 - j - m1, j1 + m2 - j))
    k_max = min(j1 + j2 - j, min(j1 - m1, j2 + m2))
    
    # Ensure k_min and k_max are integers for the range
    k_min = int(np.ceil(k_min))
    k_max = int(np.floor(k_max))
    
    for k in range(k_min, k_max + 1):
        denominator = (quantum_factorial(k) * 
                      quantum_factorial(j1 + j2 - j - k) * 
                      quantum_factorial(j1 - m1 - k) * 
                      quantum_factorial(j2 + m2 - k) * 
                      quantum_factorial(j - j2 + m1 + k) * 
                      quantum_factorial(j - j1 - m2 + k))
        if np.isclose(denominator, 0, atol=1e-10):
            continue
        term = (-1)**k / denominator
        sum_result += term
    
    return prefactor * sum_result

