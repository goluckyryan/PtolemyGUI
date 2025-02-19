#!/usr/bin/env python3

import numpy as np
from scipy.special import gamma

def pochhammer(x, n):
    """Compute Pochhammer symbol (x)_n"""
    if n == 0:
        return 1.0
    result = 1.0
    for i in range(n):
        result *= (x + i)
    return result

def hyp1f1_series(a, b, z, max_terms=1000, tol=1e-10):
    """
    Compute _1F_1(a, b, z) using series expansion.
    
    :param a, b: Parameters of the hypergeometric function
    :param z: Complex argument
    :param max_terms: Maximum number of terms to sum
    :param tol: Tolerance for convergence
    :return: Approximation of _1F_1(a, b, z)
    """
    sum = 0.0 + 0.0j
    last_term = 0.0 + 0.0j
    for n in range(max_terms):
        term = pochhammer(a, n) / pochhammer(b, n) * z**n / np.math.factorial(n)
        sum += term
        if abs(term - last_term) < tol * abs(sum):  # Check for convergence
            return sum
        last_term = term
    return sum  # If we reach here, we've used all terms without converging to desired tolerance

def coulomb_wave_function(L, eta, rho):
    """
    Compute the regular Coulomb wave function F_L(eta, rho).
    
    :param L: Angular momentum quantum number
    :param eta: Sommerfeld parameter
    :param rho: Radial coordinate scaled by k (wavenumber)
    :return: F_L(eta, rho)
    """
    # Compute normalization constant C_L(eta)
    C_L = (2**L * np.exp(-np.pi * eta / 2) * np.abs(gamma(L + 1 + 1j * eta))) / gamma(2*L + 2)
    
    # Compute _1F_1 using our series expansion
    hyp1f1_value = hyp1f1_series(L + 1 - 1j * eta, 2*L + 2, 2j * rho)
    
    # Return the Coulomb wave function
    return C_L * (rho ** (L+1)) * np.exp(-1j * rho) * hyp1f1_value

# Example usage
#L, eta, rho = 1, 0.5, 1.0
#result = coulomb_wave_function(L, eta, rho)
#print(f"F_{L}({eta}, {rho}) = {result}")