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

def KroneckerDelta(i, j):
  if i == j:
    return 1
  else:
    return 0

def obeys_triangle_rule(j1, j2, j3):
    """Check if j1, j2, j3 obey the vector summation rules."""
    # Ensure non-negativity (optional if inputs are guaranteed positive)
    if j1 < 0 or j2 < 0 or j3 < 0:
        return False
    # Triangle inequalities
    if (j3 < abs(j1 - j2) or j3 > j1 + j2):
        return False
    # Check if j1 + j2 + j3 is an integer (for half-integer j, this is automatic)
    if (j1 + j2 + j3) % 1 != 0:
        return False
    return True

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


#============ don;t use, very slow, use the sympy package
def threej(j1, m1, j2, m2, j3, m3):
    if m1 + m2 + m3 != 0:
        return 0
    if obeys_triangle_rule(j1, j2, j3) == False:
        return 0
    
    cg = clebsch_gordan(j1, m1, j2, m2, j3, -m3)
    norm = pow(-1, j1-j2-m3)/(2*j3+1)**0.5
    return norm * cg

def sixj(j1, j2, j3, j4, j5, j6):
    """Compute the 6j symbol using Clebsch-Gordan coefficients."""
    # Check triangle conditions
    if not (obeys_triangle_rule(j1, j2, j3) and
            obeys_triangle_rule(j1, j5, j6) and
            obeys_triangle_rule(j4, j2, j6) and
            obeys_triangle_rule(j4, j5, j3)):
        return 0.0
    
    sixj_value = 0.0
    
    # Ranges for m values
    m1_range = range(-j1, j1 + 1)
    m2_range = range(-j2, j2 + 1)
    m4_range = range(-j4, j4 + 1)
    m5_range = range(-j5, j5 + 1)

    # Sum over m values
    for m1 in m1_range:
        for m2 in m2_range:
            m3 = - m1 - m2
            for m4 in m4_range:
                for m5 in m5_range:
                    m6 = m2 + m4

                    if m3 + m5 not in m4_range or m1 + m6 not in m5_range:
                        continue
                    
                    # cg1 = threej(j1, -m1, j2, -m2, j3, -m3)

                    cg1 = (-1)**(j1-j2+m3) * clebsch_gordan(j1, -m1, j2, -m2, j3, m3) / (2*j3+1)**0.5

                    cg2 = threej(j1,  m1, j5, -m5, j6,  m6)
                    cg3 = threej(j4,  m4, j2,  m2, j6, -m6)
                    cg4 = threej(j4, -m4, j5,  m5, j3,  m3)
                    
                    norm = pow(-1, j1-m1 + j2-m2 + j3-m3 + j4-m4 + j5-m5 + j6-m6)

                    sixj_value += cg1 * cg2 * cg3 * cg4 * norm
    
    return sixj_value

def ninej(j1, j2, j3, j4, j5, j6, j7, j8, j9):
    """Compute the 9j symbol using 6j symbols."""
    # Check triangle conditions for rows
    if not (obeys_triangle_rule(j1, j2, j3) and
            obeys_triangle_rule(j4, j5, j6) and
            obeys_triangle_rule(j7, j8, j9)):
        return 0.0
    
    # Check triangle conditions for columns
    if not (obeys_triangle_rule(j1, j4, j7) and
            obeys_triangle_rule(j2, j5, j8) and
            obeys_triangle_rule(j3, j6, j9)):
        return 0.0
    
    ninej_value = 0.0
    
    # Determine the range of intermediate angular momentum x
    x_min = max(abs(j1 - j9), abs(j4 - j8), abs(j2 - j6))
    x_max = min(j1 + j9, j4 + j8, j2 + j6)
    
    # Sum over x (must be integer or half-integer depending on inputs)
    step = 1 if all(j % 1 == 0 for j in [j1, j2, j3, j4, j5, j6, j7, j8, j9]) else 0.5
    for x in [x_min + i * step for i in range(int((x_max - x_min) / step) + 1)]:
        # if not (obeys_triangle_rule(j1, j4, j7) and
        #         obeys_triangle_rule(j1, j9,  x) and # j1 j9
        #         obeys_triangle_rule(j8, j9, j7) and
        #         obeys_triangle_rule(j8, j4,  x) and # j8 j4
        #         obeys_triangle_rule(j2, j5, j8) and
        #         obeys_triangle_rule(j2,  x, j6) and # j2 j6
        #         obeys_triangle_rule(j4, j5, j6) and
        #         obeys_triangle_rule(j4,  x, j8) and # j4 j8
        #         obeys_triangle_rule(j3, j6, j9) and
        #         obeys_triangle_rule(j3, j1, j2) and
        #         obeys_triangle_rule( x, j6, j2) and # j2 j6
        #         obeys_triangle_rule( x, j1, j9)):   # j1 j9
        #     continue
        
        if not (obeys_triangle_rule(j1, j9,  x) and # j1 j9
                obeys_triangle_rule(j8, j4,  x) and # j8 j4
                obeys_triangle_rule(j2,  x, j6)):   # j1 j9
            continue
        
        sixj1 = sixj(j1, j4, j7, j8, j9,  x)
        sixj2 = sixj(j2, j5, j8, j4,  x, j6)
        sixj3 = sixj(j3, j6, j9,  x, j1, j2)
        
        phase = (-1) ** int(2 * x)  # Phase factor
        weight = 2 * x + 1          # Degeneracy factor
        
        ninej_value += phase * weight * sixj1 * sixj2 * sixj3
    
    return ninej_value