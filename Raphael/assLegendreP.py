#!/usr/bin/env python3

import numpy as np

def associated_legendre_array(maxL, maxM, theta_deg):
    # Convert theta from degrees to radians
    theta = np.radians(theta_deg)
    x = np.cos(theta)
    
    P = np.zeros((maxL + 1, maxM + 1))
    
    # P^m_l for m = 0, l = 0
    P[0, 0] = 1.0
    
    # P^m_l for m = 0, l > 0
    for l in range(1, maxL + 1):
        P[l, 0] = ((2*l - 1) * x * P[l-1, 0] - (l - 1) * P[l-2, 0]) / l
    
    # P^m_l for m > 0 (using recursion)
    for m in range(1, maxM + 1):
        # P^m_m
        P[m, m] = (1 - 2*m) * P[m-1, m-1]
        P[m, m] *= np.sqrt(1 - x**2)
        
        for l in range(m + 1, maxL + 1):
            # P^m_l for m < l
            P[l, m] = ((2*l - 1) * x * P[l-1, m] - (l + m - 1) * P[l-2, m]) / (l - m)
    
    return P

# Example usage
# L = 15  # Maximum l degree
# M = 3  # Maximum m order

# legendre_polynomials = associated_legendre_array(L, M, 10)

# print(legendre_polynomials)