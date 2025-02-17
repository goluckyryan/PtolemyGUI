#!/usr/bin/env python3

import numpy as np
from scipy.special import lpmv

def lgndr(mplus, lplus, thet):
    """
    Calculates Legendre polynomials Plm
    
    Parameters:
    mplus : int
        Number of m's > 0
    lplus : int
        Number of l's > 0
    thet : float
        Angle in degrees
    
    Returns:
    plm : list
        List containing Legendre polynomials
    """

    theta = np.radians(thet)
    y = np.cos(theta)
    z = np.sin(theta)
    plm = np.zeros(459, dtype=np.float64)
    
    ix = 0
    for m in range(1, mplus + 1):  # For MPLUS = 1, LPLUS = 16
        lx = m - 1  # LX = 0
        l2 = 0  # L2 = 0
        p3 = 1.0  # P3 = 1.0
        fl1 = float(lx)  # FL1 = 0
        
        if lx != 0:
            for lt in range(1, lx + 1):
                fl1 += 1.0
                p3 *= fl1 * z / 2.0
        
        p2 = 0.0  # P2 = 0.0
        fl2 = fl1 + 1.0  # FL2 = 1.0
        fl3 = 1.0  # FL3 = 1.0
        
        for lt in range(1, lplus + 1):  # Loop Lb
            ix1 = ix + lt
            
            if l2 < lx:
                plm[ix1] = 0.0
            else:
                if l2 > lx:
                    p3 = (fl2 * y * p2 - fl1 * p1) / fl3
                    fl1 += 1.0
                    fl2 += 2.0
                    fl3 += 1.0
                plm[ix1] = p3
                print(f'PLM, {lx:3d}, {l2:3d}, {ix1:3d}, {plm[ix1]:15.10f}')
                p1, p2 = p2, p3
            
            l2 += 1
        
        ix += lplus
    
    return plm


plm = lgndr(3, 16, 1)