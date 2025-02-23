#!/usr/bin/env python3

#need to use import

import numpy as np

# Woods-Saxon
v = 0
r0 = 0
a = 0
vi = 0
ri0 = 0
ai = 0
# Woods-Saxon Surface 
vsi = 0
rsi0 = 0
asi = 0
# Spin-orbit
vso = 0
rso0 = 0
aso = 0
vsoi = 0
rsoi0 = 0
asoi = 0
# Coulomb 
rc0 = 0

def AnCai(A : int, Z : int, E : float):
    global v, r0, a, vi, ri0, ai, vsi, rsi0, asi, vso, rso0, aso, vsoi, rsoi0, asoi, rc0

    A3 = A**(1./3.)
    v  = 91.85 - 0.249*E + 0.000116*pow(E,2) + 0.642 * Z / A3
    r0 = 1.152 - 0.00776 / A3
    a  = 0.719 + 0.0126 * A3

    vi  = 1.104 + 0.0622 * E
    ri0 = 1.305 + 0.0997 / A3
    ai  = 0.855 - 0.1 * A3

    vsi  = 10.83 - 0.0306 * E
    rsi0 = 1.334 + 0.152 / A3
    asi  = 0.531 + 0.062 * A3

    vso  = 3.557
    rso0 = 0.972
    aso  = 1.011

    vsoi  = 0.0
    rsoi0 = 0.0
    asoi  = 0.0

    rc0 = 1.303

def Koning(A : int, Z : int, E : float, Zproj : float):
    global v, r0, a, vi, ri0, ai, vsi, rsi0, asi, vso, rso0, aso, vsoi, rsoi0, asoi, rc0

    N   = A-Z
    A3 = A**(1./3.)
    
    vp1 = 59.3 + 21.*(N-Z)/A - 0.024*A
    vn1 = 59.3 - 21.*(N-Z)/A - 0.024*A
    
    vp2 = 0.007067 + 0.00000423*A
    vn2 = 0.007228 - 0.00000148*A
    
    vp3 = 0.00001729 + 0.00000001136 * A
    vn3 = 0.00001994 - 0.00000002 * A
    
    vp4 = 7e-9 # = vn4
    vn4 = vp4
    
    wp1 = 14.667 + 0.009629*A
    wn1 = 12.195 + 0.0167*A
    
    wp2 = 73.55 + 0.0795*A # = wn2
    wn2 = wp2
    
    dp1 = 16 + 16.*(N-Z)/A
    dn1 = 16 - 16.*(N-Z)/A
    
    dp2 = 0.018 + 0.003802/(1 + np.exp((A-156.)/8)) # = dn2
    dn2 = dp2
    
    dp3 = 11.5  # = dn3
    dn3 = dp3
    
    vso1 = 5.922 + 0.003 * A
    vso2 = 0.004
    
    wso1 = -3.1
    wso2 = 160
    
    epf = -8.4075 + 0.01378 *A
    enf = -11.2814 + 0.02646 *A
    
    rc = 1.198 + 0.697/pow(A3,2) + 12.995/pow(A3,5)
    vc = 1.73/rc * Z / A3

    v  = vp1*(1 - vp2*(E-epf) + vp3*pow(E-epf,2) - vp4*pow(E-epf,3)) + vc * vp1 * (vp2 - 2*vp3*(E-epf) + 3*vp4*pow(E-epf,2))
    #neutron
    if  Zproj == 0 :
        v  = vn1*(1 - vn2*(E-enf) + vn3*pow(E-enf,2) - vn4*pow(E-enf,3))

    r0 = 1.3039 - 0.4054 / A3
    a  = 0.6778 - 0.000148 * A

    vi  = wp1 * pow(E-epf,2)/(pow(E-epf,2) + pow(wp2,2))
    if Zproj == 0 :
        vi  = wn1 * pow(E-enf,2)/(pow(E-enf,2) + pow(wn2,2))
    
    ri0 = 1.3039 - 0.4054 / A3
    ai  = 0.6778 - 0.000148 * A

    vsi  = dp1 * pow(E-epf,2)/(pow(E-epf,2)+pow(dp3,2)) * np.exp(-dp2*(E-epf))
    if Zproj == 0 :
        vsi  = dn1 * pow(E-enf,2)/(pow(E-enf,2)+pow(dn3,2)) * np.exp(-dn2*(E-enf))

    rsi0 = 1.3424 - 0.01585 * A3
    asi  = 0.5187 + 0.0005205 * A
    if Zproj == 0:
        asi = 0.5446 - 0.0001656 * A

    vso  = vso1 * np.exp(-vso2 * (E-epf))
    if Zproj == 0:
        vso = vso1 * np.exp(-vso2 * (E-enf))

    rso0 = 1.1854 - 0.647/A3
    aso  = 0.59

    vsoi  = wso1 * pow(E-epf,2)/(pow(E-epf,2)+pow(wso2,2))
    if Zproj == 0 :  
        vsoi  = wso1 * pow(E-enf,2)/(pow(E-enf,2)+pow(wso2,2))

    rsoi0 = 1.1854 - 0.647/A3
    asoi  = 0.59

    rc0 = rc

def ConvertLSym(LSym :str) -> int:
    if LSym == "s" :
        return 0
    elif LSym == "p" :
        return 1
    elif LSym == "d" :
        return 2
    elif LSym == "f" :
        return 3
    elif LSym == "g" :
        return 4
    elif LSym == "h" :
        return 5
    elif LSym == "i" :
        return 6
    elif LSym == "j" :
        return 7
    elif LSym == "k" :
        return 8
    else :
        return -1
