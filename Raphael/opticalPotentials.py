#!/usr/bin/env python3

#need to use import

import numpy as np

class OpticalPotential:
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

    def Print(self):
        print(f"Wood-Saxon (Re) : {self.v:7.3f}, {self.r0:5.3f}, {self.a:5.3f}")
        print(f"           (Im) : {self.vi:7.3f}, {self.ri0:5.3f}, {self.ai:5.3f}")
        print(f"WS Surface (Im) : {self.vsi:7.3f}, {self.rsi0:5.3f}, {self.asi:5.3f}")
        print(f"Spin-Orbit (Re) : {self.vso:7.3f}, {self.rso0:5.3f}, {self.aso:5.3f}")
        print(f"           (Im) : {self.vsoi:7.3f}, {self.rsoi0:5.3f}, {self.asoi:5.3f}")
        print(f"Coulomb    (Re) : {'':7s}, {self.rc0:5.3f}")


class SuAndHan(OpticalPotential):
    def __init__(self,  A : int, Z :int , E : float):
        N = A - Z
        A3 = A**(1./3.)

        vsiCOND = 27.5816 - 0.0797 * E + 48.0*(N-Z)/A
        viCOND = -4.0174 + 0.1409 * E 
        
        self.v  = 175.0881 - 0.6236 * E + 0.0006*E*E + 30.*(N-Z)/A - 0.236 * Z/A3
        self.r0 = 1.3421
        self.a  = 0.6578

        self.vi  = viCOND
        if  viCOND < 0 :
            self.vi = 0.0
        self.ri0 = 1.4259
        self.ai  = 0.6578

        self.vsi  = vsiCOND
        if vsiCOND < 0 :
            self.vsi = 0.0
        self.rsi0 = 1.2928
        self.asi  = 0.6359

        self.vso  = 0.0
        self.rso0 = 1.2686
        self.aso  = 0.85

        self.vsoi  = 0.0
        self.rsoi0 = 0.0
        self.asoi  = 0.0

        self.rc0 = 1.350


class AnCai(OpticalPotential):
    def __init__(self,  A : int, Z :int , E : float):

        A3 = A**(1./3.)
        self.v  = 91.85 - 0.249*E + 0.000116*pow(E,2) + 0.642 * Z / A3
        self.r0 = 1.152 - 0.00776 / A3
        self.a  = 0.719 + 0.0126 * A3

        self.vi  = 1.104 + 0.0622 * E
        self.ri0 = 1.305 + 0.0997 / A3
        self.ai  = 0.855 - 0.1 * A3

        self.vsi  = 10.83 - 0.0306 * E
        self.rsi0 = 1.334 + 0.152 / A3
        self.asi  = 0.531 + 0.062 * A3

        self.vso  = 3.557
        self.rso0 = 0.972
        self.aso  = 1.011

        self.vsoi  = 0.0
        self.rsoi0 = 0.0
        self.asoi  = 0.0

        self.rc0 = 1.303

class Koning(OpticalPotential):
    def __init__(self,  A : int, Z :int , E : float, Zproj : float):

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

        self.v  = vp1*(1 - vp2*(E-epf) + vp3*pow(E-epf,2) - vp4*pow(E-epf,3)) + vc * vp1 * (vp2 - 2*vp3*(E-epf) + 3*vp4*pow(E-epf,2))
        #neutron
        if  Zproj == 0 :
            self.v  = vn1*(1 - vn2*(E-enf) + vn3*pow(E-enf,2) - vn4*pow(E-enf,3))

        self.r0 = 1.3039 - 0.4054 / A3
        self.a  = 0.6778 - 0.000148 * A

        self.vi  = wp1 * pow(E-epf,2)/(pow(E-epf,2) + pow(wp2,2))
        if Zproj == 0 :
            self.vi  = wn1 * pow(E-enf,2)/(pow(E-enf,2) + pow(wn2,2))
        
        self.ri0 = 1.3039 - 0.4054 / A3
        self.ai  = 0.6778 - 0.000148 * A

        self.vsi  = dp1 * pow(E-epf,2)/(pow(E-epf,2)+pow(dp3,2)) * np.exp(-dp2*(E-epf))
        if Zproj == 0 :
            self.vsi  = dn1 * pow(E-enf,2)/(pow(E-enf,2)+pow(dn3,2)) * np.exp(-dn2*(E-enf))

        self.rsi0 = 1.3424 - 0.01585 * A3
        self.asi  = 0.5187 + 0.0005205 * A
        if Zproj == 0:
            self.asi = 0.5446 - 0.0001656 * A

        self.vso  = vso1 * np.exp(-vso2 * (E-epf))
        if Zproj == 0:
            self.vso = vso1 * np.exp(-vso2 * (E-enf))

        self.rso0 = 1.1854 - 0.647/A3
        self.aso  = 0.59

        self.vsoi  = wso1 * pow(E-epf,2)/(pow(E-epf,2)+pow(wso2,2))
        if Zproj == 0 :  
            self.vsoi  = wso1 * pow(E-enf,2)/(pow(E-enf,2)+pow(wso2,2))

        self.rsoi0 = 1.1854 - 0.647/A3
        self.asoi  = 0.59

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
