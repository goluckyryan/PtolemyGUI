#!/usr/bin/env python3


from boundState import BoundState
from solveSE import WoodsSaxonPot, CoulombPotential, SpinOrbit_Pot, WS_SurfacePot
import matplotlib.pyplot as plt

boundState = BoundState(16, 8, 1, 0, 1, 0, 0.5, -3.273)
boundState.SetPotential(1.25, 0.65, -6, 1.10, 0.65, 1.30)
boundState.FindPotentialDepth(-75, -40, 0.1)
# boundState.PrintWF()
boundState.PlotBoundState()

exit()

from distortedWave import DistortedWave

# dw = DistortedWave("60Ni", "p", 30)
# dw.ClearPotential()
# dw.AddPotential(WoodsSaxonPot(-47.937-2.853j, 1.20, 0.669), False)
# dw.AddPotential(WS_SurfacePot(-6.878j, 1.28, 0.550), False)
# dw.AddPotential(SpinOrbit_Pot(-5.250 + 0.162j, 1.02, 0.590), False)
# dw.AddPotential(CoulombPotential(1.258), False)

# dw = DistortedWave("60Ni", "d", 60)
# dw.PrintInput()
# dw.ClearPotential()
# dw.AddPotential(WoodsSaxonPot(-81.919, 1.15, 0.768), False)
# dw.AddPotential(WoodsSaxonPot(-4.836j, 1.33, 0.464), False)
# dw.AddPotential(WS_SurfacePot(-8.994j, 1.373, 0.774), False)
# dw.AddPotential(SpinOrbit_Pot(-3.557, 0.972, 1.011), False)
# dw.AddPotential(CoulombPotential(1.303), False)


# dw.CalScatteringMatrix()
# dw.PrintScatteringMatrix()

# dw.PlotDCSUnpolarized(180, 1)

# for i in range(1, 19):
#   theta = 10*i
#   # ruth = dw.RutherFord(theta)
#   # coulAmp = dw.CoulombScatterintAmp(theta)
#   dw.CalLegendre(theta)
#   nuAmp1 = dw.NuclearScatteringAmp(-0.5, 0.5, 14)
#   nuAmp2 = dw.NuclearScatteringAmp(0.5, -0.5, 14)
#   # dsc = dw.DCSUnpolarized(theta, 14)
#   # print(f"{theta:3.0f}, {nuAmp1:15.5f}, {nuAmp2:15.5f}, {dsc:10.6f}, {ruth:10.6f}")
#   print(f"{theta:3.0f}, {nuAmp1:15.5f}, {nuAmp2:15.5f}")


import sys, os
import re
import numpy as np

sys.path.append(os.path.join(os.path.dirname(__file__), '../Cleopatra'))
from IAEANuclearData import IsotopeClass

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

#==========================================

nu_A = "16O"
nu_a = "d"
nu_b = "p"
nu_B = "17O"
ELabPreU = 10 # MeV/u
Ex = 0.87
J_B = 0.5
orbital = "1s1/2"

iso = IsotopeClass()

A_A, Z_A = iso.GetAZ(nu_A)
A_a, Z_a = iso.GetAZ(nu_a)
A_b, Z_b = iso.GetAZ(nu_b)
A_B, Z_B = iso.GetAZ(nu_B)

A_x = abs(A_a - A_b)
Z_x = abs(Z_a - Z_b)

mass_A = iso.GetMassFromSym(nu_A)
mass_a = iso.GetMassFromSym(nu_a)
mass_b = iso.GetMassFromSym(nu_b)
mass_B = iso.GetMassFromSym(nu_B)

mass_x = iso.GetMassFromAZ( A_x, Z_x)

if A_A < A_B : # (d,p)
    A_c = A_A
    Z_c = Z_A
    BindingEnergy = mass_B - mass_A - mass_x + Ex
else:  #(p,d)
    A_c = A_B
    Z_c = Z_B
    BindingEnergy = mass_A - mass_B - mass_x

sym_A = iso.GetSymbol(A_A, Z_A)
sym_B = iso.GetSymbol(A_B, Z_B)

if A_a == 2 and Z_a == 1:
    spin_a = 1.0
    spin_b = 0.5   
else:
    spin_a = 0.5
    spin_b = 1.0

Q_value = mass_A + mass_a - mass_b - mass_B - Ex

print(f"Q-value : {Q_value:10.6f} MeV")
print(f"Binding : {BindingEnergy:10.6f} MeV")

#=================== digest orbital
match = re.search(r'[a-zA-Z]', orbital)  # Find first letter
if match:
    index = match.start()  # Get position of the first letter
    
node = int(orbital[:index])
l_sym = orbital[index:index+1]
j_sym = orbital[index+1:]
j = eval(j_sym)
l = ConvertLSym(l_sym)

#=================== find the maximum L for partial wave
mass_I = mass_A * mass_a / (mass_A + mass_a) # reduced mass of incoming channel
hbarc = 197.3269788 # MeV.fm
k_I = np.sqrt(2*mass_I * A_a * ELabPreU)/hbarc # wave number of incoming channel
touching_Radius = 1.25*(A_A**(1./3) + A_a**(1./3)) + 10 # add 10 fm 
maxL = int(touching_Radius * k_I) # maximum partial wave
print(f"max L : {maxL}")

#================== Bound state
boundState = BoundState(A_c, Z_c, A_x, Z_x, node, l, j, BindingEnergy)
boundState.SetPotential(1.25, 0.65, -6, 1.10, 0.65, 1.30)
boundState.FindPotentialDepth(-70, -55, 0.1)
# # boundState.PrintWF()
# boundState.PlotBoundState()

#================== incoming wave function
AnCai(A_A, Z_A, A_a * ELabPreU)

dwI = DistortedWave(nu_A, nu_a, ELabPreU * A_a)
dwI.maxL = maxL
dwI.ClearPotential()
dwI.AddPotential(WoodsSaxonPot(-v, r0, a), False)
dwI.AddPotential(WoodsSaxonPot(-1j*vi, ri0, ai), False)
dwI.AddPotential(WS_SurfacePot(-1j*vsi, rsi0, asi), False)
dwI.AddPotential(SpinOrbit_Pot(-vso , rso0, aso), False)
dwI.AddPotential(SpinOrbit_Pot(- 1j* vsoi, rsoi0, asoi), False)
dwI.AddPotential(CoulombPotential(rc0), False)

sm_I, wfu_I = dwI.CalScatteringMatrix()

dwI.PrintScatteringMatrix()

#================= outgoing wave function
Koning(A_B, Z_B, A_a*ELabPreU + Q_value - Ex, Z_b)

dwO = DistortedWave(nu_B, nu_b, ELabPreU * A_a + Q_value - Ex)
dwO.maxL = maxL
dwO.ClearPotential()
dwO.AddPotential(WoodsSaxonPot(-v, r0, a), False)
dwO.AddPotential(WoodsSaxonPot(-1j*vi, ri0, ai), False)
dwO.AddPotential(WS_SurfacePot(-1j*vsi, rsi0, asi), False)
dwO.AddPotential(SpinOrbit_Pot(-vso , rso0, aso), False)
dwO.AddPotential(SpinOrbit_Pot(- 1j* vsoi, rsoi0, asoi), False)
dwO.AddPotential(CoulombPotential(rc0), False)

sm_O, wfu_O = dwO.CalScatteringMatrix()

dwO.PrintScatteringMatrix()