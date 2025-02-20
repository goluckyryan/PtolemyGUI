#!/usr/bin/env python3
import sys
import os

reaction = sys.argv[1]
JA_pi = sys.argv[2]
orbital = sys.argv[3]
JB_pi = sys.argv[4]
Ex = float(sys.argv[5])
ELab = float(sys.argv[6])

if len(sys.argv) < 7:
    print("Error: Not enough arguments provided.")
    print("Usage: ./{sys.argv[0]} reaction target_gs-spin  orbital spin-pi  Ex  ELab[Mev/u]")
    sys.exit(1)

sys.path.append(os.path.join(os.path.dirname(__file__), '../Cleopatra'))
from IAEANuclearData import IsotopeClass

#####################################################

#  only for (d,p) or (p,d) using An & Cai, Kronning 

#####################################################
import numpy as np
import re
import matplotlib.pyplot as plt

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

#================== digest reaction

nuclei = re.split(r'[(),]', reaction)

nu_A = nuclei[0]
nu_a = nuclei[1]
nu_b = nuclei[2]
nu_B = nuclei[3]

iso = IsotopeClass()

A_A, Z_A = iso.GetAZ(nu_A)
A_a, Z_a = iso.GetAZ(nu_a)
A_b, Z_b = iso.GetAZ(nu_b)
A_B, Z_B = iso.GetAZ(nu_B)

A_x = abs(A_a - A_b)
Z_x = abs(Z_a - Z_b)

#---- check mass number and charge number is balnaced
if A_A + A_a - A_b - A_B != 0 or Z_A + Z_a - Z_b - Z_B != 0 :
    print("reaction is incorrect, mass or charge not balanced.")
    exit()

#---- check is (d,p) or (p, d)
if (Z_a !=1 or Z_b != 1) or (A_a + A_b != 3) :
    print("not (d,p) or (p,d) reaction. stop.")
    exit()

mass_A = iso.GetMassFromSym(nu_A)
mass_a = iso.GetMassFromSym(nu_a)
mass_b = iso.GetMassFromSym(nu_b)
mass_B = iso.GetMassFromSym(nu_B)

mass_x = iso.GetMassFromAZ( A_x, Z_x)

#.... core 
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

#=================== outfile name
fileOutName = str(sym_A) + str(A_A) + "_" + str(nu_a) + str(nu_b) + "_" \
           + str(node) + str(l_sym) + str(int(2*j)) + "_" + str(Ex) + "_" + str(ELab) + ".in"

print(fileOutName)

#=================== find the maximum L for partial wave
mass_I = mass_A * mass_a / (mass_A + mass_a) # reduced mass of incoming channel
hbarc = 197.3269788 # MeV.fm
k_I = np.sqrt(2*mass_I * A_a * ELab)/hbarc # wave number of incoming channel
touching_Radius = 1.25*(A_A**(1./3) + A_a**(1./3)) + 10 # add 10 fm 
maxL = int(touching_Radius * k_I) # maximum partial wave
print(f"max L : {maxL}")

#=================== create outfile
with open(fileOutName, "w") as file:
    file.write("10001310500100000     " + reaction + "(" + str(Ex) + "," + orbital + ")" +  " @ " + str(ELab) + " MeV/u\n")
    file.write("+181.    +00.    +01.0\n")
    file.write(f"+{maxL}+01+{l:02d}+{int(2*j):02d}\n")
    file.write(f"{0.1:+08.4f}{15:+08.4f}\n")
#===== Block 5
    if A_a == 2 :
        AnCai(A_A, Z_A, A_a*ELab)
    else:
        Koning(A_A, Z_A, A_a*ELab, Z_a)

    file.write(f"{A_a*ELab:+08.4f}")
    file.write(f"{A_a:+08.4f}")
    file.write(f"{Z_a:+08.4f}")
    file.write(f"{A_A:+08.4f}")
    file.write(f"{Z_A:+08.4f}")
    file.write(f"{rc0:+08.4f}")
    file.write(f"{"":8s}")
    file.write(f"{"":8s}")
    file.write(f"{2*spin_a:+08.4f}\n")
    # Woods-Saxon
    file.write(f"{1:+08.4f}") 
    file.write(f"{-v:+08.4f}") # real
    file.write(f"{r0:+08.4f}") # 
    file.write(f"{a:+08.4f}") # 
    file.write(f"{"":8s}") # spin-orbit skipped
    file.write(f"{-vi:+08.4f}") # imag
    file.write(f"{ri0:+08.4f}") # 
    file.write(f"{ai:+08.4f}\n") # 
    # Woods-Saxon surface
    file.write(f"{2:+08.4f}") 
    file.write(f"{"":8s}") # real
    file.write(f"{"":8s}") # 
    file.write(f"{"":8s}") # 
    file.write(f"{"":8s}") # spin-orbit skipped
    file.write(f"{4*vsi:+08.4f}") # imag
    file.write(f"{rsi0:+08.4f}") # 
    file.write(f"{asi:+08.4f}\n") # 
    # Spin-Orbit 
    file.write(f"{-4:+08.4f}") 
    file.write(f"{-4*vso:+08.4f}") # real
    file.write(f"{rso0:+08.4f}") # 
    file.write(f"{aso:+08.4f}") # 
    file.write(f"{"":8s}") # spin-orbit skipped
    file.write(f"{-4*vsoi:+08.4f}") # imag
    file.write(f"{rsoi0:+08.4f}") # 
    file.write(f"{asoi:+08.4f}\n") # 
#===== Block 6
    if A_a == 2 :
        Koning(A_B, Z_B, A_a*ELab + Q_value - Ex, Z_b)
    else:
        AnCai(A_B, Z_B, A_a*ELab + Q_value - Ex)
    file.write(f"{Q_value:+08.4f}")
    file.write(f"{A_b:+08.4f}")
    file.write(f"{Z_b:+08.4f}")
    file.write(f"{A_B:+08.4f}")
    file.write(f"{Z_B:+08.4f}")
    file.write(f"{rc0:+08.4f}")
    file.write(f"{"":8s}")
    file.write(f"{"":8s}")
    file.write(f"{2*spin_b:+08.4f}\n")
    # Woods-Saxon
    file.write(f"{1:+08.4f}") 
    file.write(f"{-v:+08.4f}") # real
    file.write(f"{r0:+08.4f}") # 
    file.write(f"{a:+08.4f}") # 
    file.write(f"{"":8s}") # spin-orbit skipped
    file.write(f"{-vi:+08.4f}") # imag
    file.write(f"{ri0:+08.4f}") # 
    file.write(f"{ai:+08.4f}\n") # 
    # Woods-Saxon surface
    file.write(f"{2:+08.4f}") 
    file.write(f"{"":8s}") # real
    file.write(f"{"":8s}") # 
    file.write(f"{"":8s}") # 
    file.write(f"{"":8s}") # spin-orbit skipped
    file.write(f"{4*vsi:+08.4f}") # imag
    file.write(f"{rsi0:+08.4f}") # 
    file.write(f"{asi:+08.4f}\n") # 
    # Spin-Orbit 
    file.write(f"{-4:+08.4f}") 
    file.write(f"{-4*vso:+08.4f}") # real
    file.write(f"{rso0:+08.4f}") # 
    file.write(f"{aso:+08.4f}") # 
    file.write(f"{"":8s}") # spin-orbit skipped
    file.write(f"{-4*vsoi:+08.4f}") # imag
    file.write(f"{rsoi0:+08.4f}") # 
    file.write(f"{asoi:+08.4f}\n") # 
#====== bound state
    file.write(f"{BindingEnergy:+08.4f}")
    file.write(f"{A_x:+08.4f}")
    file.write(f"{Z_x:+08.4f}")
    file.write(f"{A_c:+08.4f}")
    file.write(f"{Z_c:+08.4f}")
    file.write(f"{1.30:+08.4f}") # Coulomb radius
    file.write(f"{"":8s}") # 
    file.write(f"{"":8s}") # 
    file.write(f"{1:+08.4f}\n") # neutron spin x 2
    # Woods-Saxon
    file.write(f"{-1:+08.4f}") 
    file.write(f"{-1:+08.4f}") # real
    file.write(f"{1.1:+08.4f}") # 
    file.write(f"{0.65:+08.4f}") # 
    file.write(f"{24:+8.4f}\n") # spin-orbit 
    # orbital
    file.write(f"{node:+08.4f}") 
    file.write(f"{l:+08.4f}") 
    file.write(f"{2*j:+08.4f}") 
    file.write(f"{1:+08.4f}")  # 2 x nuetron spin
    file.write(f"{58:+08.4f}\n") 
#======== end of input
    file.write("9                           end of input card")