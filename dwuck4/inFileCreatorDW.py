#!/usr/bin/env python3
import sys
import os

if len(sys.argv) < 7:
    print("Error: Not enough arguments provided.")
    print("Usage: ./{sys.argv[0]} reaction target_gs-spin  orbital spin-pi  Ex  ELab[Mev/u]")
    sys.exit(1)

reaction = sys.argv[1]
JA_pi = sys.argv[2]
orbital = sys.argv[3]
JB_pi = sys.argv[4]
Ex = float(sys.argv[5])
ELab = float(sys.argv[6])

sys.path.append(os.path.join(os.path.dirname(__file__), '../Cleopatra'))
sys.path.append(os.path.join(os.path.dirname(__file__), '../Raphael'))
from IAEANuclearData import IsotopeClass
import opticalPotentials as op
from reactionData import ReactionData

#####################################################

#  only for (d,p) or (p,d) using An & Cai, Kronning 

#####################################################
import re

#================== digest reaction

nuclei = re.split(r'[(),]', reaction)

nu_A = nuclei[0]
nu_a = nuclei[1]
nu_b = nuclei[2]
nu_B = nuclei[3]

reactionData = ReactionData(nu_A, nu_a, nu_b, JB_pi, orbital, Ex, ELab)

sym_A = reactionData.sym_A
A_A = reactionData.A_A
Z_A = reactionData.Z_A
A_a = reactionData.A_a
Z_a = reactionData.Z_a
A_B = reactionData.A_B
Z_B = reactionData.Z_B
A_b = reactionData.A_b
Z_b = reactionData.Z_b
A_x = reactionData.A_x
Z_x = reactionData.Z_x
A_c = reactionData.A_c
Z_c = reactionData.Z_c
node = reactionData.node
l_sym = reactionData.l_sym

spin_a = reactionData.spin_a
spin_b = reactionData.spin_b

l = reactionData.l
j = reactionData.j

Q_value = reactionData.Q_value
BindingEnergy = reactionData.BindingEnergy

#=================== outfile name
fileOutName = str(sym_A) + str(A_A) + "_" + str(nu_a) + str(nu_b) + "_" \
           + str(node) + str(l_sym) + str(int(2*j)) + "_" + str(Ex) + "_" + str(ELab) + "_" + JB_pi +".in"


#=================== find the maximum L for partial wave
mass_I = reactionData.mass_I # reduced mass of incoming channel
k_I = reactionData.k_I # wave number of incoming channel
touching_Radius = 1.25*(A_A**(1./3) + A_a**(1./3)) + 10 # add 10 fm 
maxL = int(touching_Radius * k_I) # maximum partial wave

print(f"file out : {fileOutName}")
print(f"   max L : {maxL}")

#=================== create outfile
with open(fileOutName, "w") as file:
    file.write("10001310500100000     " + reaction + "(" + str(Ex) + "," + orbital + ")" +  " @ " + str(ELab) + " MeV/u\n")
    file.write("+181.    +00.    +01.0\n")
    file.write(f"+{maxL}+01+{l:02d}+{int(2*j):02d}\n")
    file.write(f"{0.1:+08.4f}{0.0:+08.4f}{15:+08.4f}\n")
#===== Block 5
    if A_a == 2 :
        pot = op.AnCai(A_A, Z_A, A_a*ELab)
    if A_a == 1 :
        pot = op.Koning(A_A, Z_A, A_a*ELab, Z_a)
    if A_a == 4 :
        pot == op.SuAndHan(A_A, Z_A, A_a*ELab)

    file.write(f"{A_a*ELab:+08.4f}")
    file.write(f"{A_a:+08.4f}")
    file.write(f"{Z_a:+08.4f}")
    file.write(f"{A_A:+08.4f}")
    file.write(f"{Z_A:+08.4f}")
    file.write(f"{pot.rc0:+08.4f}")
    file.write(f"{'':8s}")
    file.write(f"{'':8s}")
    file.write(f"{2*spin_a:+08.4f}\n")
    # Woods-Saxon
    file.write(f"{1:+08.4f}") 
    file.write(f"{-pot.v:+08.4f}") # real
    file.write(f"{pot.r0:+08.4f}") # 
    file.write(f"{pot.a:+08.4f}") # 
    file.write(f"{'':8s}") # spin-orbit skipped
    file.write(f"{-pot.vi:+08.4f}") # imag
    file.write(f"{pot.ri0:+08.4f}") # 
    file.write(f"{pot.ai:+08.4f}\n") # 
    # Woods-Saxon surface
    file.write(f"{2:+08.4f}") 
    file.write(f"{'':8s}") # real
    file.write(f"{'':8s}") # 
    file.write(f"{'':8s}") # 
    file.write(f"{'':8s}") # spin-orbit skipped
    file.write(f"{4*pot.vsi:+08.4f}") # imag
    file.write(f"{pot.rsi0:+08.4f}") # 
    file.write(f"{pot.asi:+08.4f}\n") # 
    # Spin-Orbit 
    file.write(f"{-4:+08.4f}") 
    file.write(f"{-4*pot.vso:+08.4f}") # real
    file.write(f"{pot.rso0:+08.4f}") # 
    file.write(f"{pot.aso:+08.4f}") # 
    file.write(f"{'':8s}") # spin-orbit skipped
    file.write(f"{-4*pot.vsoi:+08.4f}") # imag
    file.write(f"{pot.rsoi0:+08.4f}") # 
    file.write(f"{pot.asoi:+08.4f}\n") # 
#===== Block 6
    Eout = A_a*ELab + Q_value - Ex
    if A_b == 1 :
        pot = op.Koning(A_B, Z_B, Eout, Z_b)
    if A_b == 2 :
        pot = op.AnCai(A_B, Z_B, Eout)
    if A_b == 4 :
        pot = op.SuAndHan(A_B, Z_B, Eout)

    file.write(f"{Q_value:+08.4f}")
    file.write(f"{A_b:+08.4f}")
    file.write(f"{Z_b:+08.4f}")
    file.write(f"{A_B:+08.4f}")
    file.write(f"{Z_B:+08.4f}")
    file.write(f"{pot.rc0:+08.4f}")
    file.write(f"{'':8s}")
    file.write(f"{'':8s}")
    file.write(f"{2*spin_b:+08.4f}\n")
    # Woods-Saxon
    file.write(f"{1:+08.4f}") 
    file.write(f"{-pot.v:+08.4f}") # real
    file.write(f"{pot.r0:+08.4f}") # 
    file.write(f"{pot.a:+08.4f}") # 
    file.write(f"{'':8s}") # spin-orbit skipped
    file.write(f"{-pot.vi:+08.4f}") # imag
    file.write(f"{pot.ri0:+08.4f}") # 
    file.write(f"{pot.ai:+08.4f}\n") # 
    # Woods-Saxon surface
    file.write(f"{2:+08.4f}") 
    file.write(f"{'':8s}") # real
    file.write(f"{'':8s}") # 
    file.write(f"{'':8s}") # 
    file.write(f"{'':8s}") # spin-orbit skipped
    file.write(f"{4*pot.vsi:+08.4f}") # imag
    file.write(f"{pot.rsi0:+08.4f}") # 
    file.write(f"{pot.asi:+08.4f}\n") # 
    # Spin-Orbit 
    file.write(f"{-4:+08.4f}") 
    file.write(f"{-4*pot.vso:+08.4f}") # real
    file.write(f"{pot.rso0:+08.4f}") # 
    file.write(f"{pot.aso:+08.4f}") # 
    file.write(f"{'':8s}") # spin-orbit skipped
    file.write(f"{-4*pot.vsoi:+08.4f}") # imag
    file.write(f"{pot.rsoi0:+08.4f}") # 
    file.write(f"{pot.asoi:+08.4f}\n") # 
#====== bound state
    file.write(f"{BindingEnergy:+08.4f}")
    file.write(f"{A_x:+08.4f}")
    file.write(f"{Z_x:+08.4f}")
    file.write(f"{A_c:+08.4f}")
    file.write(f"{Z_c:+08.4f}")
    file.write(f"{1.30:+08.4f}") # Coulomb radius
    file.write(f"{'':8s}") # 
    file.write(f"{'':8s}") # 
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