#!/usr/bin/env python3
import re
import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '../Cleopatra'))
from IAEANuclearData import IsotopeClass

from clebschGordan import  obeys_triangle_rule
from distortedWave import DistortedWave
import opticalPotentials as op

  
def approximate_to_half_integer(value):
  return round(value * 2) / 2

class ReactionData:
  def __init__(self, nu_A:str, nu_a:str, nu_b:str, JB:str, orbital:str, ExB:float, ELabPerU:float, JA:str = None):
    self.SpinBalanced = self.ReactionDigestion(nu_A, nu_a, nu_b, JB, orbital, ExB, ELabPerU, JA)
    
  def ReactionDigestion(self, nu_A:str, nu_a:str, nu_b:str, JB:str, orbital:str, ExB:float, ELabPerU:float, JA):
    iso = IsotopeClass()
      
    self.A_A, self.Z_A = iso.GetAZ(nu_A)
    self.A_a, self.Z_a = iso.GetAZ(nu_a)
    self.A_b, self.Z_b = iso.GetAZ(nu_b)

    self.A_B = self.A_A + self.A_a - self.A_b
    self.Z_B = self.Z_A + self.Z_a - self.Z_b

    self.ELab = self.A_a * ELabPerU

    mass_A = iso.GetMassFromSym(nu_A)
    mass_a = iso.GetMassFromSym(nu_a)
    mass_b = iso.GetMassFromSym(nu_b)
    mass_B = iso.GetMassFromAZ(self.A_B, self.Z_B)
    ExB = ExB

    self.sym_A = iso.GetSymbol(self.A_A, self.Z_A)
    self.sym_B = iso.GetSymbol(self.A_B, self.Z_B)

    nu_B = f"{self.A_B}{self.sym_B}"
    # print(nu_B)

    if JA is not None:
      spin_A_str = JA
    else:
      spin_A_str = iso.GetJpi(self.A_A, self.Z_A)

    self.spin_A = float(eval(re.sub(r'[+-]', '', spin_A_str)))
    self.spin_B = float(eval(re.sub(r'[+-]', '', JB)))

    # print("-------- spin_B",self.spin_B)

    if self.A_a == 2 and self.Z_a == 1:
        self.spin_a = 1.0
        self.spin_b = 0.5
    else:
        self.spin_a = 0.5
        self.spin_b = 1.0

    #====== transfering nucleon
    self.s = 1/2 # spin of x, neutron or proton
    self.A_x = abs(self.A_a - self.A_b)
    self.Z_x = abs(self.Z_a - self.Z_b)

    mass_x = iso.GetMassFromAZ(self.A_x, self.Z_x)

    #======== core
    if self.A_A < self.A_B : # (d,p)
        self.A_c = self.A_A
        self.Z_c = self.Z_A
        self.BindingEnergy = mass_B - mass_A - mass_x + ExB
    else:  #(p,d)
        self.A_c = self.A_B
        self.Z_c = self.Z_B
        self.BindingEnergy = mass_A - mass_B - ExB - mass_x

    #=================== digest orbital
    match = re.search(r'[a-zA-Z]', orbital)  # Find first letter
    if match:
        index = match.start()  # Get position of the first letter
        
    self.node = int(orbital[:index])
    self.l_sym = orbital[index:index+1]
    j_sym = orbital[index+1:]
    self.j = eval(j_sym)
    self.l = op.ConvertLSym(self.l_sym)

    self.j = approximate_to_half_integer(self.j)
    self.s = approximate_to_half_integer(self.s)
    self.spin_a = approximate_to_half_integer(self.spin_a)
    self.spin_b = approximate_to_half_integer(self.spin_b)

    passJ = False
    if obeys_triangle_rule(self.spin_A, self.spin_B, self.j):
      passJ = True
    else:
      print(f"the orbital spin-J ({self.j}) does not consver J({nu_A}) + J({nu_B}) = {self.spin_A} + {self.spin_B}.")

    passS = False
    if obeys_triangle_rule(self.spin_a, self.spin_b, self.s):
      passS = True
    else:
      print(f"the orbital spin-s ({self.s}) does not consver S({nu_a}) + S({nu_b}) = {self.spin_a} + {self.spin_b}.")

    passL = False
    if obeys_triangle_rule(self.j, self.s, self.l):
      passL = True
    else:
      print(f"the orbital spin-L ({self.l}) does not consver J({self.j}) + S({self.s}).")

    self.isSpinBalanced = passJ * passS * passL
    if self.isSpinBalanced == False :
      print("Fail angular momentum conservation.")
      return False
    # else:
      # print("All Spin are balance.")

    self.reactionStr = f"{nu_A}({spin_A_str})({nu_a},{nu_b}){nu_B}({ExB:.3f}|{JB}, {orbital}) @ {ELabPerU:.1f} MeV/u"

    self.Q_value = mass_A + mass_a - mass_b - mass_B - ExB
    self.dwI = DistortedWave(nu_A, nu_a, self.ELab)

    self.mass_I = self.dwI.mu
    self.k_I = self.dwI.k

    Ecm_I = self.dwI.Ecm
    Ecm_O = Ecm_I + self.Q_value
    self.Eout = ((Ecm_O + mass_b + mass_B + ExB)**2 - (mass_b + mass_B + ExB)**2)/2/(mass_B + ExB)
    self.dwO = DistortedWave(nu_B, nu_b, self.Eout)

    self.mass_O = self.dwO.mu

    # Eout2 = self.ELab + self.Q_value #this is incorrec, but used in ptolmey infileCreator

    print("==================================================")
    print(self.reactionStr)
    print(f"Transfer Orbtial : {orbital}")
    print(f"Q-value : {self.Q_value:10.6f} MeV")
    print(f"Binding : {self.BindingEnergy:10.6f} MeV")
    # print(f"   Eout : {self.Eout} MeV | {Eout2}")
    print(f"   Eout : {self.Eout} MeV ")

    return True
