#!/usr/bin/env python3
import sys, os
import re
import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import time
from sympy import S
from sympy.physics.quantum.cg import wigner_9j

sys.path.append(os.path.join(os.path.dirname(__file__), '../Cleopatra'))
from IAEANuclearData import IsotopeClass

from assLegendreP import associated_legendre_array
from clebschGordan import clebsch_gordan, quantum_factorial, obeys_triangle_rule
from boundState import BoundState
from solveSE import WoodsSaxonPot, CoulombPotential, SpinOrbit_Pot, WS_SurfacePot
from distortedWave import DistortedWave
import opticalPotentials as op

class DWBA_ZR:
  def __init__(self, nu_A:str, nu_a:str, nu_b:str, nu_B:str, JB:str, orbital:str, ExB:float, ELabPerU:float):
    start_time = time.time()
    iso = IsotopeClass()
    
    A_A, Z_A = iso.GetAZ(nu_A)
    A_a, Z_a = iso.GetAZ(nu_a)
    A_b, Z_b = iso.GetAZ(nu_b)
    A_B, Z_B = iso.GetAZ(nu_B)

    self.ELab = A_a * ELabPerU

    mass_A = iso.GetMassFromSym(nu_A)
    mass_a = iso.GetMassFromSym(nu_a)
    mass_b = iso.GetMassFromSym(nu_b)
    mass_B = iso.GetMassFromSym(nu_B)
    self.ExB = ExB

    # sym_A = iso.GetSymbol(A_A, Z_A)
    # sym_B = iso.GetSymbol(A_B, Z_B)

    spin_A_str = iso.GetJpi(A_A, Z_A)
    self.spin_A = float(eval(re.sub(r'[+-]', '', spin_A_str)))
    self.spin_B = float(eval(re.sub(r'[+-]', '', JB)))

    if A_a == 2 and Z_a == 1:
        self.spin_a = 1.0
        self.spin_b = 0.5
    else:
        self.spin_a = 0.5
        self.spin_b = 1.0

    #====== transfering nucleon
    self.s = 1/2 # spin of x, neutron or proton
    A_x = abs(A_a - A_b)
    Z_x = abs(Z_a - Z_b)

    mass_x = iso.GetMassFromAZ( A_x, Z_x)

    #======== core
    if A_A < A_B : # (d,p)
        A_c = A_A
        Z_c = Z_A
        BindingEnergy = mass_B - mass_A - mass_x + self.ExB
    else:  #(p,d)
        A_c = A_B
        Z_c = Z_B
        BindingEnergy = mass_A - mass_B - self.ExB - mass_x

    #=================== digest orbital
    match = re.search(r'[a-zA-Z]', orbital)  # Find first letter
    if match:
        index = match.start()  # Get position of the first letter
        
    node = int(orbital[:index])
    l_sym = orbital[index:index+1]
    j_sym = orbital[index+1:]
    self.j = eval(j_sym)
    self.l = op.ConvertLSym(l_sym)

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
      return
    else:
      print("All Spin are balance.")

    self.reactionStr = f"{nu_A}({spin_A_str})({nu_a},{nu_b}){nu_B}({ExB:.3f}|{JB}, {orbital}) @ {ELabPerU:.1f} MeV/u"
    print("==================================================")
    print(self.reactionStr)

    self.Q_value = mass_A + mass_a - mass_b - mass_B - ExB
    print(f"Transfer Orbtial : {orbital}")
    print(f"Q-value : {self.Q_value:10.6f} MeV")
    print(f"Binding : {BindingEnergy:10.6f} MeV")

    print("====================== Bound state ")
    self.boundState = BoundState(A_c, Z_c, A_x, Z_x, node, self.l, self.j, BindingEnergy)
    self.boundState.SetPotential(1.10, 0.65, -6, 1.25, 0.65, 1.30)
    self.boundState.FindPotentialDepth(-70, -45, 0.5)

    print("====================== Incoming wave function ")
    op.AnCai(A_A, Z_A, self.ELab)

    self.maxL = 15

    self.dwI = DistortedWave(nu_A, nu_a, self.ELab)
    self.dwI.maxL = self.maxL
    self.dwI.PrintInput()
    self.dwI.ClearPotential()
    self.dwI.AddPotential(WoodsSaxonPot(   -op.v,    op.r0,    op.a), False)
    self.dwI.AddPotential(WoodsSaxonPot(-1j*op.vi,   op.ri0,   op.ai), False)
    self.dwI.AddPotential(WS_SurfacePot(-1j*op.vsi,  op.rsi0,  op.asi), False)
    self.dwI.AddPotential(SpinOrbit_Pot(   -op.vso,  op.rso0,  op.aso), False)
    self.dwI.AddPotential(SpinOrbit_Pot(-1j*op.vsoi, op.rsoi0, op.asoi), False)
    self.dwI.AddPotential(CoulombPotential( op.rc0), False)
    self.dwI.PrintPotentials()

    self.mass_I = self.dwI.mu # reduced mass of incoming channel
    self.k_I = self.dwI.k # wave number of incoming channel
    self.maxL = self.dwI.maxL

    sm_I, wfu_I = self.dwI.CalScatteringMatrix()

    self.dwI.PrintScatteringMatrix()

    Ecm_I = self.dwI.Ecm
    Ecm_O = Ecm_I + self.Q_value
    Eout = ((Ecm_O + mass_b + mass_B + self.ExB)**2 - (mass_b + mass_B + ExB)**2)/2/mass_B

    print("====================== Outgoing wave function ")
    op.Koning(A_B, Z_B, self.ELab + self.Q_value - ExB, Z_b)

    self.dwO = DistortedWave(nu_B, nu_b, Eout)
    self.dwO.spin_A = self.spin_B
    self.dwO.maxL = self.maxL
    self.dwO.PrintInput()
    self.dwO.ClearPotential()
    self.dwO.AddPotential(WoodsSaxonPot(   -op.v,    op.r0,    op.a), False)
    self.dwO.AddPotential(WoodsSaxonPot(-1j*op.vi,   op.ri0,   op.ai), False)
    self.dwO.AddPotential(WS_SurfacePot(-1j*op.vsi,  op.rsi0,  op.asi), False)
    self.dwO.AddPotential(SpinOrbit_Pot(   -op.vso,  op.rso0,  op.aso), False)
    self.dwO.AddPotential(SpinOrbit_Pot(-1j*op.vsoi, op.rsoi0, op.asoi), False)
    self.dwO.AddPotential(CoulombPotential( op.rc0), False)

    self.dwO.PrintPotentials()

    sm_O, wfu_O_temp = self.dwO.CalScatteringMatrix()

    #============ rescale the outgoing wave 
    rpos_O_temp = self.dwO.rpos * A_B/A_A
    self.rpos_O = []
    rpos_O_filled = False
    self.wfu_O = []
    for L in range(0, self.maxL+1):
      temp_wfu_L = []
      for k in range(0, int(2*self.spin_b)+1):
        wfu_O_inter_real = interp1d(rpos_O_temp, np.real(wfu_O_temp[L][k]), kind='cubic')
        wfu_O_inter_imag = interp1d(rpos_O_temp, np.imag(wfu_O_temp[L][k]), kind='cubic')
        temp_wfu = []
        for r in self.dwI.rpos:
          if r > 20 :
             break
          if rpos_O_filled == False:
            self.rpos_O.append(r)
          temp_wfu.append(wfu_O_inter_real(r) + 1j * wfu_O_inter_imag(r))
        rpos_O_filled = True
        temp_wfu_L.append(temp_wfu)
      self.wfu_O.append(temp_wfu_L)
  
    self.dwO.PrintScatteringMatrix()

    print("====================== Calculating Radial integrals")

    self.radialInt = np.zeros((self.maxL+1, int(2*self.spin_a+1), int(2*self.l+1), int(2*self.spin_b+1)), dtype=complex)

    bs = self.boundState.GetBoundStateWF()

    for L1 in range(0, self.maxL+1):
      for index1 in range(0, len(wfu_I[L1])):
        wf1 = wfu_I[L1][index1]
        for L2 in np.arange(abs(L1 - self.l), L1 + self.l + 1, 1):
          for index2 in range(0, len(self.wfu_O[int(L2)])):
            wf2 = self.wfu_O[int(L2)][index2]
            pf1 = np.exp(1j*self.dwI.CoulombPhaseShift(L1))
            pf2 = np.exp(1j*self.dwO.CoulombPhaseShift(L2))
            integral = simpson (bs[:200]*wf1[:200]*wf2[:200], dx=self.boundState.dr)
            indexL2 = int(L2 - abs(L1-self.l))
            self.radialInt[L1][index1][indexL2][index2] = integral * pf1 * pf2

    mass_I = self.dwI.mu
    k_I = self.dwI.k
    mass_O = self.dwO.mu # reduced mass of outgoing channel
    k_O = self.dwO.k # wave number of outgoing channel

    D0 = 1.55e+4 # for (d,p)

    self.massFactor = A_B/A_A
    self.ffactor = np.sqrt(4*np.pi)/k_I /k_O * A_B/A_A

    self.xsecScalingfactor = A_B / A_A *D0 * mass_I * mass_O / np.pi / self.dwI.hbarc**4 / k_I**3 / k_O * (2*self.spin_B + 1) / (2*self.spin_A+1) / (2*self.spin_a +1)

    stop_time = time.time()
    print(f"Total time {(stop_time - start_time) * 1000:.2f} msec")

  #========== end of contructor

  def FormatSpin(self, spin : float) -> str:
    if int(2*spin) % 2 == 0 :
      return f"{int(spin):+d}"
    else:
      return f"{int(2*spin):+d}/2"


  def PrintRadialIntegral(self):   
    for index1 in range(0, int(2*self.spin_a) + 1):
      for index2 in range(0, int(2*self.spin_b) + 1):
        print(f"======================= J1 = L{self.FormatSpin(index1-self.spin_a)}, J2 = L{self.FormatSpin(index2-self.spin_b)}")
        for L1 in range(0, self.maxL+1):
          print("{", end="")
          for L2 in np.arange(abs(L1 - self.l), L1 + self.l + 1):
            J1 = L1 + index1 - self.spin_a
            J2 = int(L2) + index2 - self.spin_b
            indexL2 = int(L2 - abs(L1-self.l))
            print(f"{{{L1:2d}, {J1:4.1f}, {int(L2):2d}, {J2:4.1f}, {np.real(self.radialInt[L1][index1][indexL2][index2]):12.4e} +  {np.imag(self.radialInt[L1][index1][indexL2][index2]):12.4e}I}},", end="")
          print("},")
    print("=========================== end of Radial Integrals.")

  def PlotRadialIntegral(self):
    spin_b = self.spin_b
    spin_a = self.spin_a
    l = int(self.l)
    maxL = self.maxL

    fig, axes = plt.subplots(int(2*spin_b+1)*int(2*l+1), int(2*spin_a+1), figsize=(6*int(2*spin_a+1), 4*int(2*spin_b+1)*int(2*l+1)))

    for index2 in range(0, int(2*spin_b) + 1):
      for index1 in range(0, int(2*spin_a) + 1):
        l_list = []
        for indexL2 in range(0, int(2*l) + 1):
          haha = []
          for L1 in range(0, maxL+1):
              J1 = L1 + index1 - spin_a
              L2 = int(abs(L1-l) + indexL2)
              J2 = L2 + index2 - spin_b
              if J1 < 0 or J2 < 0 :
                  continue
              l_list.append(L1)
              haha.append(self.radialInt[L1][index1][indexL2][index2])
          axes[int(2*l+1)*index2 + indexL2, index1].plot(l_list, np.real(haha), label="Real", marker='o')
          axes[int(2*l+1)*index2 + indexL2, index1].plot(l_list, np.imag(haha), label="Imag", marker='x')
          axes[int(2*l+1)*index2 + indexL2, index1].legend()
          axes[int(2*l+1)*index2 + indexL2, index1].set_xlabel('L1')
          axes[int(2*l+1)*index2 + indexL2, index1].set_ylabel('Value')
          axes[int(2*l+1)*index2 + indexL2, index1].set_title(f'Radial Int. vs L for Spin J1 = L{self.FormatSpin(index1-spin_a)}, L2 = L1{indexL2-l:+d}, J2 = L{self.FormatSpin(index2-spin_b)}.')
          axes[int(2*l+1)*index2 + indexL2, index1].set_xlim(-1, maxL+1)
          axes[int(2*l+1)*index2 + indexL2, index1].grid()

    plt.tight_layout()
    plt.show(block=False)
    input("Press Enter to continue...")

  def PlotScatteringMatrix(self, isIncoming):
    if isIncoming :
      self.dwI.PlotScatteringMatrix()
    else:
      self.dwO.PlotScatteringMatrix()

  def PlotDistortedWave(self, isIncoming, L, J, maxR = None):
    if isIncoming:
      self.dwI.PlotDistortedWave(L, J, maxR)
    else:
      plt.plot(self.rpos_O, np.real(self.wfu_O[L][int(J-L + self.spin_b)]), label="Real")
      plt.plot(self.rpos_O, np.imag(self.wfu_O[L][int(J-L + self.spin_b)]), label="Imag")
      plt.title(f"Radial wave function for L={L} and J={J}")
      if maxR != None :
        plt.xlim(-1, maxR) 
      plt.legend()
      plt.grid()
      plt.show(block=False)
      input("Press Enter to continue...")

  def Gamma(self, L1:int, J1, L2:int, J2, m:int, ma, mb):
    if  int(L1 + L2 + self.l)%2 != 0: #check if the sum of L1 + L2 + l is even
        return 0
    else:
        fact0 = wigner_9j(S(2*self.j)/2, self.l, S(2*self.s)/2, S(2*J1)/2, L1, S(2*self.spin_a)/2, J2, S(2*L2)/2, S(2*self.spin_b)/2).evalf()
        if fact0 == 0:
            return 0
        else:
            fact1 = pow(-1, m) * np.power(1j, L1-L2-self.l) * (2*L2+1) * np.sqrt((2*self.l+1)*(2*self.s+1)*(2*L1+1)*(2*J2+1))
            fact2 = np.sqrt( quantum_factorial(L2-m) / quantum_factorial(L2+m) )
            fact3 = clebsch_gordan(J2, mb-m,      self.j, m-mb+ma, J1,   ma)
            fact4 = clebsch_gordan(L1,    0, self.spin_a,      ma, J1,   ma)
            fact5 = clebsch_gordan(L2,   -m, self.spin_b,      mb, J2, mb-m)
            fact6 = clebsch_gordan(L1,    0, int(self.l),       0, L2,    0)
            return fact0 * fact1 * fact2 * fact3 * fact4 * fact5 * fact6

  def Beta(self, m:int, ma, mb):
    result = 0
    for L1 in np.arange(0, self.maxL+1):
      for J1 in np.arange(abs(L1 - self.spin_a), L1 + self.spin_a + 1, 1):
        for L2 in np.arange(abs(L1 - self.l), L1 + self.l + 1, 1):
          for J2 in np.arange(abs(L2 - self.spin_b), L2 + self.spin_b + 1, 1):

            if not obeys_triangle_rule(J1, self.j, J2):
                continue
            if not(abs(m) <= L2):
                continue
            if int(L1 + L2 + self.l) % 2 != 0:
                continue
            
            index1 = int(J1 - L1 + self.spin_a)
            index2 = int(J2 - L2 + self.spin_b)
            indexL2 = int(L1 - abs(L1 - self.l))

            gg = self.Gamma(L1, J1, L2, J2, m, ma, mb)
            if gg == 0:
               continue
            lp = self.legendrePArray[int(L2)][int(abs(m))] 
            if m < 0 :
              lp *= (-1)**m * quantum_factorial(int(L2)+m)/ quantum_factorial(int(L2)-m)
            ri = self.radialInt[int(L1)][index1][indexL2][index2]
            # print(f"{L1:2d}, {J1:4.1f}({index1:d}), {L2:2d}({indexL2:d}), {J2:4.1f}({index2:d}), {gg:10.6f}, {ri *self.ffactor :.10f}, {lp:10.6f}")

            result +=  gg * lp * ri 

    return result  
  
  def PreCalLegendreP(self,  theta_deg:float, maxL:int = None, maxM:int = None):
    if maxL is None:
       maxL = self.maxL
    if maxM is None:
      maxM = int(self.j + self.spin_b + self.spin_a)
    self.legendrePArray = associated_legendre_array(maxL, maxM, theta_deg)
     
  def AngDist(self, theta_deg):
    xsec = 0
    self.PreCalLegendreP(theta_deg)
    for ma in np.arange(-self.spin_a, self.spin_a + 1, 1):
      for mb in np.arange(-self.spin_b, self.spin_b + 1, 1):
        for m in np.arange(-self.j + mb - ma, self.j + mb -ma + 1, 1):
          haha = self.Beta(m, ma, mb)
          xsec += np.abs(haha)**2
    
    return xsec * self.xsecScalingfactor * 10 # factor 10 for fm^2 = 10 mb







