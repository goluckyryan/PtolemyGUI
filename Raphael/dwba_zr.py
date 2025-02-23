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

    self.j = self.approximate_to_half_integer(self.j)
    self.s = self.approximate_to_half_integer(self.s)
    self.spin_a = self.approximate_to_half_integer(self.spin_a)
    self.spin_b = self.approximate_to_half_integer(self.spin_b)

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
    

    print("====================== Incoming wave function ")
    op.AnCai(A_A, Z_A, self.ELab)


    self.dwI = DistortedWave(nu_A, nu_a, self.ELab)
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
    self.maxL1 = self.dwI.maxL

    self.maxL2 = self.maxL1 + self.l

    
    Ecm_I = self.dwI.Ecm
    Ecm_O = Ecm_I + self.Q_value
    Eout = ((Ecm_O + mass_b + mass_B + self.ExB)**2 - (mass_b + mass_B + ExB)**2)/2/mass_B

    print("====================== Outgoing wave function ")
    op.Koning(A_B, Z_B, self.ELab + self.Q_value - ExB, Z_b)

    self.dwO = DistortedWave(nu_B, nu_b, Eout)
    self.dwO.spin_A = self.spin_B
    self.dwO.maxL = self.maxL2
    self.dwO.PrintInput()
    self.dwO.ClearPotential()
    self.dwO.AddPotential(WoodsSaxonPot(   -op.v,    op.r0,    op.a), False)
    self.dwO.AddPotential(WoodsSaxonPot(-1j*op.vi,   op.ri0,   op.ai), False)
    self.dwO.AddPotential(WS_SurfacePot(-1j*op.vsi,  op.rsi0,  op.asi), False)
    self.dwO.AddPotential(SpinOrbit_Pot(   -op.vso,  op.rso0,  op.aso), False)
    self.dwO.AddPotential(SpinOrbit_Pot(-1j*op.vsoi, op.rsoi0, op.asoi), False)
    self.dwO.AddPotential(CoulombPotential( op.rc0), False)

    self.dwO.PrintPotentials()

    self.radialInt = None

    mass_I = self.dwI.mu
    k_I = self.dwI.k
    mass_O = self.dwO.mu # reduced mass of outgoing channel
    k_O = self.dwO.k # wave number of outgoing channel

    D0 = 1.55e+4 # for (d,p)

    self.massBoverMassA = A_B/A_A
    self.ffactor = np.sqrt(4*np.pi)/k_I /k_O 

    self.xsecScalingfactor = D0 * mass_I * mass_O / np.pi / self.dwI.hbarc**4 / k_I**3 / k_O * (2*self.spin_B + 1) / (2*self.spin_A+1) / (2*self.spin_a +1)

    self.PreCalNineJ()

  #========== end of contructor

  def FormatSpin(self, spin : float) -> str:
    if int(2*spin) % 2 == 0 :
      return f"{int(spin):+d}"
    else:
      return f"{int(2*spin):+d}/2"

  def FindBoundState(self):
    self.boundState.FindPotentialDepth(-70, -45, 0.5)

  def ConvertLJ2RadialIndex(self, L1:int, J1, L2:int, J2):
    index1 = int(J1 - L1 + self.spin_a)
    indexL2 = int(L2 - L1 + self.l)
    index2 = int(J2 - L2 + self.spin_b)
    return [L1, index1, indexL2, index2]

  def ConvertRadialIndex2LJ(self, in1:int, in2:int, in3:int, in4:int):
    L1 = in1
    J1 = L1 + in2 - self.spin_a
    L2 = in3 + L1 - self.l
    J2 = L2 + in4 - self.spin_b
    return [L1, J1, L2, J2] 

###########################################################
  def CalRadialIntegral(self):
    start_time = time.time()
    sm_I, wfu_I = self.dwI.CalScatteringMatrix()
    self.dwI.PrintScatteringMatrix()
    
    sm_O, wfu_O_temp = self.dwO.CalScatteringMatrix()
    self.dwO.PrintScatteringMatrix()

    #============ rescale the outgoing wave 
    print("====================== Scaling the outgoing wave")
    rpos_O_temp = self.dwO.rpos * self.massBoverMassA
    self.rpos_O = []
    rpos_O_filled = False
    self.wfu_O = []
    for L2 in range(0, self.maxL2 + 1):
      temp_wfu_L = []
      for k in range(0, int(2*self.spin_b)+1):
        wfu_O_inter_real = interp1d(rpos_O_temp, np.real(wfu_O_temp[L2][k]), kind='cubic')
        wfu_O_inter_imag = interp1d(rpos_O_temp, np.imag(wfu_O_temp[L2][k]), kind='cubic')
        temp_wfu = []
        for r in self.dwI.rpos:
          if r > max(rpos_O_temp) :
             break
          if rpos_O_filled == False:
            self.rpos_O.append(r)
          temp_wfu.append(wfu_O_inter_real(r) + 1j * wfu_O_inter_imag(r))
        rpos_O_filled = True
        temp_wfu_L.append(temp_wfu)
      self.wfu_O.append(temp_wfu_L)

    print("====================== Calculating Radial integrals")
    self.radialInt = np.zeros((self.maxL1+1, int(2*self.spin_a+1), int(2*self.l+1), int(2*self.spin_b+1)), dtype=complex)
    bs = self.boundState.GetBoundStateWF()

    for L1 in range(0, self.maxL1+1):
      for J1 in np.arange(L1-self.spin_a, L1 + self.spin_a + 1, 1):
        if J1 < 0 :
          continue
        index1 = int(J1 - L1 + self.spin_a)
        wf1 = wfu_I[L1][index1]
        for L2 in np.arange(L1 - self.l, L1 + self.l + 1, 1):
          if L2 < 0 :
            continue
          for J2 in np.arange(L2 - self.spin_b, L2 + self.spin_b + 1, 1):
            if J2 < 0:
              continue
            index2 = int(J2 - L2 + self.spin_b)
            wf2 = self.wfu_O[int(L2)][index2]
            pf1 = np.exp(1j*self.dwI.CoulombPhaseShift(L1))
            pf2 = np.exp(1j*self.dwO.CoulombPhaseShift(L2))
            integral = simpson (bs*wf1*wf2, dx=self.boundState.dr)
            indexL2 = int(L2 - L1 + self.l)
            product = integral * pf1 * pf2 * self.massBoverMassA
            self.radialInt[L1][index1][indexL2][index2] = product
            # if J1 == L1 + self.spin_a and L2 == L1 + 1  and J2 == L2 - self.spin_b:
              # print(f"{L1:2d}, {J1:4.1f}({index1}), {L2:2d}({indexL2}), {J2:4.1f}({index2}), {integral * pf1 * pf2 * self.massBoverMassA:.6f}")

    stop_time = time.time()
    print(f"Total time for distorted wave and radial intergal {(stop_time - start_time) * 1000:.2f} msec")
###########################################################

  def PrintRadialIntegral(self):
    for index1 in range(0, int(2*self.spin_a) + 1):
      for index2 in range(0, int(2*self.spin_b) + 1):
        print(f"======================= J1 = L{self.FormatSpin(index1-self.spin_a)}, J2 = L{self.FormatSpin(index2-self.spin_b)}")
        for L1 in range(0, self.maxL1+1):
          print("{", end="")
          for L2 in np.arange(L1 - self.l, L1 + self.l + 1, 1):
            J1 = L1 + index1 - self.spin_a
            J2 = int(L2) + index2 - self.spin_b
            indexL2 = int(L2 - L1  + self.l)
            print(f"{{{L1:2d}, {J1:4.1f}, {int(L2):2d}, {J2:4.1f}, {np.real(self.radialInt[L1][index1][indexL2][index2]):11.4e}+{np.imag(self.radialInt[L1][index1][indexL2][index2]):11.4e}I}},  ", end="")
          print("},")
    print("=========================== end of Radial Integrals.")

  def PlotRadialIntegral(self):
    if self.radialInt is None:
      print("Radial integral not computed.")
      return
    spin_b = self.spin_b
    spin_a = self.spin_a
    l = self.l
    maxL1 = self.maxL1

    fig, axes = plt.subplots(int(2*spin_b+1)*int(2*l+1), int(2*spin_a+1), figsize=(6*int(2*spin_a+1), 4*int(2*spin_b+1)*int(2*l+1)))

    for index2 in range(0, int(2*spin_b) + 1):
      for index1 in range(0, int(2*spin_a) + 1):
        for indexL2 in range(0, int(2*l) + 1):
          haha = []
          l_list = []
          for L1 in range(0, maxL1+1):
              # J1 = L1 + index1 - spin_a
              # L2 = int(L1 - l + indexL2)
              # J2 = L2 + index2 - spin_b
              # [dummy, J1, L2, J2] = self.ConvertRadialIndex2LJ(L1, index2, indexL2, index2)
              l_list.append(L1)
              haha.append(self.radialInt[L1][index1][indexL2][index2])
          axes[int(2*l+1)*index2 + indexL2, index1].plot(l_list, np.real(haha), label="Real", marker='o')
          axes[int(2*l+1)*index2 + indexL2, index1].plot(l_list, np.imag(haha), label="Imag", marker='x')
          axes[int(2*l+1)*index2 + indexL2, index1].legend()
          axes[int(2*l+1)*index2 + indexL2, index1].set_xlabel('L1')
          axes[int(2*l+1)*index2 + indexL2, index1].set_ylabel('Value')
          axes[int(2*l+1)*index2 + indexL2, index1].set_title(f'Radial Int. vs L for Spin J1 = L1{self.FormatSpin(index1-spin_a)}, L2 = L1{indexL2-l:+d}, J2 = L2{self.FormatSpin(index2-spin_b)}.')
          axes[int(2*l+1)*index2 + indexL2, index1].set_xlim(-1, maxL1+1)
          axes[int(2*l+1)*index2 + indexL2, index1].set_xticks(np.arange(0, maxL1+1, 2))
          axes[int(2*l+1)*index2 + indexL2, index1].grid()

    plt.tight_layout()
    plt.show(block=False)
    input("Press Enter to continue...")

  def PlotRadialIntegralSigle(self, dJ1, dL, dJ2):
    if self.radialInt is None:
      print("Radial integral not computed.")
      return
    haha = []
    l_list = []
    for L1 in range(0, self.maxL1+1):
      l_list.append(L1)
      [dummy, index1, indexL2, index2] = self.ConvertLJ2RadialIndex(L1, L1 + dJ1, L1 + dL, L1 + dL + dJ2)
      haha.append(self.radialInt[L1][index1][indexL2][index2])
      print(f"{L1:2d}, {L1 + dJ1:4.1f}({index1}), {L1 + dL:.0f}({indexL2}), {L1 + dL + dJ2:4.1f}({index2}), {haha[-1]:.6f}")
    plt.plot(l_list, np.real(haha), label="Real", marker='o')
    plt.plot(l_list, np.imag(haha), label="Imag", marker='x')
    plt.xlabel("L1")
    plt.title(f'Radial Int. vs L for Spin J1 = L1{self.FormatSpin(dJ1)}, L2 = L1{dL:+d}, J2 = L2{self.FormatSpin(dJ2)}.')
    plt.grid()
    plt.show(block=False)
    input("Press Enter to continue...")

  def PlotScatteringMatrix(self, isIncoming):
    if isIncoming :
      self.dwI.PlotScatteringMatrix()
    else:
      self.dwO.PlotScatteringMatrix()

  def PlotIncomingDistortedWave(self,  L, J, maxR = None):
    self.dwI.PlotDistortedWave(L, J, maxR)
    
  
  def PlotOutgoingDistortedWave(self, L, J, maxR = None):
    plt.plot(self.rpos_O, np.real(self.wfu_O[L][int(J-L + self.spin_b)]), label="Real")
    plt.plot(self.rpos_O, np.imag(self.wfu_O[L][int(J-L + self.spin_b)]), label="Imag")
    plt.title(f"Radial wave function for L={L} and J={J}")
    if maxR != None :
      plt.xlim(-1, maxR) 
    plt.legend()
    plt.grid()
    plt.show(block=False)
    input("Press Enter to continue...")
  
  def approximate_to_half_integer(self, value):
    return round(value * 2) / 2
  
  def PreCalNineJ(self):
    self.NineJ = np.zeros((self.maxL1+1, int(2*self.spin_a+1), (2*self.l+1), int(2*self.spin_b+1)), dtype=float)
    for L1 in range(0, self.maxL1+1):
      for ind1 in range(0,  int(2*self.spin_a+1)):
        for indL2 in range(0, 2*self.l+1):
          for ind2 in range(0, int(2*self.spin_b+1)):
            J1 = self.approximate_to_half_integer(L1 + ind1 - self.spin_a)
            L2 = int(L1 + indL2 - self.l)
            J2 = self.approximate_to_half_integer(L2 + ind2 - self.spin_b)
            self.NineJ[L1, ind1, indL2, ind2] = wigner_9j(self.j, self.l, self.s, J1, L1, self.spin_a, J2, L2, self.spin_b)

  def GetPreCalNineJ(self, L1:int, J1, L2:int, J2):
    [dummy, ind1, indL2, ind2] = self.ConvertLJ2RadialIndex(L1, J1, L2, J2)
    return self.NineJ[L1, ind1, indL2, ind2]

  def Gamma(self, L1:int, J1, L2:int, J2, m:int, ma, mb):
    if  int(L1 + L2 + self.l)%2 != 0: #check if the sum of L1 + L2 + l is even
        return 0
    else:
        fact0 = self.GetPreCalNineJ(L1, J1, L2, J2)
        if fact0 == 0:
            return 0
        else:
            fact1 = pow(-1, m) * np.power(1j, L1-L2-self.l) * (2*L2+1) * np.sqrt((2*self.l+1)*(2*self.s+1)*(2*L1+1)*(2*J2+1))
            fact2 = np.sqrt( quantum_factorial(L2-abs(m)) / quantum_factorial(L2 + abs(m)) )
            fact3 = clebsch_gordan(J2, mb-m,      self.j, m-mb+ma, J1,   ma)
            fact4 = clebsch_gordan(L1,    0, self.spin_a,      ma, J1,   ma)
            fact5 = clebsch_gordan(L2,   -m, self.spin_b,      mb, J2, mb-m)
            fact6 = clebsch_gordan(L2,    0,      self.l,       0, L1,    0)
            # print(f"{fact1:.5f}, {fact2:.5f}, {fact3:.5f}, {fact4:.5f}, {fact5:.5f}, {fact6:.5f}")
            return fact0 * fact1 * fact2 * fact3 * fact4 * fact5 * fact6

  def Beta(self, m:int, ma, mb):
    if self.radialInt is None :
      return 
    result = 0
    for L1 in np.arange(0, self.maxL1+1):
      for J1 in np.arange(L1 - self.spin_a, L1 + self.spin_a + 1, 1):
        if J1 < 0 :
          continue
        for L2 in np.arange(L1 - self.l, L1 + self.l + 1, 1):
          if L2 < 0:
            continue
          for J2 in np.arange(L2 - self.spin_b, L2 + self.spin_b + 1, 1):
            if J2 < 0 :
              continue
            if not obeys_triangle_rule(J1, self.j, J2):
                continue
            if not(abs(m) <= L2):
                continue
            if int(L1 + L2 + self.l) % 2 != 0:
                continue
            

            gg = self.Gamma(L1, J1, L2, J2, m, ma, mb)
            if gg == 0:
               continue
            lp = self.legendrePArray[L2][int(abs(m))] 

            [dummy, index1, indexL2, index2] = self.ConvertLJ2RadialIndex(L1, J1, L2, J2)
            ri = self.radialInt[int(L1)][index1][indexL2][index2]
            # print(f"{L1:2d}, {J1:4.1f}({index1:d}), {L2:2d}({indexL2:d}), {J2:4.1f}({index2:d}), {gg:10.6f}, {ri *self.ffactor :.10f}, {lp:10.6f}")

            result +=  gg * lp * ri 

    return result  
  
  def PreCalLegendreP(self,  theta_deg:float, maxL:int = None, maxM:int = None):
    if maxL is None:
       maxL = max(self.maxL1, self.maxL2)
    if maxM is None:
      maxM = int(self.j + self.spin_b + self.spin_a)
    self.legendrePArray = associated_legendre_array(maxL, maxM, theta_deg)
     
  def AngDist(self, theta_deg:float) -> float:
    xsec = 0
    self.PreCalLegendreP(theta_deg)
    for ma in np.arange(-self.spin_a, self.spin_a + 1, 1):
      for mb in np.arange(-self.spin_b, self.spin_b + 1, 1):
        for m in np.arange(-self.j + mb - ma, self.j + mb -ma + 1, 1):
          haha = self.Beta(m, ma, mb)
          xsec += np.abs(haha)**2
    
    return xsec * self.xsecScalingfactor * 10 # factor 10 for fm^2 = 10 mb

  def CalAngDistribution(self, angMin:float = 0, angMax:float = 180, angStep:float = 1):
    self.angMin = angMin
    self.angMax = angMax
    self.angList = []
    self.angDist = []
    print(f"======== Calcalating Angular distribution from {angMin:.1f} to {angMax:.1f}, step {angStep:.1f}...")
    start_time = time.time()
    progress_time = time.time()
    for i in np.arange(angMin, angMax + angStep, angStep):
      self.angList.append(i)
      self.angDist.append(self.AngDist(i))
      if time.time() - progress_time > 1:
        elapsed_time = time.time() - start_time
        print(f"\r Time elapsed: {elapsed_time:.2f} sec, Progress: {100 * (i - angMin) / (angMax - angMin):.1f}%", end="")
        progress_time = time.time()
    stop_time = time.time()
    print(f"\nTotal time {(stop_time - start_time) :.2f} sec")

  def PrintAngDist(self):
    for th, xs in zip(self.angList, self.angDist):
      print(f"{th:6.1f}, {xs:13.10f}")

  def PlotAngDist(self, angMin = None, angMax = None):
    plt.plot(self.angList, self.angDist)
    plt.title(self.reactionStr)
    if angMin is None and angMax is None:
      plt.xlim(-1 + self.angMin, self.angMax + 1)
    if angMin is None and angMax != None:
      plt.xlim(-1 + self.angMin, angMax + 1)
    if angMin != None and angMax is None :
      plt.xlim(-1 + angMin, self.angMax + 1)

    plt.yscale("log")
    plt.grid()
    plt.show(block=False)
    input("Press Enter to continue...")



