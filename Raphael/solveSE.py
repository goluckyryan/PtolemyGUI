#!/usr/bin/env python3

import math
import numpy as np
import re
import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), '../Cleopatra'))
from IAEANuclearData import IsotopeClass

class PotentialForm:
  def __init__(self):
    self.V0 = -10
    self.r0 = 1.3
    self.a0 = 0.65
    self.id = -1

  def setCharge(self, Z):
    self.Charge = Z

  # def setA(self, A):
  #   self.R0 = self.r0 * math.pow(A, 1/3)

  def setAa(self, A, a):
    self.R0= self.r0 * (math.pow(A, 1/3) + math.pow(a, 1/3))

  def output(self, x):
    return 0
  
  def printPot(self, msg:str):
    print(f"{msg:20s} : V0 {np.real(self.V0):7.3f} + {np.imag(self.V0):7.3f} I, R0 {self.R0:6.4f}({self.r0:6.4f}), a0 {self.a0:6.4f}, ")

class CoulombPotential(PotentialForm):
  def __init__(self, rc):
    self.V0 = 0
    self.r0 = rc
    self.a0 = 0
    self.id = 0
    self.ee = 1.43996 # MeV.fm

  def output(self, x):
    if self.Charge == 0 :
      return 0
    else:
      if x >self.R0:
          return (self.Charge * self.ee) / (x + 1e-20)  # Add a small value to avoid division by zero
      else:
          return (self.Charge * self.ee) / (2 * self.R0) * (3 - (x / self.R0)**2)
    
  def printPot(self):
    return super().printPot("Coulomb")

class WoodsSaxonPot(PotentialForm):
  def __init__(self, V0, r0, a0) :
    self.V0 = V0
    self.r0 = r0
    self.a0 = a0
    self.id = 1

  def output(self, x):
    if self.V0 == 0.0:
      return 0
    else:
      return self.V0/(1 + math.exp((x-self.R0)/self.a0))

  def printPot(self):
    return super().printPot("Woods-Saxon")

class SpinOrbit_Pot(PotentialForm):
  def __init__(self, VSO, rSO, aSO) :
    # the LS factor is put in the SolvingSE Class
    self.V0 = VSO
    self.r0 = rSO
    self.a0 = aSO
    self.id = 2

  def output(self, x):
    if self.V0 == 0.0 :
      return 0
    else:
      if x > 0 :
        return 4*(self.V0 * math.exp((x-self.R0)/self.a0))/(self.a0*math.pow(1+math.exp((x-self.R0)/self.a0),2))/x
      else :
        return 4*1e+19

  def printPot(self):
    return super().printPot("Spin-Orbit")

class WS_SurfacePot(PotentialForm):
  def __init__(self, V0, r0, a0):
    self.V0 = V0
    self.r0 = r0
    self.a0 = a0
    self.id = 3

  def output(self, x):
    if self.V0 == 0 :
      return 0
    else:
      exponent = (x - self.R0) / self.a0
      return 4* self.V0 * math.exp(exponent) / (1 + math.exp(exponent))**2

  def printPot(self):
    return super().printPot("Woods-Saxon Surface")

#========================================
class SolvingSE:
  #grid setting
  rStart = 0.0
  dr = 0.05
  nStep = 600*5
  rpos = np.arange(rStart, rStart+nStep*dr, dr)
  solU = [] # raidal wave function
  maxSolU = 0.0

  #constant
  mn = 939.56539 #MeV/c2
  amu = 931.494 #MeV/c2
  hbarc = 197.326979 #MeV.fm
  ee = 1.43996 # MeV.fm

  #RK4 constants
  parC = [0, 1./2, 1./2, 1.]
  parD = [1./6, 2./6, 2./6, 1./6]

  #inital condition
  solu0 = 0.0
  dsolu0 = 0.0001

  potential_List = []

  def PrintInput(self):
    print(f"     A : ({self.A_A:3d}, {self.Z_A:3d}), spin : {self.spin_A},")
    print(f"     a : ({self.A_a:3d}, {self.Z_a:3d}), spin : {self.spin_a},")
    print(f"  Elab : {self.Energy : 10.5f} MeV")
    print(f"    mu : {self.mu: 10.5f} MeV/c2")
    print(f"   Ecm : {self.Ecm: 10.5f} MeV")
    print(f"     k : {self.k: 10.5f} MeV/c")
    print(f"   eta : {self.eta: 10.5f}")
    print(f"     L : {self.L},  maxL : {self.maxL}")
    print(f"    dr : {self.dr} fm, nStep : {self.nStep}")
    print(f"rStart : {self.rStart} fm,  rMax : {self.nStep * self.dr} fm")

  def __init__(self, A_or_SymA = None, ZA_or_Syma = None \
               , a_or_ELabPerA = None, Za_or_none = None, ELabPerA_or_none = None):
    if Za_or_none is None :
      self.ConstructUsingSymbol(A_or_SymA, ZA_or_Syma, a_or_ELabPerA)
    else:
      self.ConstructUsingAZ(A_or_SymA, ZA_or_Syma, a_or_ELabPerA, Za_or_none, ELabPerA_or_none)

    haha = IsotopeClass()
    self.mass_A = haha.GetMassFromAZ(self.A_A, self.Z_A)
    self.mass_a = haha.GetMassFromAZ(self.A_a, self.Z_a)
    self.spin_A = float(eval(re.sub(r'[+-]', '', haha.GetJpi(self.A_A, self.Z_A))))
    self.spin_a = float(eval(re.sub(r'[+-]', '', haha.GetJpi(self.A_a, self.Z_a))))
    self.S = self.spin_a

    self.mu = (self.mass_A * self.mass_a)/(self.mass_A + self.mass_a)
    self.Ecm = 0.0

  def ConstructUsingAZ(self, A, ZA, a, Za, ELabPerA):
    print(f"ConstructUsingAZ : {A}, {ZA}, {a}, {Za}, {ELabPerA}")
    self.A_A = A
    self.A_a = a
    self.Z_A = ZA
    self.Z_a = Za
    self.Z = ZA * Za

    self.L = 0
    self.S = 0
    self.J = 0
    self.Energy = ELabPerA

  def ConstructUsingSymbol(self, Sym_A, Sym_a, ELabPerA):
    print(f"ConstructUsingSymbol : {Sym_A}, {Sym_a}, {ELabPerA}")
    self.L = 0
    self.S = 0
    self.J = 0

    haha = IsotopeClass()
    self.A_A, self.Z_A = haha.GetAZ(Sym_A)
    self.A_a, self.Z_a = haha.GetAZ(Sym_a)
    self.Z = self.Z_A * self.Z_a
    self.Energy = ELabPerA


  def CalCMConstants(self, useELabAsEcm = False):
    if useELabAsEcm:
      self.E_tot = self.Energy # total energy in CM
      self.Ecm = self.Energy # KE in cm
    else:
      self.E_tot = math.sqrt(math.pow((self.mass_a+self.mass_A),2) + 2 * self.mass_A * self.Energy)
      self.Ecm = self.E_tot - (self.mass_a + self.mass_A) 
    
    self.k = math.sqrt(self.mu * 2 * abs(self.Ecm)) / self.hbarc
    if self.Z == 0 :
      self.eta = 0
    else:
      self.eta = self.Z * self.ee * self.k /2 /self.Ecm
    self.maxL = int(self.k * (1.4 * (self.A_A**(1/3) + self.A_a**(1/3)) + 5))

  # def SetA_ExSpin(self, ExA, sA ):
  #   self.ExA = ExA
  #   self.spin_A = sA

  def SetLJ(self, L, J):
    self.L = L
    self.J = J

  def LS(self, L = None, J = None) :
    if L is None:
      L = self.L
    if J is None:
      J = self.J
    return (J*(J+1)-L*(L+1)-self.S*(self.S))/2.

  # set the range in fm
  def SetRange(self, rStart, dr, nStep):
    self.rStart = rStart
    self.dr  = dr
    self.nStep = nStep
    self.rpos = np.arange(self.rStart, self.rStart+self.nStep*dr, self.dr)
    self.solU = []
    self.maxSolU = 0.0

  def ClearPotential(self):
    self.potential_List = []

  def AddPotential(self, pot : PotentialForm, useBothMass : bool = False):
    if isinstance(pot, PotentialForm):
      if pot.id == 0:
        pot.setCharge(self.Z)
      if useBothMass:
        pot.setAa(self.A_A, self.A_a)
      else:
        pot.setAa(self.A_A, 0)
      self.potential_List.append(pot)

  def __PotentialValue(self, x):
    value = 0
    for pot in self.potential_List:
      if isinstance(pot, PotentialForm):
        if pot.id == 2 and self.L > 0:
          value = value + self.LS() * pot.output(x)
        else:
          value = value + pot.output(x)
    return value
  
  def PrintPotentials(self):
    for pot in self.potential_List:
      if isinstance(pot, PotentialForm):
        pot.printPot()

  def GetPotentialValue(self, x):
    return self.__PotentialValue(x)

  # The G-function, u''[r] = G[r, u[r], u'[r]]
  def __G(self, x, y, dy):
    #return  -2*x*dy -2*y  # solution of gaussian
    if x > 0 :
      return 2*self.mu/math.pow(self.hbarc,2)*(self.__PotentialValue(x) - self.Ecm)*y + self.L*(1+self.L)/x/x*y
    else:
      return 0

  # Using Rungu-Kutta 4th method to solve u''[r] = G[r, u[r], u'[r]]
  def SolveByRK4(self):
    #initial condition
    self.solU = [self.solu0]
    dSolU = [self.dsolu0]

    dyy = np.array([1., 0., 0., 0., 0.], dtype= complex)
    dzz = np.array([1., 0., 0., 0., 0.], dtype= complex)

    self.maxSolU = 0.0

    for i in range(self.nStep-1):
      r = self.rStart + self.dr * i
      y = self.solU[i]
      z = dSolU[i]

      for j in range(4):
        dyy[j + 1] = self.dr * (z + self.parC[j] * dzz[j])
        dzz[j + 1] = self.dr * self.__G(r + self.parC[j] * self.dr, y + self.parC[j] * dyy[j], z + self.parC[j] * dzz[j])

      dy = sum(self.parD[j] * dyy[j + 1] for j in range(4))
      dz = sum(self.parD[j] * dzz[j + 1] for j in range(4))

      self.solU.append(y + dy)
      dSolU.append(z + dz)

      if np.real(self.solU[-1]) > self.maxSolU:
        self.maxSolU = abs(self.solU[-1])

    return self.solU
  
  def NearestPosIndex(self, r):
    return min(len(self.rpos)-1, int((r - self.rStart) / self.dr))

