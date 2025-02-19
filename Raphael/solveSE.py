#!/usr/bin/env python3

import math
import numpy as np

import sys, os

sys.path.append(os.path.join(os.path.dirname(__file__), '../Cleopatra'))
from IAEANuclearData import IsotopeClass

class CoulombPotential:
  def __init__(self, rc):
    self.rc = rc
    self.id = 0
    self.ee = 1.43996 # MeV.fm

  def setA(self, A):
    self.Rc = self.rc * math.pow(A, 1/3)

  def setAa(self, A, a):
    self.Rc = self.rc * (math.pow(A, 1/3) + math.pow(a, 1/3))

  def setCharge(self, Z):
    self.Charge = Z

  def output(self, x):
    if x >self.Rc:
        return (self.Charge * self.ee) / (x + 1e-20)  # Add a small value to avoid division by zero
    else:
        return (self.Charge * self.ee) / (2 * self.Rc) * (3 - (x / self.Rc)**2)

class WSPotential:
  def __init__(self, V0, r0, a0) :
    self.V0 = V0
    self.r0 = r0
    self.a0 = a0
    self.id = 1

  def setA(self, A):
    self.R0 = self.r0 * math.pow(A, 1/3)

  def setAa(self, A, a):
    self.R0 = self.r0 * (math.pow(A, 1/3) + math.pow(a, 1/3))

  def output(self, x):
    return self.V0/(1 + math.exp((x-self.R0)/self.a0))

class SpinOrbitPotential:
  def __init__(self, VSO, rSO, aSO) :
    # the LS factor is put in the SolvingSE Class
    self.VSO = VSO
    self.rSO = rSO
    self.aSO = aSO
    self.id = 2

  def setA(self, A):
    self.RSO = self.rSO * math.pow(A, 1/3)

  def setAa(self, A, a):
    self.RSO = self.rSO * (math.pow(A, 1/3) + math.pow(a, 1/3))

  def output(self, x):
    if x > 0 :
      return 4*(self.VSO * math.exp((x-self.RSO)/self.aSO))/(self.aSO*math.pow(1+math.exp((x-self.RSO)/self.aSO),2))/x
    else :
      return 4*1e+19

class WSSurface:
  def __init__(self, V0, r0, a0):
    self.V0 = V0
    self.r0 = r0
    self.a0 = a0
    self.id = 3

  def setA(self, A):
    self.R0 = self.r0 * math.pow(A, 1/3)

  def setAa(self, A, a):
    self.R0 = self.r0 * (math.pow(A, 1/3) + math.pow(a, 1/3))

  def output(self, x):
    exponent = (x - self.R0) / self.a0
    return self.V0 * math.exp(exponent) / (1 + math.exp(exponent))**2

#========================================
class SolvingSE:
  #grid setting
  rStart = 0.0
  dr = 0.05
  nStep = 600*5
  rpos = np.arange(rStart, rStart+nStep*dr, dr)
  SolU = [] # raidal wave function
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
    print(f"     A : ({self.A_A:3d}, {self.Z_A:3d})")
    print(f"     a : ({self.A_a:3d}, {self.Z_a:3d})")
    print(f"  Elab : {self.Energy : 10.3f} MeV")
    print(f"    mu : {self.mu: 10.3f} MeV/c2")
#    print(f"   Ecm : {self.Ecm: 10.3f} MeV")
#    print(f"     k : {self.k: 10.3f} MeV/c")
#    print(f"   eta : {self.eta: 10.3f}")
#    print(f"     L : {self.L},  maxL : {self.maxL}")
    print(f"    dr : {self.dr} fm, nStep : {self.nStep}")
    print(f"rStart : {self.rStart} fm,  rMax : {self.nStep * self.dr} fm")
    print(f"spin-A : {self.sA}, spin-a : {self.sa} ")


  def __init__(self, A, ZA, a, Za, Energy):
    self.A_A = A
    self.A_a = a
    self.Z_A = ZA
    self.Z_a = Za
    self.Z = ZA * Za

    self.sA = 0
    self.sa = 0

    self.L = 0
    self.S = 0
    self.J = 0

    haha = IsotopeClass()
    self.mass_A = haha.GetMassFromAZ(self.A_A, self.Z_A)
    self.mass_a = haha.GetMassFromAZ(self.A_a, self.Z_a)

    #self.mu = (A * a)/(A + a) * self.amu
    self.mu = (self.mass_A * self.mass_a)/(self.mass_A + self.mass_a)

    self.Energy = Energy
#    self.E_tot = math.sqrt(math.pow((a+A)*self.amu,2) + 2 * A * self.amu * eng_Lab)
#    self.Ecm = self.E_tot - (a + A) * self.amu
#    self.k = math.sqrt(self.mu * 2 * abs(self.Ecm)) / self.hbarc
#    self.eta = self.Z * self.ee * math.sqrt( self.mu/2/self.Ecm ) / self.hbarc

#    self.maxL = int(self.k * (1.4 * (self.A**(1/3) + self.a**(1/3)) + 3))

  def SetSpin(self, sA, sa):
    self.sA = sA
    self.sa = sa
    self.S = self.sa

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
    self.SolU = []
    self.maxSolU = 0.0

  def ClearPotential(self):
    self.potential_List = []

  def AddPotential(self, pot, useBothMass : bool = False):
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
      if pot.id == 2 and self.L > 0:
        value = value + self.LS() * pot.output(x)
      else:
        value = value + pot.output(x)
    return value

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
    self.SolU = [self.solu0]
    dSolU = [self.dsolu0]

    dyy = np.array([1., 0., 0., 0., 0.], dtype= complex)
    dzz = np.array([1., 0., 0., 0., 0.], dtype= complex)

    self.maxSolU = 0.0

    for i in range(self.nStep-1):
      r = self.rStart + self.dr * i
      y = self.SolU[i]
      z = dSolU[i]

      for j in range(4):
        dyy[j + 1] = self.dr * (z + self.parC[j] * dzz[j])
        dzz[j + 1] = self.dr * self.__G(r + self.parC[j] * self.dr, y + self.parC[j] * dyy[j], z + self.parC[j] * dzz[j])

      dy = sum(self.parD[j] * dyy[j + 1] for j in range(4))
      dz = sum(self.parD[j] * dzz[j + 1] for j in range(4))

      self.SolU.append(y + dy)
      dSolU.append(z + dz)

      if abs(self.SolU[-1]) > self.maxSolU:
        self.maxSolU = abs(self.SolU[-1])

    return self.SolU
  
  def NearestPosIndex(self, r):
    return min(len(self.rpos)-1, int((r - self.rStart) / self.dr))

