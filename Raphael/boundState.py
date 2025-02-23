#!/usr/bin/env python3

from solveSE import WoodsSaxonPot, CoulombPotential, SpinOrbit_Pot,  SolvingSE

from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
from scipy.integrate import simpson
import numpy as np
import matplotlib.pyplot as plt
from mpmath import whitw
import mpmath
import time

mpmath.mp.dps = 15  # Decimal places of precision

class BoundState(SolvingSE):
  def __init__(self, A, ZA, a, Za, node,  L, J, BE):
    super().__init__(A, ZA, a, Za, BE)
    self.CalCMConstants( True)
    self.SetRange(0, 0.1, 300) # default range
    self.SetLJ(L, J)
    self.PrintInput()
    self.node = node # number of nodes of the wave function r > 0
    self.FoundBounfState = False
    self.wf = None

  def SetPotential(self, r0, a0, Vso, rso, aso, rc = 0.0):
    self.r0 = r0
    self.a0 = a0
    self.Vso = Vso
    self.rso = rso
    self.aso = aso
    self.rc = rc
    self.ClearPotential()
    self.AddPotential(WoodsSaxonPot(-60, r0, a0), False) # not use mass number of a
    self.AddPotential(SpinOrbit_Pot(Vso, rso, aso), False) # not use mass number of a
    if rc > 0 and self.Z_a > 0:
      self.AddPotential(CoulombPotential(rc), False) # not use mass number of a

  def FindPotentialDepth(self, Vmin, Vmax, Vstep=1, isPathWhittaker = True):
    start_time = time.time()  # Start the timer
    V0List = np.arange(Vmin, Vmax, Vstep)
    lastSolU = []
    minLastSolU = 0
    maxLastSolU = 0
    for V0 in V0List:
      self.ClearPotential()
      self.AddPotential(WoodsSaxonPot(V0, self.r0, self.a0), False)
      self.AddPotential(SpinOrbit_Pot(self.Vso, self.rso, self.aso), False) # not use mass number of a
      if self.rc > 0 and self.Z_a > 0:
        self.AddPotential(CoulombPotential(self.rc), False)
      wf = self.SolveByRK4()
      lastSolU.append(np.real(wf[-1]))
      if lastSolU[-1] < minLastSolU:
        minLastSolU = lastSolU[-1]
      if lastSolU[-1] > maxLastSolU:
        maxLastSolU = lastSolU[-1]

    if minLastSolU >= 0 or maxLastSolU <= 0:
      print("No bound state found in the range")
      print("min Last SolU = ", minLastSolU)
      print("max Last SolU = ", maxLastSolU)
      plt.plot(V0List, lastSolU, marker="x")
      plt.grid()
      plt.xlabel('V0 (MeV)')
      plt.show(block=False)
      input("Press Enter to exit.")

      return

    f = interp1d(lastSolU, V0List, kind='cubic')
    self.V_BS = f(0)
    print(f"Potential Depth = {self.V_BS:15.11f} MeV")

    # recaluclate the wave function for V_BS
    self.ClearPotential()
    self.AddPotential(WoodsSaxonPot(self.V_BS, self.r0, self.a0), False)
    self.AddPotential(SpinOrbit_Pot(self.Vso, self.rso, self.aso), False) # not use mass number of a
    if self.rc > 0 and self.Z_a > 0:
      self.AddPotential(CoulombPotential(self.rc), False)

    self.SolveByRK4()
    self.solU = np.real(self.solU)

    # find number of node in self.SolU with in to potenital range
    # find how many time self.SolU change sign
    detNode = 0
    for i, r  in enumerate(self.rpos):
      if r > self.r0 * (pow(self.A_A, 1/3) + pow(self.A_a, 1/3)) + 5 * self.a0:
        break
      if self.solU[i] * self.solU[i-1] < 0:
        detNode += 1

    if detNode != self.node:
      print(f"\033[91mNumber of node is not matched, expected : {self.node}, found : {detNode}\033[0m")
      return

    self.FoundBounfState = True

    # normalize the wave function
    norm = simpson(self.solU**2) * self.dr
    self.solU /= np.sqrt(norm)
    self.wf = np.zeros_like(self.solU)
    self.wf[1:] = self.solU[1:] / self.rpos[1:]
    self.wf[0] = self.solU[0]  # Handle the first element separately if needed

    #extrapolate wf with quadrotic from 0.2, 0.1 to 0.0
    def func(x, a, b, c):
        return a * x**2 + b * x + c
    popt, pcov = curve_fit(func, self.rpos[1:6], self.wf[1:6])
    self.wf[0] = func(0, *popt)
    
    if isPathWhittaker:
      # patch the long range function to Whittaker function
      R = min(self.r0 * (pow(self.A_A, 1/3) + pow(self.A_a, 1/3)) + 5 * self.a0, max(self.rpos))
      rIndex = self.NearestPosIndex(R)
      print(f"replacing wave function with Whittaker after R : {R:5.1f}, index : {rIndex}")
      if rIndex < len(self.rpos):
        R = self.rpos[rIndex]
        W_values = [np.real(whitw(-1j*self.eta, self.L + 0.5,  2*self.k * r))/(self.k * r) for r in self.rpos[1:]]
        W_values.insert(0, 0)

        self.ANC = float(self.wf[rIndex] / W_values[rIndex])
        print(f"ANC : {self.ANC:10.6e}")
        self.wf[rIndex:] = self.ANC * np.array(W_values[rIndex:])

    end_time = time.time()  # End the timer
    print(f"Finding Potential Depth and Bound state took {(end_time - start_time) * 1000:.2f} milliseconds")

  def GetBoundStateWF(self):
    return self.wf

  def PlotBoundState(self, maker=None):
    if not self.FoundBounfState:
      plt.plot(self.rpos[1:], self.solU[1:]/self.rpos[1:])
    else:
      if maker is None:
        plt.plot(self.rpos, self.wf)
      else:
        plt.plot(self.rpos, self.wf, marker=maker)
    plt.grid()
    plt.xlabel('r (fm)')
    plt.ylabel('u(r)/r')
    plt.show(block=False)
    input("Press Enter to exit.")

  def PrintWF(self, func = None):
    if func is None:
      for r, wf in zip(self.rpos, self.wf):
        print(f"{r:5.1f}, {wf:10.6e}")
    else:
      for r, wf in zip(self.rpos, func):
        print(f"{r:5.1f}, {wf:10.6e})")
    
