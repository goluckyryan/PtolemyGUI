#!/usr/bin/env python3

from boundState import BoundState

from solveSE import WoodsSaxonPot, CoulombPotential, SpinOrbit_Pot, WS_SurfacePot, SolvingSE
from mpmath import coulombf, coulombg

import numpy as np
from scipy.special import gamma, factorial
import matplotlib.pyplot as plt

from assLegendreP import associated_legendre_array

def SevenPointsSlope(data, n):
  return (-data[n + 3] + 9 * data[n + 2] - 45 * data[n + 1] + 45 * data[n - 1] - 9 * data[n - 2] + data[n - 3]) / 60

def FivePointsSlope(data, n):
  return ( data[n + 2] - 8 * data[n + 1] + 8 * data[n - 1] -  data[n - 2] ) / 12

from clebschGordan import clebsch_gordan, KroneckerDelta
import time

############################################################
class DistortedWave(SolvingSE):
  def __init__(self, target, projectile, ELab):
    super().__init__(target, projectile, ELab)
    self.SetRange(0, 0.1, 300)
    self.CalCMConstants()

    self.ScatMatrix = []
    self.distortedWaveU = []

    self.legendreArray = []

  def SetLJ(self, L, J):
    self.L = L
    self.J = J
    self.dsolu0 = pow(0.1, 2*L+1)

  def CoulombPhaseShift(self, L = None, eta = None):
    if L is None:
      L = self.L
    if eta is None:
      eta = self.eta
    return np.angle(gamma(L+1+1j*eta))
  
  def CoulombHankel(self, sign, rho, L = None, eta = None) -> complex:
    if L is None:
      L = self.L
    if eta is None:
      eta = self.eta
    return float(coulombg(L, eta, rho)) + sign * 1j * float(coulombf(L, eta, rho)) 

  def CalScatteringMatrix(self, normTo1 = False, maxL = None, verbose = False):
    start_time = time.time()  # Start the timer

    if maxL is None:
      maxL = self.maxL

    self.ScatMatrix = []
    self.distortedWaveU = []

    tempZeroList = np.zeros(len(self.rpos), dtype=np.complex128)

    for L in range(0, maxL+1):

      temp_ScatMatrix = []
      temp_distortedWaveU = []

      for J in np.arange(L-self.S, L + self.S+1, 1):
        if J < 0 or (L==0 and J != self.S):
          temp_ScatMatrix.append(0)
          temp_distortedWaveU.append(tempZeroList)
          continue

        self.SetLJ(L, J)
        self.SolveByRK4()

        #============ 2-point using Hankel
        # r1 = self.rpos[-2]
        # u1 = self.solU[-2]
        # Hp1 = self.CoulombHankel( 1, self.k * r1)
        # Hm1 = self.CoulombHankel(-1, self.k * r1)

        # r2 = self.rpos[-1]
        # u2 = self.solU[-1]
        # Hp2 = self.CoulombHankel( 1, self.k * r2)
        # Hm2 = self.CoulombHankel(-1, self.k * r2)

        # ScatMatrix = (u1 * Hm2 - Hm1 * u2) / (u1 * Hp2 - Hp1 * u2)
        # norm = 2j * (u1 * Hp2 - Hp1 * u2 ) / (Hp1 * Hm2 - Hm1 * Hp2)

        #=========== 2 point G, F, this is fastest
        r1 = self.rpos[-2]
        u1 = self.solU[-2]
        f1 = float(coulombf(self.L, self.eta, self.k*r1))
        g1 = float(coulombg(self.L, self.eta, self.k*r1))

        r2 = self.rpos[-1]
        u2 = self.solU[-1]
        f2 = float(coulombf(self.L, self.eta, self.k*r2))
        g2 = float(coulombg(self.L, self.eta, self.k*r2))

        print(f"{L:2d}, {J:4.1f} | {r1:.3f}, {f1:10.6f}, {g1:10.6f} | {r2:.3f}, {f2:10.6f}, {g2:10.6f}")
        
        det = f2*g1 - f1*g2
        A = (f2*u1 - u2*f1) / det
        B = (u2*g1 - g2*u1) / det
        ScatMatrix = (B + A * 1j)/(B - A * 1j)
        norm = B - A * 1j

        #============ 5 points slopes
        # r_list = self.rpos[-5:]
        # u_list = self.solU[-5:]
        # f_list = [float(coulombf(self.L, self.eta, self.k * r)) for r in r_list]
        # g_list = [float(coulombg(self.L, self.eta, self.k * r)) for r in r_list]
        # u0 = u_list[2]
        # f0 = f_list[2]
        # g0 = g_list[2]
        # du = FivePointsSlope(u_list, 2)
        # df = FivePointsSlope(f_list, 2)
        # dg = FivePointsSlope(g_list, 2)

        # det = df*g0 - dg*f0
        # A = (df*u0 - du*f0) /det
        # B = (du*g0 - dg*u0) /det
        # ScatMatrix = (B + A * 1j)/(B - A * 1j)
        # norm = B - A * 1j

        # #============ 100 points fitting
        # r_list = self.rpos[-5:]
        # u_list = self.solU[-5:]
        # f_list = [float(coulombf(self.L, self.eta, self.k * r)) for r in r_list]
        # g_list = [float(coulombg(self.L, self.eta, self.k * r)) for r in r_list]

        # # Fit u_list to A * g_list + B * f_list
        # A, B = np.linalg.lstsq(np.vstack([g_list, f_list]).T, u_list, rcond=None)[0]
        # ScatMatrix = (B + A * 1j) / (B - A * 1j)
        # norm = B - A * 1j

        if verbose:
          print(f"{{{L},{J}, {np.real(ScatMatrix):10.6f} +  {np.imag(ScatMatrix):10.6f}I}}")

        temp_ScatMatrix.append(ScatMatrix)

        dwU = np.array(self.solU, dtype=np.complex128)
        if normTo1 :
          dwU /= self.maxSolU
        else:
          dwU *= 1./norm
        temp_distortedWaveU.append(dwU)
      
      self.ScatMatrix.append(temp_ScatMatrix)
      self.distortedWaveU.append(temp_distortedWaveU)

    end_time = time.time()  # End the timer
    print(f"Calculate Scattering Matrixes took {(end_time - start_time) * 1000:.2f} milliseconds")

    return [self.ScatMatrix, self.distortedWaveU]

  def PrintScatteringMatrix(self):
    print("======================= Scattering Matrix")
    for L in range(0, len(self.ScatMatrix)):

      if L == 0 :
        print(" ", end="")
        for i in range(0, len(self.ScatMatrix[L])):
            print(f"{{{'L':>2s},{'J':>4s}, {'Real':>10s} +   {'Imaginary':>10s}}}, ", end="")
        print("")

      print("{", end="")
      for i in range(0, len(self.ScatMatrix[L])):
        print("{", end="")
        print(f"{L:2d},{L+i-self.S:4.1f}, {np.real(self.ScatMatrix[L][i]):10.6f} +  {np.imag(self.ScatMatrix[L][i]):10.6f}I", end="")
        if i < len(self.ScatMatrix[L])-1 :
          print("}, ", end="")
        else:
          print("}", end="")
      print("},")

  def GetScatteringMatrix(self, L, J):
    return self.ScatMatrix[L][J-L+self.S]

  def GetDistortedWave(self, L, J):
    return self.distortedWaveU[L][int(J-L+self.S)]

  def PlotDistortedWave(self, L, J, maxR = None):
    plt.plot(self.rpos, np.real(self.GetDistortedWave(L, J)), label="Real")
    plt.plot(self.rpos, np.imag(self.GetDistortedWave(L, J)), label="Imaginary")
    plt.title(f"Radial wave function for L={L} and J={J}")
    if maxR != None :
      plt.xlim(-1, maxR) 
    plt.legend()
    plt.grid()
    plt.show(block=False)
    input("Press Enter to continue...")

  def PlotScatteringMatrix(self):
    nSpin = int(self.S*2+1)
    fig, axes = plt.subplots(1, nSpin, figsize=(6*nSpin, 4))

    for i in range(0, nSpin):
      sm = []
      l_list = []
      for L in range(0, len(self.ScatMatrix)):
        if i == 0 and L == 0 :
          continue
        l_list.append(L)
        sm.append(self.ScatMatrix[L][i])

      axes[i].plot(l_list, np.real(sm), label="Real", marker='o')
      axes[i].plot(l_list, np.imag(sm), label="Imaginary", marker='x')
      axes[i].legend()
      axes[i].set_xlabel('L')
      axes[i].set_xlim(-1, self.maxL+1)
      axes[i].set_xticks(np.arange(0, self.maxL + 1, 2))
      axes[i].set_ylabel('Value')
      if self.S*2 % 2 == 0 :
        str = f'{int(i-self.S):+d}'
      else:
        str = f'{int(2*(i-self.S)):+d}/2' 
      axes[i].set_title(f'Real and Imaginary Parts vs L for Spin J = L{str}')
      axes[i].grid()

    plt.tight_layout()
    plt.show(block=False)
    input("Press Enter to continue...")

  def RutherFord(self, theta_deg):
    sin_half_theta = np.sin(np.radians(theta_deg + 1e-20) / 2)
    result = self.eta**2 / (4 * (self.k**2) * (sin_half_theta**4))
    return result

  def CoulombScatterintAmp(self, theta_deg):
    sin_sq = pow(np.sin(np.radians(theta_deg + 1e-20)/2), 2)
    coulPS = self.CoulombPhaseShift(0)
    return - self.eta / (2 * self.k * sin_sq) * np.exp(1j * (2*coulPS - self.eta * np.log(sin_sq)))
  
  def GMatrix1Spin(self, v, v0, l ) -> complex:
    if self.S == 0 :
      return self.ScatMatrix[l][0] - KroneckerDelta(v, v0)
    else:
      Jmin = l - self.S
      Jmax = l + self.S
      value = 0
      for J in  np.arange(Jmin, Jmax + 1, 1):
        index = int(J - Jmin)
        cg1 = clebsch_gordan(l, 0, self.S, v0, J, v0)
        cg2 = clebsch_gordan(l, v0 - v, self.S, v, J, v0)
        value += cg1 * cg2 * self.ScatMatrix[l][index]
      return value - KroneckerDelta(v, v0)

  def CalLegendre(self, theta_deg, maxL = None, maxM = None):
    if maxL is None:
      maxL = self.maxL
    if maxM is None:
      maxM = int(2*self.S)
    self.legendreArray = associated_legendre_array(maxL, maxM, theta_deg)
    return self.legendreArray

  def GetPreCalLegendre(self, L, M):
    if abs (M) <= int(2*self.S):
      return self.legendreArray[L][int(abs(M))]
    else :
      return 0

  def NuclearScatteringAmp(self, v, v0, maxL = None ) -> complex:
    value = 0
    if maxL is None:
      maxL = self.maxL
    for l in range(0, maxL+1):
      if abs(v0-v) > l :
        value += 0
      else:
        coulPS = self.CoulombPhaseShift(l)
        fact = pow(-1, v0-v) * np.sqrt(factorial(l - abs(v0-v))/factorial(l + abs(v0-v)))
        value += (2*l+1) * fact * self.GetPreCalLegendre(l, v0-v) * np.exp(2j * coulPS)* self.GMatrix1Spin(v, v0, l)

    return value / 2j / self.k
  
  def DCSUnpolarized(self, theta_deg, maxL = None):
    value = 0
    self.CalLegendre(theta_deg)
    jaja = self.CoulombScatterintAmp(theta_deg)
    for v in np.arange(-self.S, self.S + 1, 1):
      for v0 in np.arange(-self.S, self.S + 1, 1):
        value += abs( jaja * KroneckerDelta(v, v0) +  self.NuclearScatteringAmp(v, v0, maxL))**2

    value = value / (2 * self.S + 1) 
    return value

  def CalAngDistribution(self,  thetaRange = 180, thetaStepDeg = 0.2, maxL = None, verbose = False):
    self.theta_list = np.linspace(0, thetaRange, int(thetaRange/thetaStepDeg)+1)
    self.angDist = []
    for theta in self.theta_list:
      if theta == 0:
        self.angDist.append(1)
      else:
        self.angDist.append(self.DCSUnpolarized(theta, maxL)/ self.RutherFord(theta))
        if verbose :
          print(f"{theta:6.2f}, {self.angDist[-1]:10.6f}")

  def PrintAngDistribution(self):
    print("======================= Angular Distribution")
    for theta, dist in zip(self.theta_list, self.angDist):
      print(f"{theta:5.1f}, {dist:10.3f}")

  def PlotDCSUnpolarized(self):
    thetaTick = 30
    thetaRange = max(self.theta_list)
    if thetaRange < 180:
      thetaTick = 10

    plt.figure(figsize=(8, 6))
    # plt.plot(theta_values, y_values, marker='o', linestyle='-', color='blue')    
    plt.plot(self.theta_list, self.angDist, linestyle='-', color='blue')
    plt.title("Differential Cross Section (Unpolarized)")
    plt.xlabel("Angle [deg]")
    plt.ylabel("D.C.S / Ruth")
    plt.yscale("log")
    plt.xticks(np.arange(0, thetaRange + 1, thetaTick))
    plt.xlim(0, thetaRange)
    plt.grid()
    plt.show(block=False)
    input("Press Enter to continue...")