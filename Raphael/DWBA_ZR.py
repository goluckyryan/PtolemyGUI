#!/usr/bin/env python3

from boundState import BoundState

from solveSE import WoodsSaxonPot, CoulombPotential, SpinOrbit_Pot, WS_SurfacePot, SolvingSE
from mpmath import coulombf, coulombg

import numpy as np
from scipy.special import gamma

# boundState = BoundState(16, 8, 1, 0, 1, 0, 0.5, -4.14)
# boundState.SetPotential(1.10, 0.65, -6, 1.25, 0.65, 1.25)
# boundState.FindPotentialDepth(-75, -60, 0.1)
# # boundState.PrintWF()
# boundState.PlotBoundState()

def SevenPointsSlope(data, n):
  return (-data[n + 3] + 9 * data[n + 2] - 45 * data[n + 1] + 45 * data[n - 1] - 9 * data[n - 2] + data[n - 3]) / 60

def FivePointsSlope(data, n):
  return ( data[n + 2] - 8 * data[n + 1] + 8 * data[n - 1] -  data[n - 2] ) / 12


dw = SolvingSE("60Ni", "p", 30)
dw.SetRange(0, 0.1, 300)
dw.SetLJ(0, 0.5)
dw.dsolu0 = 1
dw.CalCMConstants()
dw.PrintInput()

dw.ClearPotential()
dw.AddPotential(WoodsSaxonPot(-47.937-2.853j, 1.20, 0.669), False)
dw.AddPotential(WS_SurfacePot(-6.878j, 1.28, 0.550), False)
dw.AddPotential(SpinOrbit_Pot(-5.250 + 0.162j, 1.02, 0.590), False)
dw.AddPotential(CoulombPotential(1.258), False)

rpos = dw.rpos
solU = dw.SolveByRK4()

solU /= dw.maxSolU

# for r, u in zip(rpos, solU):
#     print(f"{r:.3f} {np.real(u):.6f} {np.imag(u):.6f}")

def CoulombPhaseShift(L, eta):
  return np.angle(gamma(L+1+1j*eta))

sigma = CoulombPhaseShift(dw.L, dw.eta)
# find pahse shift by using the asymptotic behavior of the wave function
r1 = rpos[-2]
f1 = float(coulombf(dw.L, dw.eta, dw.k*r1))
g1 = float(coulombg(dw.L, dw.eta, dw.k*r1))
u1 = solU[-2]

r2 = rpos[-1]
f2 = float(coulombf(dw.L, dw.eta, dw.k*r2))
g2 = float(coulombg(dw.L, dw.eta, dw.k*r2))
u2 = solU[-1]

det = f2*g1 - f1*g2
A = (f2*u1 - u2*f1) / det
B = (u2*g1 - g2*u1) / det

print(f"A = {np.real(A):.6f} + {np.imag(A):.6f} I")
print(f"B = {np.real(B):.6f} + {np.imag(B):.6f} I")

ScatMatrix = (B + A * 1j)/(B - A * 1j)
print(f"Scat Matrix = {np.real(ScatMatrix):.6f} + {np.imag(ScatMatrix):.6f} I")

solU = np.array(solU, dtype=np.complex128)
solU *= np.exp(1j * sigma)/(B-A*1j)

from matplotlib import pyplot as plt

plt.plot(rpos, np.real(solU), label="Real")
plt.plot(rpos, np.imag(solU), label="Imaginary")
plt.legend()
plt.show(block=False)

input("Press Enter to continue...")