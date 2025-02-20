#!/usr/bin/env python3

from boundState import BoundState
from solveSE import WoodsSaxonPot, CoulombPotential, SpinOrbit_Pot, WS_SurfacePot
import matplotlib.pyplot as plt

# boundState = BoundState(16, 8, 1, 0, 1, 0, 0.5, -4.14)
# boundState.SetPotential(1.10, 0.65, -6, 1.25, 0.65, 1.25)
# boundState.FindPotentialDepth(-75, -60, 0.1)
# # boundState.PrintWF()
# boundState.PlotBoundState()

from distortedWave import DistortedWave

# dw = DistortedWave("60Ni", "p", 30)
# dw.ClearPotential()
# dw.AddPotential(WoodsSaxonPot(-47.937-2.853j, 1.20, 0.669), False)
# dw.AddPotential(WS_SurfacePot(-6.878j, 1.28, 0.550), False)
# dw.AddPotential(SpinOrbit_Pot(-5.250 + 0.162j, 1.02, 0.590), False)
# dw.AddPotential(CoulombPotential(1.258), False)

dw = DistortedWave("60Ni", "d", 60)
dw.PrintInput()
dw.ClearPotential()
dw.AddPotential(WoodsSaxonPot(-81.919, 1.15, 0.768), False)
dw.AddPotential(WoodsSaxonPot(-4.836j, 1.33, 0.464), False)
dw.AddPotential(WS_SurfacePot(-8.994j, 1.373, 0.774), False)
dw.AddPotential(SpinOrbit_Pot(-3.557, 0.972, 1.011), False)
dw.AddPotential(CoulombPotential(1.303), False)




dw.CalScatteringMatrix()
dw.PrintScatteringMatrix()

dw.PlotDCSUnpolarized(180, 1)

# for i in range(1, 19):
#   theta = 10*i
#   # ruth = dw.RutherFord(theta)
#   # coulAmp = dw.CoulombScatterintAmp(theta)
#   dw.CalLegendre(theta)
#   nuAmp1 = dw.NuclearScatteringAmp(-0.5, 0.5, 14)
#   nuAmp2 = dw.NuclearScatteringAmp(0.5, -0.5, 14)
#   # dsc = dw.DCSUnpolarized(theta, 14)
#   # print(f"{theta:3.0f}, {nuAmp1:15.5f}, {nuAmp2:15.5f}, {dsc:10.6f}, {ruth:10.6f}")
#   print(f"{theta:3.0f}, {nuAmp1:15.5f}, {nuAmp2:15.5f}")
