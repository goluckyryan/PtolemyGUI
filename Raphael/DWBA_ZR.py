#!/usr/bin/env python3

from boundState import BoundState
from solveSE import WoodsSaxonPot, CoulombPotential, SpinOrbit_Pot, WS_SurfacePot

# boundState = BoundState(16, 8, 1, 0, 1, 0, 0.5, -4.14)
# boundState.SetPotential(1.10, 0.65, -6, 1.25, 0.65, 1.25)
# boundState.FindPotentialDepth(-75, -60, 0.1)
# # boundState.PrintWF()
# boundState.PlotBoundState()

from distortedWave import DistortedWave

dw = DistortedWave("60Ni", "p", 30)

dw.ClearPotential()
dw.AddPotential(WoodsSaxonPot(-47.937-2.853j, 1.20, 0.669), False)
dw.AddPotential(WS_SurfacePot(-6.878j, 1.28, 0.550), False)
dw.AddPotential(SpinOrbit_Pot(-5.250 + 0.162j, 1.02, 0.590), False)
dw.AddPotential(CoulombPotential(1.258), False)

dw.CalScatteringMatrix()

dw.PrintScatteringMatrix()

# dw.PlotScatteringMatrix()

dw.PlotDCSUnpolarized()