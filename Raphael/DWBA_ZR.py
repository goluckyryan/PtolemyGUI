#!/usr/bin/env python3

from solveSE import WS, Coulomb, SO, WSSurface, SolvingSE


boundState = SolvingSE(16, 8, 1, 0, -4.14)
boundState.SetRange(0, 0.1, 300)
boundState.PrintInput()

boundState.ClearPotential()
boundState.AddPotential(WS(-40, 1.10, 0.65))

boundState.SetLJ(0, 0.5)

print(boundState.GetPotentialValue(0))