#!/usr/bin/env python3

import time
import matplotlib.pyplot as plt
from distortedWave import WoodsSaxonPot, CoulombPotential, SpinOrbit_Pot, WS_SurfacePot, DistortedWave
from dwba_zr import DWBA_ZR

# haha = DWBA_ZR("16O", "d", "p", "17O", "1/2+", "1s1/2", 0.0, 10)
# haha = DWBA_ZR("16O", "d", "p", "17O", "5/2+", "0d5/2", 0.0, 10)


# haha.FindBoundState()
# haha.boundState.PlotBoundState()

# haha.CalRadialIntegral()

# haha.PlotScatteringMatrix(False)

# haha.PlotIncomingDistortedWave(2, 1, 20)
# haha.PlotOutgoingDistortedWave(15, 15.5, 20) 

# haha.PrintRadialIntegral()
# haha.PlotRadialIntegral()
# haha.PlotRadialIntegralSigle(1, 1, -0.5)


# haha.CalAngDistribution(0, 180, 1)
# haha.PrintAngDist()
# haha.PlotAngDist()


####################################### Simple distorted wave calculation
kaka = DistortedWave("11C", "p", 60)
kaka.SetRange(0, 0.01, 2000)
kaka.maxL = 14
kaka.PrintInput()
kaka.ClearPotential()
kaka.AddPotential(WoodsSaxonPot(-34.714+6.749j, 1.122, 0.676), False) # False = only use 11C for radius calculation
kaka.AddPotential(WS_SurfacePot(       -3.194j, 1.307, 0.524), False)
kaka.AddPotential(SpinOrbit_Pot( -4.532+0.477j, 0.894, 0.500), False)
kaka.AddPotential(CoulombPotential(1.578), False)
kaka.PrintPotentials()

kaka.CalScatteringMatrix(True)
kaka.PrintScatteringMatrix()
kaka.PlotScatteringMatrix()

# kaka.PlotDistortedWave(1, 1.5)

# kaka.PlotDCSUnpolarized(180, 1, None, True)