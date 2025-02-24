#!/usr/bin/env python3

import time
import matplotlib.pyplot as plt
from dwba_zr import DWBA_ZR

haha = DWBA_ZR("16O", "d", "p", "17O", "1/2+", "1s1/2", 0.0, 10)
# haha = DWBA_ZR("16O", "d", "p", "17O", "5/2+", "0d5/2", 0.0, 10)
haha.FindBoundState()

# haha.boundState.PlotBoundState()

haha.CalRadialIntegral()

# haha.PlotScatteringMatrix(False)

# haha.PlotIncomingDistortedWave(2, 1, 20)
# haha.PlotOutgoingDistortedWave(15, 15.5, 20) 

# haha.PrintRadialIntegral()
# haha.PlotRadialIntegral()
# haha.PlotRadialIntegralSigle(1, 1, -0.5)


haha.CalAngDistribution(0, 180, 1)
haha.PrintAngDist()
haha.PlotAngDist()


######################################

# print(haha.GetPreCalNineJ(1, 1, 3, 3.5))

# print(haha.Gamma(1, 1, 3, 3.5, 0, 1, 0.5))



# haha.PreCalLegendreP(5)
# haha.Beta(0, 1, 0.5)
