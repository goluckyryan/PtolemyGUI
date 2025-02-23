#!/usr/bin/env python3

import time
import matplotlib.pyplot as plt
from dwba_zr import DWBA_ZR

haha = DWBA_ZR("16O", "d", "p", "17O", "1/2+", "1s1/2", 0.0, 10)
haha.FindBoundState()
haha.CalRadialIntegral()

# haha.PrintRadialIntegral()

# haha.boundState.PlotBoundState()

# haha.PlotRadialIntegral()

# haha.PlotDistortedWave(True, 2, 1, 20)
# haha.PlotDistortedWave(False, 2, 1.5, 20)

haha.CalAngDistribution()
haha.PlotAngDist()


