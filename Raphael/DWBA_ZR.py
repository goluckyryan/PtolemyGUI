#!/usr/bin/env python3

import time
from dwba_zr import DWBA_ZR

haha = DWBA_ZR("16O", "d", "p", "17O", "1/2+", "1s1/2", 0.0, 10)

# haha.boundState.PlotBoundState()

haha.PrintRadialIntegral()
# haha.PlotRadialIntegral()

# haha.PlotDistortedWave(True, 2, 1, 20)
# haha.PlotDistortedWave(False, 2, 1.5, 20)

j = 1/2
sa = 1
sb = 1/2

# hehe = haha.Gamma(0, 1, 0, 0.5, 0, 1, 0.5)
# print(hehe)

# haha.PreCalLegendreP(10)
# print(haha.legendrePArray)

# jaja = haha.Beta(-1, 1, -0.5) * haha.ffactor
# print(jaja)

# lala = haha.AngDist(10)
# print(lala)

angList = []
xsec = []

start_time = time.time()
for i in range(0, 181, 5):
  angList.append(i)
  kaka = haha.AngDist(i)
  xsec.append(kaka)
  print(i, kaka)

stop_time = time.time()
print(f"Total time {(stop_time - start_time) :.2f} sec")

import matplotlib.pyplot as plt

plt.plot(angList, xsec)
plt.xlim(-1, 181)
plt.yscale("log")
plt.grid()
plt.show(block=False)

input("Press Enter to continue...")
