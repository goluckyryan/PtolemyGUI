#!/usr/bin/env python3

import numpy as np
from mpmath import coulombf, coulombg

L = 20
eta = 2.0
rho = 30

f_values = coulombf(L, eta, rho)
g_values = coulombg(L, eta, rho)

print(f_values)
print(g_values)