#!/usr/bin/env python3

import sys

filename = sys.argv[1]

if len(sys.argv) < 2:
    print("Error: Not enough arguments provided.")
    print("Usage: ./{sys.argv[0]} filename")
    sys.exit(1)


import numpy as np


# Load data while skipping the first two lines
data = np.loadtxt(filename, skiprows=2)

# Extract columns
angles, values = data[:, 0], data[:, 1]

for a, b in zip(angles, values):
  print(f"{{{a:5.1f}, {b:7.3f}}},", end="")

print("\n")