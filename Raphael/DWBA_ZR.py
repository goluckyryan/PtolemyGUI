#!/usr/bin/env python3

from boundState import BoundState
    
boundState = BoundState(16, 8, 1, 0, 0, 2, 2.5, -4.14)
boundState.SetPotential(1.10, 0.65, -6, 1.25, 0.65)
boundState.FindPotentialDepth(-70, -60, 0.5)
# boundState.PrintWF()
boundState.PlotBoundState()