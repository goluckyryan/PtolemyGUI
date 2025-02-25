# Background

This is a python implementation of Distorted Wave Born Approximation calculation for direct nuclear reaction. The theory is stated on [here](https://nukephysik101.wordpress.com/2025/02/21/handbook-of-direct-nuclear-reaction-for-retarded-theorist-v1/)

The motivation is two folds: 1) for my own curiousity, 2) All of DWBA code I know of is using Fortran and some of the code cannot be compiled anymore, and for those still can be compiled, the input card is a mess. So, i decided to have a mordern version of DWBA code. 

# Requirement
* numpy
* scipy - for gamma function, curve_fit, interp1d, simpson
* matplotlib
* mpmath - for coulombf, couombg and whittaker
* sympy - for winger_9j

the version should not matter. 

# Start

Open the DWBA.py, there are some examples 

# Code Components

The foundation of the code base are 
* assLegendreP.py - for associate Legendre polynomial for positive m
* clebschGordan.py - for custom build CG, which is faster
* opticalPotential.py - for optical potential, only have An & Cai for deuteron and Kronning for proton now
* ../Cleopatra/IAEANuclearData.py - for getting nuclear data like mass and spin-partiy
* coulombWave.py - attemp to make fast CoulombWave....

Next, we have 
* solveSE.py 

which have __PotentialFrom__ class, from this class, we have WoodsSaxon, WSSurface, Spin-Orbit, and Coulomb potential. Also, the __SolveSE__ class can get the mass from IAEA, calculated the Kinematics, and Solve the Radial wave function. 

based on the __SolveSE__ class, we have two classes
* boundState.py - for solving boudstate wave function
* distortedWave.py - for solving the distorted wave and extract the Scattering Matrix

the __DistortedWave__ class can also calculated the elastics differential cross section in a ratio with RutherFord scattering. However, the code is not optimized and it may take minutes to have the angualr distribution.

Finally, we have the 
* dwba_zr.py

which is a Zero-raneg DWBA class for calculating the radial integral and do the argular summation to obtain the differential cross section for transfer reaction. The code is kind of optimized. The optimization is done by precalculated the radial integral, the associated Legendre polaynomail, the Clebsch-Gordon coefficients, and the 9j-symbol. There may be some small tweets can be done, but it is python, no machines code level optimization.

# Performanace

for s-orbtial transfer, it takes 10 secs. for d-orbital, it takes like half minute.

# Limitation & future work
* This is only for zero-range approximation, the angular distribution agree with experiment.
* Only for (d,p) or (p,d) (not and need test), should be work for all single nucleon transfer reaction. 

Need to 
* work for inelastic scattering
* implement the Finite-Range calcualtion
* polarization
* work on the 2-nucleons transfers
* work on the Coupled-Channels
* have a C++ version for speed sake (or use multiple cores, or both)
