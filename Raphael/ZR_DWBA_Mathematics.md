# Zero-Range DWBA: Mathematical Formulation and Implementation

This document details the mathematical formulation and computational implementation of the Zero-Range Distorted Wave Born Approximation (ZR-DWBA) for direct nuclear reactions, as implemented in the `dwba_zr.py` code.

## Table of Contents
1. [Introduction](#introduction)
2. [Theoretical Framework](#theoretical-framework)
3. [Mathematical Formulation](#mathematical-formulation)
4. [Computational Implementation](#computational-implementation)
5. [Optimization Strategies](#optimization-strategies)
6. [Code Structure](#code-structure)

---

## Introduction

The Zero-Range DWBA calculates differential cross sections for single-nucleon transfer reactions such as (d,p) and (p,d). The "zero-range" approximation assumes the interaction between transferred nucleon and projectile/ejectile occurs at a single point, simplifying the radial integrals while maintaining good agreement with experimental angular distributions.

### Reaction Schema

For a transfer reaction: **A(a,b)B**
- **A**: Target nucleus (mass A_A, charge Z_A, spin J_A)
- **a**: Projectile (mass A_a, charge Z_a, spin s_a)
- **b**: Ejectile (mass A_b, charge Z_b, spin s_b)
- **B**: Residual nucleus (mass A_B, charge Z_B, spin J_B)
- **x**: Transferred nucleon (binding energy BE, orbital (n,l,j))

Conservation: `A + a = B + b` and `A_A + A_a = A_B + A_b`

---

## Theoretical Framework

### DWBA Transition Amplitude

The DWBA transition amplitude for a transfer reaction is:

$$T_{fi} = \langle \chi_b^{(-)}(\mathbf{k}_b) \Phi_B | V | \Phi_A \chi_a^{(+)}(\mathbf{k}_a) \rangle$$

Where:
- $\chi_a^{(+)}$: Incoming distorted wave (scattering wave function for $a+A$)
- $\chi_b^{(-)}$: Outgoing distorted wave (time-reversed scattering for $b+B$)
- $\Phi_A, \Phi_B$: Internal wave functions of target and residual
- $V$: Interaction potential (zero-range approximation)

### Zero-Range Approximation

In the zero-range limit, the interaction is approximated as:

$$V(\mathbf{r}_a, \mathbf{r}_b, \mathbf{r}_x) \rightarrow D_0 \, \delta(\mathbf{r}_a - \mathbf{r}_x) \, \delta(\mathbf{r}_b - \mathbf{r}_x)$$

Where:
- $D_0$: Zero-range normalization constant ($D_0 = 1.55 \times 10^4$ MeV·fm³ for (d,p))
- $\mathbf{r}_a, \mathbf{r}_b, \mathbf{r}_x$: Position vectors of projectile, ejectile, and transferred nucleon

This simplifies the transition matrix element to a product of:
1. **Radial integral**: Overlap of bound state with distorted waves
2. **Angular coupling coefficients**: Clebsch-Gordan and 9j-symbols

---

## Mathematical Formulation

### Differential Cross Section

The differential cross section is:

$$\frac{d\sigma}{d\Omega} = \frac{\mu_I \mu_O}{2\pi \hbar^4} \cdot \frac{k_O}{k_I^3} \cdot S \cdot \sum |T_{fi}|^2$$

**Implementation in code** ([dwba_zr.py:146-147](file:///home/ryan/PtolemyGUI/Raphael/dwba_zr.py#L146-L147)):
```python
self.xsecScalingfactor = mass_I * mass_O / np.pi / self.dwI.hbarc**4 / k_I**3 / k_O 
                         * self.spinFactor * self.transMatrixFactor
```

Where:
- $\mu_I, \mu_O$: Reduced masses (incoming/outgoing channels)
- $k_I, k_O$: Wave numbers (incoming/outgoing channels)
- $S$: Spectroscopic factor (encoded in `transMatrixFactor` and `spinFactor`)
- **Spin factor**: $\frac{1}{(2J_A+1)(2s_a+1)}$ - initial state averaging
- **Trans. Matrix factor**: $(2J_B+1) \times A_{lsj}^2$ where $A_{lsj}^2 = D_0 \times 3/2$ for (d,p)

### Transition Matrix Element Decomposition

The transition amplitude is decomposed using angular momentum coupling:

$$T_{fi} = \sum_{L_1,J_1,L_2,J_2,m,m_a,m_b} \beta(m, m_a, m_b) \times Y_{L_2}^m(\theta, \phi)$$

Where the summation runs over:
- $L_1, J_1$: Orbital and total angular momentum of incoming channel
- $L_2, J_2$: Orbital and total angular momentum of outgoing channel
- $m$: Magnetic quantum number projection
- $m_a, m_b$: Spin projections of projectile and ejectile

### Beta Coefficient

The core quantity $\beta(m, m_a, m_b)$ combines radial and angular parts:

$$\beta(m, m_a, m_b) = \sum_{L_1,J_1,L_2,J_2} \Gamma(L_1,J_1,L_2,J_2,m,m_a,m_b) \times P_{L_2}^{|m|}(\cos \theta) \times I(L_1,J_1,L_2,J_2)$$

**Implementation** ([dwba_zr.py:461-494](file:///home/ryan/PtolemyGUI/Raphael/dwba_zr.py#L461-L494)):
```python
def Beta(self, m:int, ma, mb):
  result = 0
  for L1, J1, L2, J2 in [all valid quantum number combinations]:
    gg = self.Gamma(L1, J1, L2, J2, m, ma, mb)
    lp = self.legendrePArray[L2][int(abs(m))]  # Associated Legendre polynomial
    ri = self.radialInt[L1][index1][indexL2][index2]  # Radial integral
    result += gg * lp * ri
  return result
```

### Gamma Function - Angular Coupling

The $\Gamma$ function encodes all angular momentum coupling via:

$$\begin{align}
\Gamma(L_1,J_1,L_2,J_2,m,m_a,m_b) &= (-1)^m \times i^{(L_1-L_2-l)} \times \sqrt{(2l+1)(2s+1)(2L_1+1)(2J_2+1)(2L_2+1)} \\
&\times \sqrt{\frac{(L_2-|m|)!}{(L_2+|m|)!}} \\
&\times \begin{Bmatrix} J_2 & l & J_1 \\ L_2 & s & L_1 \\ s_b & j & s_a \end{Bmatrix} \quad \text{(9j-symbol)} \\
&\times \langle J_2, m_b-m; j, m-m_b+m_a | J_1, m_a \rangle \quad \text{(CG)} \\
&\times \langle L_1, 0; s_a, m_a | J_1, m_a \rangle \quad \text{(CG)} \\
&\times \langle L_2, -m; s_b, m_b | J_2, m_b-m \rangle \quad \text{(CG)} \\
&\times \langle L_2, 0; l, 0 | L_1, 0 \rangle \quad \text{(CG)}
\end{align}$$

**Implementation** ([dwba_zr.py:444-459](file:///home/ryan/PtolemyGUI/Raphael/dwba_zr.py#L444-L459)):
```python
def Gamma(self, L1:int, J1, L2:int, J2, m:int, ma, mb):
  # Selection rule: L1 + L2 + l must be even
  if int(L1 + L2 + self.l) % 2 != 0:
    return 0
    
  fact0 = self.GetPreCalNineJ(L1, J1, L2, J2)  # 9j-symbol
  fact1 = pow(-1, m) * np.power(1j, L1-L2-self.l) * (2*L2+1) 
          * np.sqrt((2*self.l+1)*(2*self.s+1)*(2*L1+1)*(2*J2+1))
  fact2 = np.sqrt(quantum_factorial(L2-abs(m)) / quantum_factorial(L2 + abs(m)))
  fact3 = self.GetPreCalCG(J2, mb-m,      self.j, m-mb+ma, J1,   ma)
  fact4 = self.GetPreCalCG(L1,    0, self.spin_a,      ma, J1,   ma)
  fact5 = self.GetPreCalCG(L2,   -m, self.spin_b,      mb, J2, mb-m)
  fact6 = self.GetPreCalCG(L2,    0,      self.l,       0, L1,    0)
  
  return fact0 * fact1 * fact2 * fact3 * fact4 * fact5 * fact6
```

**Selection Rules** (enforced in code):
1. $L_1 + L_2 + l$ must be even (parity conservation)
2. Triangle inequalities: $|J_1 - j| \leq J_2 \leq J_1 + j$
3. $|m| \leq L_2$ (magnetic quantum number constraint)
4. Various spin coupling constraints via Clebsch-Gordan coefficients

---

## Computational Implementation

### 1. Reaction Data Processing

**Class**: `ReactionData` ([reactionData.py:21-184](file:///home/ryan/PtolemyGUI/Raphael/reactionData.py#L21-L184))

**Purpose**: Parse reaction notation and compute kinematics

**Input**: `ReactionData(nu_A, nu_a, nu_b, JB, orbital, ExB, ELabPerU, JA)`

Example: `ReactionData("40Ca", "d", "p", "3/2+", "1d5/2", 2.0, 10.0)`

**Key computations**:
1. Extract nuclear masses from IAEA database
2. Identify transferred nucleon: `x = a - b` (or `b - a`)
3. Calculate binding energy: `BE = M_B - M_A - M_x + Ex_B`
4. Parse orbital: `"1d5/2"` → node=1, l=2, j=5/2
5. Compute Q-value: `Q = M_A + M_a - M_b - M_B - Ex_B`
6. Calculate outgoing energy using relativistic kinematics
7. Verify angular momentum conservation (triangle rules)

### 2. Bound State Wave Function

**Class**: `BoundState` ([boundState.py:16-154](file:///home/ryan/PtolemyGUI/Raphael/boundState.py#L16-L154))

**Purpose**: Solve for the bound state radial wave function `u_{\nlj}(r)` of the transferred nucleon

**Method**: 
1. **Potential**: Woods-Saxon + Spin-Orbit + Coulomb
   $$V(r) = V_0 f(r) + V_{SO} \frac{1}{r}\frac{d}{dr}f(r) \, \mathbf{L}\cdot\mathbf{S} + V_C(r)$$
   $$f(r) = \frac{1}{1 + \exp\left(\frac{r - R_0}{a_0}\right)}$$

2. **Depth Search**: Vary V₀ until radial wave function has correct number of nodes and vanishes at infinity
   - Uses Simpson integration to find zero-crossing
   - Cubic interpolation to find exact V₀

3. **Normalization**: Normalize such that $\int u^2(r) \, dr = 1$

4. **Asymptotic Matching**: At large $r$, match to Whittaker function $W_{-i\eta, l+1/2}(2kr)$
   - Extracts **Asymptotic Normalization Coefficient (ANC)**

**Radial Schrödinger Equation** (solved by RK4):
$$\frac{d^2u}{dr^2} + \left[\frac{2\mu}{\hbar^2} (E - V(r)) - \frac{l(l+1)}{r^2}\right] u = 0$$

### 3. Distorted Waves

**Class**: `DistortedWave` ([distortedWave.py:22-323](file:///home/ryan/PtolemyGUI/Raphael/distortedWave.py#L22-L323))

**Purpose**: 
1. Solve for incoming/outgoing channel wave functions using optical potential
2. Extract scattering matrix elements `S_{LJ}`

**Optical Potentials**:
- **Deuterons**: An-Cai global parametrization
- **Protons**: Koning-Delaroche global parametrization

**Form**:
$$U(r) = -V f(r) - iW f_I(r) - i4W_s \frac{d}{dr}f_s(r) - V_{SO} \frac{1}{r}\frac{d}{dr}f_{SO}(r) \, \mathbf{L}\cdot\mathbf{S} + V_C(r)$$

**Scattering Matrix Extraction**:
The code matches the numerical solution to the asymptotic Coulomb-modified spherical Hankel functions:

$$u_{LJ}(r) \xrightarrow{r \to \infty} H_L^{(+)}(kr) \left[\delta_{LJ} - S_{LJ} \frac{H_L^{(-)}(kr)}{H_L^{(+)}(kr)}\right]$$

Numerically extracted by comparing logarithmic derivatives at large r.

### 4. Radial Integral Calculation

**Key function**: `CalScatMatrixAndRadialIntegral()` ([dwba_zr.py:179-271](file:///home/ryan/PtolemyGUI/Raphael/dwba_zr.py#L179-L271))

**Radial integral**:
$$I(L_1,J_1,L_2,J_2) = \int_0^\infty u_{nlj}(r) \times u_{L_1J_1}^{\text{in}}(r) \times u_{L_2J_2}^{\text{out}}(r) \times e^{i(\sigma_1+\sigma_2)} \, dr$$

Where:
- $u_{nlj}$: Bound state wave function
- $u_{L_1J_1}^{\text{in}}$: Incoming distorted wave
- $u_{L_2J_2}^{\text{out}}$: Outgoing distorted wave  
- $\sigma_1, \sigma_2$: Coulomb phase shifts

**Critical detail**: **Mass rescaling**
Since incoming and outgoing channels have different reduced masses, the radial mesh must be scaled by the mass ratio to ensure proper overlap:

```python
if self.dwI.A_a > self.dwO.A_a:  # e.g., (d,p)
  rpos_O_temp = self.rpos_I * self.massFactor  # Scale outgoing mesh
  # Interpolate outgoing wave onto incoming mesh
else:  # e.g., (p,d)
  rpos_I_temp = self.rpos_O * self.massFactor  # Scale incoming mesh
  # Interpolate incoming wave onto outgoing mesh
```

**Integration**: Simpson's rule
```python
integral = simpson(bs[:min_length] * wf1[:min_length] * wf2[:min_length], 
                   dx=self.boundState.dr)
product = integral * pf1 * pf2 * self.massFactor  # Include phase factors
```

The radial integrals are stored in a 4D array:
```python
self.radialInt[L1][index_J1][index_L2][index_J2]
```

### 5. Angular Distribution Calculation

**Function**: `AngDist(theta_deg)` ([dwba_zr.py:503-513](file:///home/ryan/PtolemyGUI/Raphael/dwba_zr.py#L503-L513))

**Algorithm**:
```python
def AngDist(self, theta_deg) -> float:
  self.PreCalLegendreP(theta_deg)  # Pre-calculate P_L^m(cos θ)
  xsec = 0
  for ma in [-s_a, ..., +s_a]:
    for mb in [-s_b, ..., +s_b]:
      for m in [range based on ma, mb, j]:
        beta = self.Beta(m, ma, mb)  # Compute amplitude
        xsec += |beta|²               # Sum incoherent spins
  
  return xsec * self.xsecScalingfactor * 10  # Convert fm² to mb
```

**Sum structure**:
1. Sum over initial spin projections `m_a` (unpolarized beam averaging)
2. Sum over final spin projections `m_b` (unobserved)
3. Sum over magnetic quantum numbers `m` (angular dependence)
4. Square amplitude and sum incoherently

**Conversion factor**: `10 fm² = 1 mb` (millibarn)

---

## Optimization Strategies

The code employs several pre-calculation strategies to avoid redundant computation:

### 1. Pre-calculate Clebsch-Gordan Coefficients

**Function**: `PreCalClebschGordan()` ([dwba_zr.py:375-419](file:///home/ryan/PtolemyGUI/Raphael/dwba_zr.py#L375-L419))

**Storage**: 6D array $\text{CG}[2j_1, 2m_1+\text{offset}, 2j_2, 2m_2+\text{offset}, 2j_3, 2m_3+\text{offset}]$
- Uses integer indexing by storing $2j$ and $2m$ (handles half-integers)
- Pre-computes all needed CG coefficients once during initialization
- Lookup in $O(1)$ time during angular distribution calculation

**Typical size**: ~1000-10000 coefficients depending on max $L$

### 2. Pre-calculate Wigner 9j-Symbols

**Function**: `PreCalNineJ()` ([dwba_zr.py:426-438](file:///home/ryan/PtolemyGUI/Raphael/dwba_zr.py#L426-L438))

**Storage**: 4D array indexed by $[L_1, \text{index}_{J_1}, \text{index}_{L_2}, \text{index}_{J_2}]$

The 9j-symbol structure is:
$$\begin{Bmatrix} j & l & s \\ J_1 & L_1 & s_a \\ J_2 & L_2 & s_b \end{Bmatrix}$$

**Calculation**: Uses `sympy.physics.quantum.cg.wigner_9j` (symbolic → numerical)

**Benefit**: 9j-symbols are computationally expensive; pre-calculation saves ~90% of runtime

### 3. Pre-calculate Associated Legendre Polynomials

**Function**: `PreCalLegendreP(theta_deg)` ([dwba_zr.py:496-501](file:///home/ryan/PtolemyGUI/Raphael/dwba_zr.py#L496-L501))

For each angle $\theta$, pre-compute:
$$P_L^m(\cos \theta) \quad \text{for all} \quad L \in [0, \text{maxL2}], \; m \in [0, L]$$

**Storage**: 2D array $\text{legendrePArray}[L][m]$ (only positive $m$; use symmetry for negative)

**Module**: Uses custom `associated_legendre_array()` from [assLegendreP.py](file:///home/ryan/PtolemyGUI/Raphael/assLegendreP.py)

### Performance Summary

From README: 
- **s-orbital transfer**: ~10 seconds
- **d-orbital transfer**: ~30 seconds

Performance breakdown:
- Bound state solution: ~1-5 sec
- Distorted wave calculation: ~5-10 sec  
- Radial integrals: ~1-5 sec
- Angular distribution (180 angles): ~10-20 sec
  - Pre-calculation overhead: ~0.1-0.5 sec
  - Per-angle calculation: ~50-100 ms

---

## Code Structure

### Main Class Hierarchy

```
SolvingSE (solveSE.py)
  ├─ Solves radial Schrödinger equation by RK4
  ├─ Manages potentials (Woods-Saxon, Coulomb, Spin-Orbit)
  ├─ Handles nuclear data (IAEA masses, spins)
  │
  ├── BoundState (boundState.py)
  │     └─ Finds bound state wave functions
  │
  └── DistortedWave (distortedWave.py)
        └─ Calculates scattering wave functions and S-matrix

DWBA_ZR (dwba_zr.py)
  ├─ Uses ReactionData for kinematics
  ├─ Creates BoundState for transferred nucleon
  ├─ Creates two DistortedWave objects (incoming/outgoing)
  ├─ Calculates radial integrals
  ├─ Pre-computes angular coefficients
  └─ Computes differential cross section
```

### Key Data Structures

**Radial Integrals** (4D complex array):
```python
radialInt[L1][index_J1][index_L2][index_J2]
# Dimensions: [maxL1+1][2*s_a+1][2*l+1][2*s_b+1]
```

**Clebsch-Gordan** (6D float array):
```python
CG[2*j1][2*m1+offset][2*j2][2*m2+offset][2*j3][2*m3+offset]
```

**Wigner 9j** (4D float array):
```python
NineJ[L1][index_J1][index_L2][index_J2]
```

**Distorted Waves** (nested lists):
```python
wfu_I[L][index_J]  # List of complex arrays (r-dependent)
wfu_O[L][index_J]
```

### Key Methods Flow

```
DWBA_ZR.__init__()
  ├─ Parse reaction via ReactionData
  ├─ Setup bound state (BoundState)
  ├─ Setup incoming channel (DistortedWave + optical potential)
  ├─ Setup outgoing channel (DistortedWave + optical potential)
  ├─ PreCalNineJ()
  └─ PreCalClebschGordan()

FindBoundState()
  └─ BoundState.FindPotentialDepth()

CalScatMatrixAndRadialIntegral()
  ├─ dwI.CalScatteringMatrix() → (S_matrix, distorted_waves)
  ├─ dwO.CalScatteringMatrix() → (S_matrix, distorted_waves)
  ├─ Rescale radial meshes by mass factor
  ├─ Simpson integration: ∫ bound × incoming × outgoing dr
  └─ Store in radialInt[L1][J1][L2][J2]

CalAngDistribution(angMin, angMax, angStep)
  └─ for each angle:
       └─ AngDist(theta)
            ├─ PreCalLegendreP(theta)
            └─ ∑_{ma,mb,m} |Beta(m,ma,mb)|²

Beta(m, ma, mb)
  └─ ∑_{L1,J1,L2,J2} Gamma(...) × Legendre × RadialInt

Gamma(L1, J1, L2, J2, m, ma, mb)
  └─ 9j-symbol × phase factors × 4 CG coefficients
```

---

## Selection Rules and Constraints

### Implemented Selection Rules

1. **Parity Conservation**: $L_1 + L_2 + l$ must be even
   ```python
   if int(L1 + L2 + self.l) % 2 != 0: return 0
   ```

2. **Triangle Inequalities**:
   - $|J_A - J_B| \leq j \leq J_A + J_B$ (total angular momentum)
   - $|s_a - s_b| \leq s \leq s_a + s_b$ (channel spin)
   - $|l - s| \leq j \leq l + s$ (orbital coupling)
   - $|J_1 - j| \leq J_2 \leq J_1 + j$ (channel coupling)
   ```python
   if not obeys_triangle_rule(J1, self.j, J2): continue
   ```

3. **Magnetic Quantum Number Constraint**: $|m| \leq L_2$
   ```python
   if not(abs(m) <= L2): continue
   ```

4. **Valid Spin Projections**:
   ```python
   if not(abs(m-mb+ma) <= self.j): continue
   if not(abs(mb-m) <= J2): continue
   ```

### Range of Summations

The code automatically determines summation ranges:

**Incoming channel**: $L_1 \in [0, \text{maxL1}]$ where maxL1 depends on convergence (typically 10-20)

**Outgoing channel**: $L_2 \in [L_1-l, L_1+l]$ (limited by orbital angular momentum transfer)

**Total angular momenta**: 
- $J_1 \in [|L_1 - s_a|, L_1 + s_a]$
- $J_2 \in [|L_2 - s_b|, L_2 + s_b]$

**Magnetic quantum numbers**: Determined by $m \in [m_b - m_a - j, m_b - m_a + j]$

---

## Physical Constants

From [dwba_zr.py](file:///home/ryan/PtolemyGUI/Raphael/dwba_zr.py):

- $D_0 = 1.55 \times 10^4$ MeV·fm³ (Zero-range normalization constant)
- $\hbar c = 197.326979$ MeV·fm (Reduced Planck constant $\times$ speed of light)
- $m_n = 939.56539$ MeV (Neutron mass)
- $\text{amu} = 931.494$ MeV (Atomic mass unit)
- $e^2 = 1.43996$ MeV·fm (Coulomb constant, $e^2/4\pi\epsilon_0$)

**Unit Conversions**:
- Cross section: $\text{fm}^2 \to \text{mb}$ (multiply by 10)
- Energy: Laboratory $\to$ Center-of-mass frame
- Wave number: $k = \sqrt{2\mu E}/\hbar$

---

## References

Theory is detailed in the thesis/handbook referenced in [README.md](file:///home/ryan/PtolemyGUI/Raphael/README.md):
- [Handbook of Direct Nuclear Reaction for Retarded Theorist v1](https://nukephysik101.wordpress.com/2025/02/21/handbook-of-direct-nuclear-reaction-for-retarded-theorist-v1/)

### Recommended Reading

1. **Satchler, "Direct Nuclear Reactions"** - Classic DWBA reference
2. **Thompson, "Coupled Reaction Channels Theory"** - Advanced formalism
3. **Austern, "Direct Nuclear Reaction Theories"** - Mathematical foundations
4. **Koning & Delaroche (2003)** - Optical potential parametrization
5. **An & Cai (2006)** - Deuteron global optical potential

---

## Limitations and Future Extensions

### Current Limitations (from README)
1. **Zero-range only**: Finite-range effects not included (affects absolute normalization)
2. **Single nucleon transfer**: Two-nucleon transfers not supported
3. **No inelastic scattering**: Only working for transfer reactions
4. **No coupled-channels**: Treats channels independently
5. **No polarization observables**: Only unpolarized cross sections

### Planned Extensions
1. Finite-range DWBA calculation
2. Inelastic scattering support
3. Polarization observables (analyzing powers, spin transfer)
4. Two-nucleon transfer reactions  
5. Coupled-channels formalism
6. C++ implementation or multiprocessing for speed

---

## Conclusion

This implementation provides a modern, readable Python implementation of zero-range DWBA with:
- **Modular design**: Clear separation of bound states, distorted waves, and angular coupling
- **Optimization**: Pre-calculation of expensive angular momentum coefficients
- **Flexibility**: Easy to modify optical potentials and bound state parameters
- **Transparency**: Explicit implementation of all mathematical steps

The code prioritizes clarity and maintainability over raw performance, making it suitable for:
- Educational purposes (understanding DWBA formalism)
- Prototyping new reaction mechanisms
- Cross-checking legacy Fortran codes
- Rapid analysis of experimental data

Future C++ implementation or parallelization would enable production-level throughput while maintaining this clear mathematical structure.
