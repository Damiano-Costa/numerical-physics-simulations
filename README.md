# numerical-physics-simulations
Numerical simulations in physics: ODE integrators, nonlinear dynamics, stochastic processes and Monte Carlo methods implemented in C .

## System Requirements
To build and run the simulation modules, the following environment is required:
* **C Compiler:** `gcc` or `clang`
* **Math Library:** Standard C math library (`libm`)
* **Python Environment:** Python 3.8+ (Dependencies: `numpy`, `matplotlib`)

## Scripts

### 1. Dynamical Systems Integration Core 
numerical integration of ordinary differential equations (ODEs), focusing on stability and phase-space analysis of non-linear systems.
* **`RungeKutta.c`**: Full differential equation solver for damped/driven pendulums. Includes energy convergence analysis.
* **`Eul_EulCro.c`**: Harmonic oscillator simulation engine to benchmark integration algorithms (Euler vs. Verlet vs. Midpoint).
* **`Poincare.c`**: RK4-based integrator optimized for generating Poincaré sections and continuous phase-space trajectories.
* **`Biforcazioni.c`**: Generates bifurcation diagrams to analyze system transitions from deterministic to chaotic states.
* **`BaciniAttr.c`**: Grid-based computing module to classify basins of attraction across different initial conditions.

### 2. Stochastic Simulation Engine & PRNG Benchmarking 
Tools for generating stochastic processes and evaluating the quality of pseudo-random algorithms.
* **`Gen.c`**: Benchmarking suite for evaluating 4 distinct PRNGs through 1D Random Walk simulations.
* **`RW1d.c`**: Configurable 1D Random Walk simulator. Supports macro-driven execution to output single trajectories, MSD, or probability distributions.
* **`BoxMuller.c`**: Implementation of the Box-Muller transform to efficiently generate normally distributed random numbers.

### 3. Lattice Gas Dynamics 
* **`Latgas.c`**: 2D Lattice Gas simulation utilizing Random Walk dynamics with excluded volume. Highly configurable via compile-time directives.
* **`staticfluct.py`**: Python data pipeline script to compute and plot statistical fluctuations and standard deviations.

## Build Instructions
The C modules utilize standard libraries but require the math library flag (`-lm`). Some modules also rely on preprocessor directives to switch runtime modes.

**Standard Compilation:**

```bash
gcc -O3 RungeKutta.c -o RungeKutta -lm
