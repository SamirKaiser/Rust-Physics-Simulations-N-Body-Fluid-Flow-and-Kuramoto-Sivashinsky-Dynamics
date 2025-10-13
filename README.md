# Rust-Physics-Simulations-N-Body-Fluid-Flow-and-Kuramoto-Sivashinsky-Dynamics
High-performance physics simulations implemented in Rust. Includes an N-Body gravitational system using Barnes–Hut approximation, a 2D incompressible fluid flow simulation generating a Kármán vortex street, and a chaotic Kuramoto–Sivashinsky equation solver using the ETDRK4 scheme with real-time and static visualization.
# N-Body Gravitational Simulation (Barnes–Hut Approximation)

## Overview

This project simulates the gravitational dynamics of many interacting point masses under Newtonian gravity.  
Each particle experiences the combined gravitational force of all others, governed by:

$$
\ddot{\mathbf{r}}_i = G \sum_{j \ne i} \frac{m_j (\mathbf{r}_j - \mathbf{r}_i)}{|\mathbf{r}_j - \mathbf{r}_i|^3},
$$

where $$G$$ is the gravitational constant and  $$\mathbf{r}_i$$ denotes the position vector of the $$i$$-th body.

## Numerical Method

- **Force computation:** The **Barnes–Hut approximation** groups distant particles into hierarchical nodes, achieving $$O(N \log N)$$ scaling.  
- **Time integration:** The **leapfrog symplectic method** ensures stability and long-term energy conservation.  
- **Parallelism:** The `rayon` crate enables multi-threaded force evaluation.  
- **Visualization:** Real-time display using `eframe` and `egui`.

## Physical Behavior

The simulation illustrates:
- Gravitational clustering and orbital motion.
- Dynamic equilibrium between attraction and momentum.
- Emergent galactic-like patterns.


---

# Incompressible Fluid Flow — 2D Navier–Stokes Solver

## Overview

This simulation solves the two-dimensional incompressible **Navier–Stokes equations** in the vorticity–streamfunction formulation:

$$
\frac{\partial \omega}{\partial t} + \mathbf{u} \cdot \nabla \omega = \nu \nabla^2 \omega,
$$

$$
\nabla^2 \psi = -\omega, \qquad \mathbf{u} = (\partial_y \psi, -\partial_x \psi),
$$

where $$\omega$$ is the vorticity, $$\psi$$ the streamfunction, and $$\nu$$ the kinematic viscosity.

## Numerical Method

- **Advection:** Semi-Lagrangian scheme for unconditional stability.  
- **Diffusion:** Central finite differences for the Laplacian term.  
- **Pressure projection:** Solved using a Poisson equation for incompressibility.  
- **Boundary conditions:** Periodic edges and a no-slip condition on a circular obstacle.  
- **Visualization:** Vorticity field rendered via `egui`.

## Physical Behavior

A uniform inflow from left to right passes a circular obstacle, generating a **Kármán vortex street** — a sequence of alternating vortices in the cylinder wake.  
This phenomenon is a fundamental example of laminar instability and vortex shedding.


---

# Kuramoto–Sivashinsky (KS) Equation Simulation

## Overview

The **Kuramoto–Sivashinsky (KS) equation** models nonlinear wave interactions and chaotic behavior in dissipative systems such as thin film flows and flame fronts:

$$
\frac{\partial u}{\partial t} + u \frac{\partial u}{\partial x}
+ \frac{\partial^2 u}{\partial x^2} + \frac{\partial^4 u}{\partial x^4} = 0.
$$

This system exhibits **spatiotemporal chaos**, where small perturbations grow, interact, and decay in a self-organizing but unpredictable manner.

## Numerical Method

- **Spatial discretization:** Fourier spectral differentiation using FFTs.  
- **Time integration:** Exponential Time-Differencing Runge–Kutta 4 (ETDRK4) scheme for stiff systems.  
- **Domain:** Periodic boundaries on $$x \in [0, L]$$.  
- **Visualization:** Spatio-temporal color plot of $$u(x,t)$$ using `plotters`.

## Physical Behavior

The KS equation demonstrates:
- Pattern formation and decay in nonlinear PDEs.  
- Onset of chaotic wave dynamics in one-dimensional media.  
- Sensitive dependence on initial conditions and system size.




