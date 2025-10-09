# Rust-Physics-Simulations-N-Body-Fluid-Flow-and-Kuramoto-Sivashinsky-Dynamics
High-performance physics simulations implemented in Rust. Includes an N-Body gravitational system using Barnes–Hut approximation, a 2D incompressible fluid flow simulation generating a Kármán vortex street, and a chaotic Kuramoto–Sivashinsky equation solver using the ETDRK4 scheme with real-time and static visualization.
# Rust Physics Simulations  
### Numerical Models for N-Body Systems, Incompressible Fluid Flow, and the Kuramoto–Sivashinsky Equation

This repository contains three independent numerical solvers written in **Rust**, each implementing a distinct class of physical dynamics.  
The projects demonstrate computational modeling using efficient and stable numerical methods, combined with modern visualization frameworks.  

All simulations are implemented in two spatial dimensions and emphasize clarity, numerical stability, and reproducibility.

---

## 1. N-Body Gravitational Simulation

This program models the gravitational interaction between a collection of point masses under Newtonian gravity.  
The motion of each particle is governed by the coupled second-order differential equations:

$$
\ddot{\mathbf{r}}_i = G \sum_{j \ne i} \frac{m_j (\mathbf{r}_j - \mathbf{r}_i)}{|\mathbf{r}_j - \mathbf{r}_i|^3},
$$

where \( G \) is the gravitational constant and \( \mathbf{r}_i \) denotes the position of the \( i \)-th body.  

To achieve computational efficiency, the simulation employs the **Barnes–Hut tree algorithm**, which approximates distant clusters of bodies as single composite nodes.  
This reduces the computational scaling from \(O(N^2)\) to \(O(N \log N)\).  

The integration of the equations of motion uses a **symplectic leapfrog scheme**, which conserves total energy over long time spans.  
A real-time visualization displays the evolving particle trajectories and cluster formation.

**Run:**
```bash
cd N_Body_Simulation
cargo run --release
