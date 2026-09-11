# Computational N-Body Gravity Simulator

A Python-based computational physics project for simulating gravitational dynamics in interacting two-body, three-body, and N-body systems using Newtonian mechanics and numerical integration.

The simulator computes the gravitational acceleration of each body from the combined influence of all other bodies and evolves the system through discrete time steps. It supports orbital simulations, multi-body interactions, trajectory visualization, animation, and basic conservation diagnostics.

## Key Features

- Two-body gravitational simulations
- Three-body orbital dynamics
- General N-body gravitational interactions
- Newtonian gravitational force model
- Velocity Verlet numerical integration
- 2D trajectory visualization
- 3D trajectory visualization
- Animated orbital dynamics
- Custom initial conditions
- Multiple predefined physical systems
- Energy and angular-momentum diagnostics
- Center-of-mass tracking
- Randomized N-body configurations
- SI-unit based physical parameters

## Physical Model

The gravitational interaction between two bodies is described by Newton's law of universal gravitation:

\[
\mathbf{F}_{ij}
=
-G\frac{m_i m_j}{|\mathbf{r}_i-\mathbf{r}_j|^3}
(\mathbf{r}_i-\mathbf{r}_j)
\]

The acceleration of body \(i\) is obtained from the resultant gravitational influence of all other bodies:

\[
\mathbf{a}_i
=
-G\sum_{j\neq i}
m_j
\frac{\mathbf{r}_i-\mathbf{r}_j}
{|\mathbf{r}_i-\mathbf{r}_j|^3}
\]

The equations of motion are then integrated numerically to obtain the positions and velocities of the bodies as functions of time.

## Numerical Method

The simulator uses the **Velocity Verlet integration scheme**, which is well suited to classical dynamical systems and orbital simulations because of its favorable long-term energy behavior compared with simple Euler integration.

At each time step, the algorithm updates:

1. Position
2. Gravitational acceleration
3. Velocity

This allows the system to evolve dynamically rather than calculating an orbit from a predefined analytical path.

## Coordinate System

The simulator represents each body using Cartesian coordinates:

- \(x, y\) for 2D simulations
- \(x, y, z\) for 3D simulations

Initial conditions can be defined using radial distance and angular parameters before being converted to Cartesian coordinates.

## Simulations

The project can be used to investigate:

### Two-Body Systems
Examples include:

- Sun–Earth
- Binary-star systems

### Three-Body Systems
Examples include:

- Sun–Earth–Mars
- Binary stars with a smaller third body

### N-Body Systems
Multiple bodies interact simultaneously through their mutual gravitational fields, allowing the study of complex and chaotic dynamical behavior.

## Visualization

Simulation results can be examined through:

- Orbital trajectory plots
- 2D position plots
- 3D trajectory plots
- Animated N-body motion
- Energy evolution
- Angular momentum
- Center-of-mass motion

## Running the Project

### Requirements

Python 3.x

Install the required libraries:

```bash
pip install numpy matplotlib
