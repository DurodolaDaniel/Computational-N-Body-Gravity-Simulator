# Computational N-Body-Gravity-Simulator
Computational N-Body Gravity Simulator

A Python-based computational physics project that simulates gravitational interactions in two-body, three-body, and N-body systems using numerical methods.

Features

* Two-body simulations: Sun–Earth and Sun–Mars
* Three-body simulations with stars and planets
* General N-body system with customizable/randomized initial conditions
* 2D and 3D trajectory visualization
* Animated gravitational motion
* Distance and position analysis
* Numerical integration of the equations of motion
* Gravitational force and acceleration calculations

Physics

The simulator is based on Newton’s law of universal gravitation:

[
F = G\frac{m_1m_2}{r^2}
]

The equations of motion are solved numerically to determine the position and velocity of each body as the system evolves with time.

Requirements

* Python
* NumPy
* Matplotlib
* SciPy

Install dependencies with:

pip install -r requirements.txt

Running the Project

The project is provided as a Jupyter Notebook (.ipynb) and is designed to be executed cell-by-cell.

For interactive plots and animations, include:

%matplotlib qt

with the library imports when running in Jupyter Notebook, VS Code notebooks, or Spyder.

Running cells individually prevents multiple animations from running simultaneously and improves performance.

Output

The simulator produces:

* Orbital trajectory plots
* 3D gravitational trajectories
* Earth–Sun and Earth–Sun–Mars visualizations
* Animated three-body systems
* Animated N-body gravitational systems

Author

Durodola Daniel
Physics Undergraduate | Computational Physics & Scientific Computing