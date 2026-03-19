# Orbital Dynamics Simulator
A quick research project (and study) on the different numeric integrators used to simulate orbital dynamics. Something I want to do is to conclude the most accurate, the computationally faster, and the overall best one that can be implemented.

_The main aim of this research is going to be to figure out if developing a custom integrator offers tangible advantages over existing methods (mainly for my project [Fizix Mech](https://github.com/fyzl329/fizixmech))_

## Project Structure
```
_orbitalDynSim/
├── simulators/     # Integration method implementations (euler.py, verlet.py, rk4.py)
├── data/           # Simulation outputs (CSV logs, PNG plots)
├── analysis/       # Analysis and visualization scripts
├── README.md
└── .gitignore
```

## Usage

### Run a simulation
```bash
python simulators/euler.py
```
Outputs saved to `data/` folder.

### Analyze simulation data
```bash
python analysis/csv-plot.py
```
Reads CSV from `data/` and saves plots to `data/` folder.

## Core Comparison
- Euler Method
- Verlet
- Runge-Kutta 4 (RK4)

## Advanced Method (to explore)
- Runge-Kutta-Fehlberg 45 (RK45)
- Leapfrog

## Tech Stack

**Core:**
Python
numpy

**Plotting**
matplotlib

**Data / Logging**
pandas

**Version Control / Portfolio:**
GitHub

**Probably will add at the latter end of the research**
scipy
numba