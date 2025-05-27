# Simple Rebound

A simplified N-body simulation package based on the REBOUND framework.

## Installation

To install in development mode:

```bash
pip install -e .
```

## Usage

```python
import simple_rebound

# Create a simulation
sim = simple_rebound.Simulation()

# Add particles and run simulation
sim.add(...)
sim.integrator = 'ias15' # leapfrog, whfast, mercurius
sim.dt = ...
sim.integrate(t_max = ...)
```

## Requirements

- Python >= 3.6
- NumPy

## License

MIT License
