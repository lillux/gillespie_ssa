# Gillespie Stochastic Simulation Algorithm (SSA)
 
Python implementation of the Gillespie stochastic simulation algorithm for chemical reaction networks, plus a small evolutionary simulation module. The library focuses on simulating reaction systems defined by state-change vectors and combinatorial propensity functions, with utilities to compute propensities and trajectories.

## Scope

This repository provides:

- **Core SSA simulation** for reaction networks with user-defined propensities.
- **Dynamic/state-switching SSA** to alternate between reaction regimes over time. **Currently supports 2 states**
- **Evolutionary simulations** for simple wild-type vs. standing-variation population dynamics.
- **Notebook tutorials** that demonstrate common models and plotting workflows.

## Main functionalities

### Gillespie SSA (static reactions)
Use `gillespie.gillespie.gillespie_ssa` for classic SSA with fixed reactions.

```python
from gillespie.gillespie import gillespie_ssa

# Example: two reactions with two species
reagent_quantity = [50, 0]
state_change_vectors = [
    [-1, 1],  # A -> B
    [1, -1],  # B -> A
]
combinatorics = [
    lambda a, b: 0.1 * a,
    lambda a, b: 0.05 * b,
]

sim = gillespie_ssa(
    reagent_quantity=reagent_quantity,
    state_change_vectors=state_change_vectors,
    combinatorics=combinatorics,
    iteration=1000,
)

# Results
trajectory = sim.molecular_species_history
timepoints = sim.timestep_list
```

### Dynamic SSA (state switching)
Use `gillespie.gillespie_dynamic.gillespie_dynamic` to alternate between reaction sets (e.g., environmental regimes).

```python
from gillespie.gillespie_dynamic import gillespie_dynamic

reagent_quantity = [100, 0]
state_change_vectors = {
    "state_a": [[-1, 1]],
    "state_b": [[1, -1]],
}
combinatorics = {
    "state_a": [[lambda a, b: 0.2 * a]],
    "state_b": [[lambda a, b: 0.1 * b]],
}

sim = gillespie_dynamic(
    reagent_quantity=reagent_quantity,
    state_change_vectors=state_change_vectors,
    combinatorics=combinatorics,
    max_time=50.0,
    stop_condition="time",
    oscillate=True,
    oscillation_interval={"state_a": 10.0, "state_b": 10.0},
    start_with="state_a",
)

trajectory = sim.molecular_species_history
timepoints = sim.timestep_list
```

### Evolutionary simulations
Use `gillespie.evolution.EvoSim` and `gillespie.evolution.MultiSim` for population evolution experiments.

```python
from gillespie.evolution import EvoSim, MultiSim

sim = EvoSim(wt=1000, sv=10, r=0.2, s=0.5, u=1e-4, max_gen=200)

params = {"wt": 1000, "sv": 10, "r": 0.2, "s": 0.5, "u": 1e-4, "max_gen": 200}
multisim = MultiSim(epochs=100, params=params, ncore=-1, replicable=True)
```

### Stochastic utilities
Use `gillespie.stochastic_backend` for direct access to propensity and time-step helpers:

- `calculate_propensity_funct`
- `calculate_tau`
- `calculate_mu`
 
## Installation
 
The `gillespie` package can be installed with:
```bash
pip install git+https://github.com/lillux/gillespie_ssa.git@main#egg=gillespie
```
 
## Tutorials

The following notebooks provide guided examples and plots:
 

- **Schlögl model SSA walkthrough:** `gillespie.ipynb`: [link](https://github.com/lillux/gillespie_ssa/blob/main/gillespie.ipynb)
- **Dynamic SSA walkthrough:** `gillespie_dynamic_work.ipynb` [link](https://github.com/lillux/gillespie_ssa/blob/main/gillespie_dynamic_work.ipynb)
- **Lotka–Volterra example:** `lotka_volterra.ipynb`: [link](https://github.com/lillux/gillespie_ssa/blob/main/lotka_volterra.ipynb)
 

## Project layout

```
gillespie/
  gillespie.py            # Static SSA implementation
  gillespie_dynamic.py    # Dynamic/state-switching SSA
  stochastic_backend.py   # Propensity and sampling utilities
  evolution.py            # Evolutionary simulations
```
