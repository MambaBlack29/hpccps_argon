# Argon MD Report

## Code Structure Change

The current codebase introduces a dedicated Verlet neighbor-list layer in [verlet.f90](verlet.f90), which changes the flow of the molecular dynamics program from direct pairwise force evaluation to a neighbor-list-driven design.

In the previous structure at commit 88c7dc3fcda37f8c1815edd33976efbd2abf3cb4, the simulation was organized more simply:

- [main.f90](main.f90) handled the main simulation loop.
- [integrate.f90](integrate.f90) advanced positions and velocities.
- [force.f90](force.f90) computed Lennard-Jones forces with a full pairwise scan over all atoms.
- There was no separate verlet.f90 file or reusable neighbor-list stage.

In the current structure, the force calculation is separated from the neighbor search:

- [init.f90](init.f90) initializes the atom positions and velocities.
- [verlet.f90](verlet.f90) builds the neighbor list with new_verlet, using the minimum-image convention with periodic boundaries.
- [main.f90](main.f90) allocates and refreshes the neighbor-list arrays and rebuilds the Verlet list every 15 steps.
- [integrate.f90](integrate.f90) advances the system state and then calls the force update using the current neighbor list.
- [force.f90](force.f90) evaluates the Lennard-Jones interaction only for neighbors stored in the list.

This restructuring reduces repeated all-pairs work and makes the neighbor-list logic explicit and reusable. The fixed-size list buffers and periodic rebuild cadence also make the data flow easier to follow during the simulation.

## Results From plot_results.py

The plotting script compares serial wall time for three Verlet update intervals. The table below reproduces the values used in the plot. The baseline without Verlet is defined in the script as well, but that line is currently commented out in the figure generation.

| Method | 1600 points | 5400 points | 12800 points |
| --- | ---: | ---: | ---: |
| Without Verlet | 7.98 | 77.56 | 427.49 |
| With Verlet, update every 10 steps | 2.37 | 10.75 | 50.64 |
| With Verlet, update every 15 steps | 1.79 | 8.34 | 34.88 |
| With Verlet, update every 20 steps | 1.60 | 6.92 | 28.22 |

## Figures

Figure 1 shows the full comparison, including the non-Verlet baseline. Figure 2 zooms in on the Verlet-enabled runs so the timing differences between update intervals are easier to see.

![Figure 1](Figure_1.png)

![Figure 2](Figure_2.png)