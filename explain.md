# Argon Molecular Dynamics: What This Code Is Doing

This project is a small classical molecular dynamics simulator for argon. It moves a set of atoms through time by solving Newton's equations numerically. The code is organized into a few Fortran source files, each with a single job: read the input, initialize coordinates and velocities, compute forces from a Lennard-Jones potential, advance the system in time, and compute the instantaneous temperature.

The rest of this document explains the program in two parallel ways:

1. In physics terms: what the simulation is modeling.
2. In code terms: how the Fortran files work together.

## 1. The big picture

At a high level, the program simulates a box of argon atoms.

Physics view:
- Each atom feels forces from every other atom.
- Those forces come from the Lennard-Jones potential, which is a simple model of how neutral atoms attract and repel each other.
- The atoms move according to Newton's second law, $F = ma$.
- Periodic boundary conditions are used, so the box behaves like it repeats forever in all directions.
- The simulation starts with a thermostat-like velocity scaling step so the atoms begin near a target temperature.
- During the equilibration period, the code keeps rescaling velocities so the temperature stays near the target value.
- After equilibration, the thermostat is turned off and the system evolves more naturally.

Code view:
- `main.f90` reads input, allocates arrays, calls initialization, computes the first forces, and then runs the time-stepping loop.
- `init.f90` reads coordinates from an XYZ file and generates initial velocities.
- `force.f90` computes forces and potential energy.
- `integrate.f90` advances positions and velocities with a Velocity Verlet scheme.
- `temp.f90` computes the current temperature from the velocities.

## 2. The physics model

### 2.1 Why argon is a good test case

Argon is a monatomic noble gas. That makes it much simpler than a molecule like water or methane because there are no bonds, angles, or internal vibrations to track. The atoms are treated as point particles with a mass and a position. That is exactly the kind of system where classical molecular dynamics is easiest to understand.

### 2.2 The Lennard-Jones potential

The interaction model used here is the Lennard-Jones potential:

$$
V(r) = 4\epsilon \left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^6\right]
$$

This potential has two parts:
- The $r^{-12}$ term is strongly repulsive at short distances. It prevents atoms from collapsing into each other.
- The $r^{-6}$ term is attractive at longer distances. It models weak van der Waals attraction.

The parameters mean:
- $\sigma$ is the distance scale. It is roughly the size of the atom.
- $\epsilon$ is the energy scale. It tells you how deep the attractive well is.

The force is the negative derivative of the potential, so when the atoms are too close the force pushes them apart, and when they are moderately separated the force pulls them together.

### 2.3 Periodic boundary conditions

The atoms are not simulated in open space. Instead, the code puts them in a cubic box and uses periodic boundary conditions. If an atom leaves one side of the box, it re-enters from the opposite side.

Physics meaning:
- The box is treated as a tiny piece of a much larger bulk material.
- This avoids wall effects and makes the system better represent a liquid or dense gas.

Code meaning:
- After positions are updated, any coordinate smaller than zero gets wrapped back by adding the box length.
- Any coordinate larger than the box gets wrapped back by subtracting the box length.
- In the force calculation, the minimum image convention is used, meaning the code always measures the shortest vector between two atoms through the periodic box.

### 2.4 The force cutoff

The code does not compute Lennard-Jones interactions out to infinite distance. Instead, it uses a cutoff radius `Rcut`.

Why this is done:
- The Lennard-Jones force becomes very small far away.
- Ignoring very distant interactions saves time.
- The cost of the force loop is otherwise proportional to $N^2$.

To reduce the jump in potential energy at the cutoff, the code shifts the potential by subtracting the value at `Rcut`. That makes the potential zero at the cutoff radius.

Important detail:
- The energy is shifted, but the force is still abruptly set to zero beyond the cutoff. That is common in simple MD codes, although it is not perfectly smooth.

### 2.5 Time integration

The program uses the Velocity Verlet method.

This matters because Newton's equations are continuous in time, but a computer can only update the atoms in discrete steps. So the code approximates the motion over small time intervals `TimeStep`.

Velocity Verlet works by:
- moving the positions using the current velocities and forces,
- recalculating the forces at the new positions,
- then finishing the velocity update.

This method is popular because it is stable, reversible enough for many MD uses, and conserves energy reasonably well when the thermostat is off and the timestep is small enough.

### 2.6 Temperature and the thermostat

The instantaneous temperature is estimated from the kinetic energy using equipartition:

$$
T = \frac{\sum_i m v_i^2}{3N-3}
$$

The factor $3N-3$ is used because the code removes the center-of-mass momentum, so three degrees of freedom are effectively taken away.

During equilibration, the code rescales velocities to keep the system near the target temperature. This is a simple thermostat. It is useful for preparing the system, but it is not a fully physical thermostat for production sampling because it directly changes the velocities.

After equilibration, that scaling stops and the atoms evolve with their natural dynamics.

## 3. File-by-file explanation

### 3.1 `main.f90`: the driver

This is the top-level program. It is the conductor that tells the rest of the code when to do each task.

What it reads from the input file:
- number of molecules
- atoms per molecule
- timestep
- total number of MD steps
- equilibration length
- target temperature
- coordinate file name
- box size
- Lennard-Jones `sigma`
- Lennard-Jones `epsilon`
- cutoff radius
- atomic mass

What it does with those values:
- writes a summary to `md.out`
- converts the input values into atomic units
- computes the total number of atoms as `NMol * NAtom`
- allocates arrays for positions, velocities, forces, and labels
- calls `initialize` to set up the starting state
- calls `force_calc` once to get the initial forces and initial potential energy
- loops over MD steps
- calls `integrate` once per step
- writes potential, kinetic, and total energies to `md.out`

Physics meaning:
- The system starts from a saved coordinate file.
- Initial velocities are assigned and rescaled.
- Then the simulation proceeds in small increments through time.

Code meaning:
- `main` manages input/output and the outer loop.
- Almost all physics work happens in the subroutines.

One important practical detail is unit conversion. The input file appears to use ordinary simulation units such as femtoseconds, angstroms, Kelvin, amu, and kcal/mol. The code converts everything to atomic units before the actual dynamics begins.

### 3.2 `init.f90`: reading coordinates and creating initial velocities

This subroutine prepares the starting configuration.

Step by step:
1. It opens the XYZ coordinate file.
2. It checks that the number of atoms in the file matches the expected total.
3. It reads each atom label and its coordinates.
4. It converts coordinates from angstroms to atomic units.
5. It applies a simple wrap into the periodic box if any coordinate lies outside the box range.
6. It assigns random initial velocities.
7. It subtracts the average velocity so the total linear momentum is zero.
8. It rescales the velocities so the system starts near the requested temperature.
9. It prints the resulting temperature estimate to `md.out`.

Physics meaning:
- The initial atomic positions come from an input structure file.
- Random velocities imitate a thermalized starting state.
- Zeroing the center-of-mass velocity prevents the whole box from drifting through space.
- Rescaling to the target temperature gives a reasonable starting kinetic energy.

Code meaning:
- The routine directly fills the position array `r` and the velocity array `v`.
- It uses the mass and target temperature to compute a scale factor.
- It uses the velocity vector average to remove bulk motion.

One thing to notice is that the random velocities are generated with `rand()`, so the exact initial state can vary from run to run.

### 3.3 `force.f90`: pair forces and potential energy

This is the core physics kernel.

What it computes:
- the total force on each atom
- the total potential energy of the system

How it works:
1. It computes the squared cutoff distance.
2. It computes the Lennard-Jones potential value at the cutoff, so the final energy can be shifted to zero there.
3. It loops over all distinct atom pairs `i < j`.
4. For each pair, it forms the displacement vector between the atoms.
5. It applies the minimum image convention with the periodic box.
6. It computes the squared separation.
7. If the pair is inside the cutoff, it evaluates the Lennard-Jones force and potential.
8. It adds equal and opposite force contributions to the two atoms, so Newton's third law is satisfied.
9. It accumulates the shifted potential energy.

Physics meaning:
- Every atom interacts with every other atom within the cutoff.
- The force on one atom is exactly the opposite of the force on the partner atom.
- The potential is a sum over pairs, which is the standard classical approximation for simple noble-gas liquids.

Code meaning:
- The double loop makes the cost scale like $N^2$.
- The routine avoids taking square roots by working mostly with squared distances.
- The force array is reset to zero before accumulation.
- `atom1` and `atom2` are used as loop indices from the shared `general` module.

The minimum image step is especially important. Without it, atoms near opposite sides of the box would appear far apart even though they are actually neighbors through the periodic boundary.

### 3.4 `integrate.f90`: moving the system forward

This routine advances the atom positions and velocities by one time step.

The intended physical sequence is:
1. Update positions using the current velocities and forces.
2. During equilibration, rescale velocities toward the target temperature.
3. Recompute forces at the new positions.
4. Finish the velocity update using the new forces.
5. Wrap any coordinates back into the periodic box.
6. Compute kinetic energy and the current temperature estimate.

Physics meaning:
- This is the discrete approximation to the continuous motion governed by $F = ma$.
- The thermostat stage keeps the system from drifting far from the requested temperature while it is equilibrating.
- Once equilibration is over, the motion is purely Newtonian.

Code meaning:
- `Step` is computed from `t / TimeStep` to decide whether equilibration is still active.
- The position update uses the current force divided by the mass.
- The subroutine calls `force_calc` again after moving the atoms.
- The velocity update uses the new forces.
- The kinetic energy is computed from the squared velocities.

There is a subtle implementation detail worth noticing: the thermostat temperature check is performed inside the loop over atoms. That means the scaling factor is recomputed repeatedly during a single integration step, using the current state of the velocity array at that moment. The code still follows the intended idea of velocity rescaling, but this is a detail that matters when reading the implementation closely.

### 3.5 `temp.f90`: instantaneous temperature

This file contains a small helper routine that computes temperature from the velocity array.

What it does:
- sums $v_x^2 + v_y^2 + v_z^2$ over all atoms
- multiplies by the mass
- divides by $3N - 3$ to get a temperature estimate

Physics meaning:
- In classical mechanics, temperature is proportional to average kinetic energy.

Code meaning:
- The routine is simple and isolated so `integrate.f90` can ask for a temperature estimate without duplicating the formula.

## 4. The main loop in plain language

Here is the whole simulation in simple steps:

1. Read the input file.
2. Convert all quantities to atomic units.
3. Read the initial coordinates from the XYZ file.
4. Generate random starting velocities.
5. Remove center-of-mass drift.
6. Scale velocities to the target temperature.
7. Compute the initial forces and potential energy.
8. Repeat for each MD step:
   - move the atoms forward a little bit,
   - optionally apply the thermostat while equilibrating,
   - recompute the forces,
   - finish the velocity update,
   - wrap atoms back into the periodic box,
   - compute kinetic energy,
   - write energies to the output file.

## 5. Important variables and what they mean

- `r(TotAtom, 3)`: position of each atom in three dimensions.
- `v(TotAtom, 3)`: velocity of each atom.
- `Force(TotAtom, 3)`: total force acting on each atom.
- `PE`: total potential energy.
- `KE`: total kinetic energy.
- `Box`: side length of the cubic periodic box.
- `Sig`: Lennard-Jones $\sigma$ parameter.
- `Eps`: Lennard-Jones $\epsilon$ parameter.
- `Rcut`: cutoff radius.
- `Mass`: atomic mass.
- `Temp`: target temperature.
- `TimeStep`: integration step size.
- `EQMDStep`: number of steps where the thermostat remains active.

## 6. What the output means

The file `md.out` begins with a summary of the simulation settings. After that, each step prints:
- the step number
- the potential energy
- the kinetic energy
- the total energy

How to interpret this:
- During equilibration, total energy is not expected to be constant because the thermostat changes the velocities.
- After equilibration, the total energy should be much more stable if the timestep is reasonable.
- If the total energy drifts badly, the timestep may be too large, the starting configuration may be too dense, or the cutoff/integration settings may be too aggressive.

## 7. What assumptions the code makes

This is a deliberately simple molecular dynamics code. That simplicity is useful, but it also means the model has limits.

- The atoms are monatomic.
- The interaction is pairwise Lennard-Jones only.
- The box is cubic.
- Long-range interactions are ignored beyond the cutoff.
- The thermostat is a simple velocity rescaling method, mainly for equilibration.
- The force evaluation is a direct all-pairs calculation, so performance scales poorly for large systems.

## 8. Short conceptual summary

If you want the shortest possible description, it is this:

The code starts with a box of argon atoms, gives them thermal velocities, computes pairwise Lennard-Jones forces, and then repeatedly moves the atoms forward in time using Velocity Verlet. While the system is equilibrating, it rescales velocities to hold the temperature near a target value. After that, it lets the system evolve naturally and records energies along the way.

## 9. Where to look in the source

- [main.f90](main.f90)
- [init.f90](init.f90)
- [force.f90](force.f90)
- [integrate.f90](integrate.f90)
- [temp.f90](temp.f90)
- [Readme.md](Readme.md)
- [Makefile](Makefile)
