Lennard-Jones Intermolecular Potential Simulation

- This repository contains a MATLAB implementation of a Lennard Jones intermolecular potential simulation.
- The program computes Lennard-Jones potential curves, initializes particles inside a 3D cube, performs random displacement moves, evaluates system energy, and outputs final particle positions.
- The model reflects the classic 12-6 Lennard-Jones formulation used to describe short-range repulsion and long-range attraction between neutral atoms and molecules.

Overview:

- Reads Lennard-Jones parameters (epsilon, sigma) for Argon (ArAr) and Methane (CH4CH4) from lj_param.dat.
- Computes Lennard-Jones potential curves Eij(r) over the range [2, 10A].
- Prompts the user for:
- N = number of particles
- L = cube side length
- nmove = number of displacement steps
- delta = maximum displacement per step
- Randomly assigns particle positions and species.
- Computes total Lennard-Jones energy using species-specific parameters and Lorentz Berthelot mixing rules for Ar CH4 interactions.
- Applies random displacements and recalculates energy after each move.
- Writes final positions and energies to final_positions2.txt.
- Produces a 3D scatter plot of final particle positions.

1. How to Run:
Open MATLAB

2. Add the project folder to the path:
addpath('path/to/repo')

3. Run the script: LJ_simulation

4. Enter the required parameters when prompted: 
Number of particles (N), Cube side length (L), Number of moves (nmove), Maximum displacement (delta)

5. Outputs include: 
2D Lennard-Jones plots for ArAr and CH4CH4, Command-window energy logs, final_positions2.txt containing species and coordinates, 3D scatter plot of final particle positions