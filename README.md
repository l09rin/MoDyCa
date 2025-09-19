# MoDyCa
**Mo**lecular **Dy**namics and Monte **Ca**rlo algorithms for numerical simulations of Soft Matter systems.
This code has been used to produce part of the simulation data published on _Soft matter 15 (40), 8113-8128 (2019)_ and _J. Chem. Phys. 160, 224109 (2024)_.

Two additional utilities are included:
- `histogram.exe` : to compute distributions from simulation outputs.
- `calculate_GofR.exe` : to efficiently compute radial pair distribution functions **(multi-thread parallelized)**.

---

## Getting started

### Clone & Compile
```bash
git clone git@github.com:l09rin/MoDyCa.git
cd MoDyCa
make all
```
This builds the main executable (`modyca.exe`) and the two utility programs.

### Install (optional)
Move the executables into your `PATH`:
```bash
mv modyca.exe histogram.exe calculate_GofR.exe  ~/.local/bin/
```

---

## Usage

The main program accepts two optional arguments:
1. `-in <file>` : Path to input parameter file.
2. `-seed <int>` : Random seed for initialization.

Example:
```bash
./modyca.exe -in input.dat -seed 84488
```

This will create output files containing results and simulation parameters in the **output directory**.  

👉 Example inputs are provided in the subfolder [examples/](https://github.com/l09rin/MoDyCa/tree/main/examples).

---

## Input
The simulation is controlled by a plain-text input file (examples in [examples/](https://github.com/l09rin/MoDyCa/tree/main/examples)).
Each line consists of a **keyword** followed by one or more parameters.
Comments start with `#`.

### General Settings
- **`DIMENSION <2|3>`**  :  Set simulation dimensionality (default = 3).
- **`PARTICLE_CLASS <particle_3D|particle_2D|patchy_2D <npatches>>`**  :  Choose particle class, storing different per-particle attributes (default = `particle_3D`).
- **`ENSEMBLE <type> <parameters>`**  :  Select simulation ensemble:
  - `NVE` (default, constant energy MD)
  - `NVT_NHC` (MD with Nosè–Hoover chains thermostat)
  - `BROWNIAN <friction>` (Brownian dynamics)
  - `LANGEVIN <friction>` (Langevin dynamics)
  - `MC_NVT <MC_step_length>` (Monte Carlo, translational/rotational moves)
  - `MC_NVT_CMfixed <MC_step_length>` (MC, constrained center of mass)
  - `MC_NPT <MC_step_length>` (MC, constant pressure)

### Monte Carlo Parameters
- **`MC_PARAMETER translational_probability <float>`**  :  probability to attempt a particle translation over a rotation.
- **`MC_PARAMETER max_translation <value> [AUTO_ADJUST <min> <max>]`**  :  maximum extent per each component of the random translation.
- **`MC_PARAMETER max_rotation <value> [AUTO_ADJUST <min> <max>] [JUMP <sym_factor>]`**  :  maximum value for randomly chosen rotation angle (`<sym_factor>` is the rotational symmetry factor of the particle, for example `4` for a tetra-valent patchy particle, used to attempt jump moves too).
- **`MC_PARAMETER max_boxside_change <value> [AUTO_ADJUST <min> <max>]`**  :  minimum and maximum linear variation in box side move attempts, in units of `SIGMA`.

### Physical Units
- **`SIGMA <float>`**  :  Length unit (default = 1, not explicitly used in code).  
- **`EPSILON <float>`**  :  Energy unit (default = 1, not explicitly used).  
- **`PARTICLE_MASS <float>`**  :  Particle mass (default = 1).  
- **`SIDE <float>`**  :  Size of a cubic/square box in units of `SIGMA`.  
- **`SIDES <float> <float>`**  :  For anisotropic boxes: specify each dimension separately.  
- **`DENSITY <float>`**  :  Number density (alternative to explicit side length).  
- **`KbT <float>`**  :  Thermal energy (in `EPSILON`), fixed in NVT/NPT ensembles.  
- **`PRESSURE <float>`**  :  Pressure (in `EPSILON/SIGMA^3`), fixed in NPT ensemble.  

### Initial Configuration
1. **From file**  
   ```
   INITIAL_CONFIGURATION file <xyz|nico_bonds|nico_bonds_pbc|nico_rings|nico_rings_pbc> <file(s)> <energy_per_particle>
   ```
   - `xyz` : _XYZ_-like format
   - `sph` : sph data file as dumped by this [Event Driven Molecular Dynamics code](https://github.com/FSmallenburg/EDMD)
   - `lmp` : LAMMPS data file
   - `nico_bonds`, `nico_bonds_pbc` : format encoding molecules with bonded particles
   - `nico_rings`, `nico_rings_pbc` : format encoding ring polymers with polydispersity

2. **Generated**
   ```
   INITIAL_CONFIGURATION <random|fcc|bcc|hardspheres|lattice> <N> <energy_per_particle> [options]
   ```
   Options include:
   - `hardspheres <HS_radius>`  :  for non-overlapping particles
   - `lattice square|hex`  :  for 2D square or hexagonal lattices
   - `fcc|bcc`  :  for 3D face-centered or body-centered cubic lattices
   - `random`  :  for disordered configurations (energy is calculated according to a truncated and shifted Lennard-Jones)

- **`CHARGES  all|all_noXl|file [random <fraction_of_charges>] [surface|gauss <fraction_of_charges> <sigma_of_distrib_as_fraction_of_max_mol_radial_extent>]`**  :  to setup various charge distributions, to use before enabling Debye-Huckel interactions.

- **`INITIAL_VELOCITIES_BY_FILE <flag> <file>`** : If `flag=1`, load velocities from `file`.

### Interaction Potentials
Several models are supported (parameters vary):
- **`LENNARD_JONES <cutoff> <shift|auto> <verlet_rad_delta>`**  :  defines Lennard–Jones potential
  - `cutoff` : cutoff distance in `SIGMA` units
  - `shift` : energy shift in `EPSILON` units (or `auto` for automatic)
  - `verlet_rad_delta` : Verlet list skin in `SIGMA` units
- **`WCA_ALPHA <cutoff> <solvophobic_parameter> <verlet_rad_delta>`**  :  defines repulsive WCA + solvophobic attractive potential
  - `solvophobic_parameter` : solvophobic parameter α as in the interaction model used in _Soft matter 15 (40), 8113-8128 (2019)_
- **`KERN_FRENKEL  <core_radius> <patch_radial_thickness> <deltatheta_patch> <bond_energy> <verlet_radius>`**  :  defines Kern-Frenkel model of patchy particle potential (only in 2D)
  - `core_radius` : hard core radius in `SIGMA` units
  - `patch_radial_thickness` : radial extent of the attractive patch from the hard core surface in `SIGMA` units
  - `deltatheta_patch` : angular half-opening of the attractive patch in radians units
  - `bond_energy` : patch-patch binding energy in `EPSILON` units
  - `verlet_radius` : cutoff radius of the Verlet list in `SIGMA` units
- **`HARD_SPHERE  <core_radius> <verlet_radius>`**  :  defines excluded volume hard sphere interactions
- **`HS_SQUARE_WELL <core_radius> <square_well_radius> <square_well_energy> <verlet_radius>`**
  - `square_well_radius` : radial extent of the attractive/repulsive well in `SIGMA` units
  - `square_well_energy` : square-well attractive/repulsive energy in `EPSILON` units
- **`HARMONIC_SPRING_FIELD <spring_constant>`**  :  defines harmonic spring attractions towards initial lattice sites.
  - `spring_constant` : spring elastic constant in `EPSILON/SIGMA^2` units
- **`HARMONIC_SPRING_FIELD_WIGNERSEITZ <spring_constant> <neigh_cutoff>`**  :  defines harmonic spring attractions towards initial lattice sites, with further constraint to preserve local connectivity.
  - `neigh_cutoff` : cutoff used to build connectivity network and respective Wigner-Seitz-like cells, in `<debye_length>` units
- **`ANGULAR_HARMONIC_SPRING_FIELD <spring_constant> <prefactor> <npatches>`**  :  defines harmonic angular spring attractions towards initial orientations (only in 2D)
  - `spring_constant` : interaction pre-factor in `EPSILON/SIGMA^2` units
  - `prefactor` : global prefactor, used in thermodynamic integration calculations
  - `npatches` : rotational symmetry factor of the particle
- **`FENE <force_constant> <max_extension>`**  :  Finitely Extensible Non-linear Elastic bonds
  - `force_constant` : interaction pre-factor in `EPSILON/SIGMA^2` units
  - `max_extension` : maximum bond extension in `SIGMA` units
- **`DEBYE_HUCKEL <bjerrum_length> <debye_length> <cutoff> <shift|auto> <verlet_rad_delta>`**  :  defines a screened Coulomb interaction
  - `cutoff` : cutoff distance specified in `<debye_length>` units

### MD Integration & Simulation Control
- **`INTEGRATION_TIME_STEP <float>`** : Time step (in `SIGMA * sqrt(MASS/EPSILON)`).
- **`EQUILIBRATION_STEPS <int>`** : Number of equilibration steps.
- **`PRODUCTION_STEPS <int>`** : Number of production steps.
- **`CM_RESET_INTERVAL <int>`** : Reset interval for center-of-mass (NVT).

### Output Control
- **`SAVING_INTERVALS <global_quantities> <configurations> <backup>`**  :  Output frequency for thermodynamic observables, snapshots, and position/velocity backup files.
- **`SAVE_FORMAT <xyz>`**  :  Format for particle trajectory (`xyz` supported).
- **`SAVE_ATTRIBUTE <name> <eq_flag> <prod_flag>`**  :  Control saving of attributes in equilibration / production runs (flags can be set to `0` or `1`). Options:
  - `particles`
  - `forces`
  - `velocities`
  - `energy`
  - `virial`
  - … (others listed in example file)
- **`SAVE_POTENTIAL_ENERGY <eq_flag> <prod_flag>`**
- **`SAVE_PARTIAL_POTENTIAL_ENERGIES <eq_flag> <prod_flag>`**
- **`SAVE_KINETIC_ENERGY <eq_flag> <prod_flag>`**
- **`SAVE_TEMPERATURE <eq_flag> <prod_flag>`**
- **`SAVE_VOLUME <eq_flag> <prod_flag>`**
- **`SAVE_MSD <eq_flag> <prod_flag> <format>`**
  - `format = all` (per-particle)
  - `format = mols_cm` (molecule centers of mass)
- **`SAVE_VIRIAL_PRESSURE <eq_flag> <prod_flag>`**

- **`SAVING_DIRECTORY_NAME <string>`** : Directory for outputs.

---

## Output
The simulation produces:
- **Trajectory files** (particle positions, velocities, etc.).
- **Thermodynamic logs** (energy, temperature, pressure, etc.).
- **Configuration backups**  (bond network, if present, and position and velocities with machine precision).
- **Runtime analysis outputs** (if enabled).

Outputs are placed in the directory specified by `SAVING_DIRECTORY_NAME`  (default: current working directory).

---

## Utilities
- `histogram.exe` : to compute distributions from simulation outputs.
- `calculate_GofR.exe` : to efficiently compute radial pair distribution functions **(multi-thread parallelized)**.

To get information on their usage:
```bash
./histogram.exe --help
./calculate_GofR.exe --help
```

---
