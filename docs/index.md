[![Build Status](https://travis-ci.org/kbsezginel/IPMOF.svg?branch=webpage)](https://travis-ci.org/kbsezginel/IPMOF)
# IPMOF
### Discovering interpenetration in MOFs
![alt text][Fig0]

## Installation
IPMOF uses Python 3.5.1 with required libraries listed in requirements.txt file.

You can install IPMOF by cloning the repository and running setup.py as follows:

`git clone https://github.com/kbsezginel/IPMOF.git`

`cd IPMOF`

`python setup.py install`

If you wish to use HPC capabilities you need to install other dependencies:

`pip install -r requirements_hpc.txt`

## Usage
IPMOF reads structure files from the mof folder in root directory.

The default file reader is ASE and IPMOF supports all the file formats supported by ASE.
A list of the formats supported by ASE is given [here][ASE-formats].

In order to test IPMOF there is one MOF file (MOF-5 / REFCODE: SAHYIK) in the 'mof' directory.
Using default simulation parameters homo-interpenetrated MOF-5 structure can be generated in 'results' directory.

## Algorithm Description
![alt text][Animation1]

IPMOF tests whether two given MOFs can interpenetrate each other by rapidly trying different relative orientations of the two frameworks and reports the plausibly energetically favorable ones. The algorithm tries many different orientations of two given MOFs by performing rotation and translation operations according to user configurable parameters. After an orientation is chosen, its energetic favorability is calculated based on the pairwise interactions between each atom on one framework with every atom on the other framework. Overall, IPMOF can rapidly detect cases where interpenetration is impossible, and suggest ones where it may be plausible. A flowchart is provided below.

## Energy Map
![alt text][Fig1]

An energy map is a regular grid of points representing the potential energy inside the crystal unit cell, from the perspective of an atom being inserted into that space. The potential energy for each point in the grid is calculated using [Lennard-Jones potential][LJ-Wikipedia] with built-in parameters provided for [Universal Force Field (UFF)][UFF-ref] and [DREIDING][DRE-ref] force fields.

Energy map is exported as a numpy array consisting of 3D coordinates and energy values for different types of atoms. There are four options that can be used to generate the atom list for a particular energy map:
- _full_: Full atom list in the force field parameters database. (103 atoms)
- _uniq_: Unique atoms for a given list of MOFs
- _dummy_: Single dummy atom with predefined force field parameters
- _qnd_: Simplified atom list consisting of 1 dummy atom and 10 most frequent atoms in CoRE MOF database

The runtime of the energy map generation and interpenetration tests depend on the selection of atom list in the following
fashion: _full_ > _uniq_ > _qnd_ > _dummy_. All list types except for _uniq_ allows checking interpenetration of the passive MOF with any other MOF. The _uniq_ energy map uses only the atoms in the given MOF set, therefore any other MOF that has other types of atoms requires new energy map. _qnd_ and _dummy_ enery map types were designed to work with large scale screening calculations to have one type of energy map for large amount of MOF combinations. _full_ energy map can also be used in the same application to increase accuracy however that would result in slower runtime and bigger file sizes for energy maps.

**To generate energy map type following in a command-line window:**

`python ipmof_energymap.py`

By default this will create energy maps for each MOF file in _~/mof_ directory. The atom list and energy map are stored as a numpy array. This can be changed to a human readable format (yaml) by changing the _energy_map_type_ simulation parameter to _yaml_.

## Interpenetration
![alt text][Fig2]
After generating the energy map for the _passive_ MOF, an _active_ MOF is selected to test interpenetration. The _active_ unit cell is first aligned with the _passive_ unit cell and then rotated around the global _x, y,_ and _z_ axes by increments defined by the user. Then the _active_ unit cell is translated in _x, y,_ and _z_ directions, also by increments defined by the user, within the _passive_ unit cell. After testing all orientations the most favorable ones (if any) are reported.

Energy map is required only for the _passive_ MOF (should be in _~/energymap_) and structure files (such as _cif_) are required for both MOFs (should be in _~/mof_).

**To test intepenetration type following in a command-line window:**

`python ipmof_interpenetration.py`

Simulation summary, simulation parameters, and information on discovered structures will be exported to _~/results/*Structure1_Structure2*/results.yaml_. Methods to analyze results are included in _~ipmof/analysis.py_ library. The structures discovered will be exported to _~/results/*Structure1_Structure2*_.

### Simulation Parameters
![alt text][Fig5]

**Grid size:** The grid size if defined as the distance between grid points of the energy map in each dimension. The energy map is generated as a 3D rectangular grid that surrounds the unit cell. For high-throughput screening we used grid size of 1 Å for the generation of energy map.

**Interpolation:** When an energy value is needed at a position that does not fall exactly on a grid point in the energy map, which is the usual case, interpolation is used. For high-throughput screening we used trilinear interpolation which approximates the value of a point in a regular grid linearly using data in the lattice points.

**Cut-off distance:** The cut-off distance used in the calculation of interatomic potentials. Any atom pair that are further away from this distance were assumed to have no interaction. For high-throughput screening cut-off distance was taken as 12 Å

**Extension cut-off distance:** The extension cut-off distance is used to determine the number of extension that will be performed on the unit cell for multiple unit cell collision test. For a given distance the unit cell is extended in each dimension to make sure that it covers an area at least the size of a sphere with radius of that distance. For high-throughput screening extension cut-off was taken as 50 Å.

**Rotational freedom:** The increments of rotation performed on the active unit cell to obtain different orientations. For high-throughput screening only 90° rotations were considered.

**Rotational limit:** The number of times a rotation is performed for same initial coordinate. For high-throughput screening all possible 90° rotations (24 total) were performed for each point.

**Energy limit:** There are three types of energy limit defined in the algorithm: atom energy limit, structure energy limit, and energy density limit.
- Atom energy limit is the maximum allowed energy for insertion of a single atom.
- Structure energy limit is the maximum allowed energy for insertion of a collection of atoms that constitute the active structure. It is calculated by summing the insertion energies for each atom in the unit cell.
- Energy density limit is similar to structure energy limit however it is divided by the unit cell volume to make it more comparable among different structures. For unit cells with different sizes constituting different number of atoms we chose to have a universal energy limit scaled by the volume of the cell. The energy density, _ρ<sub>energy</sub>_, is calculated according to equation below,

   <img src="https://github.com/kbsezginel/IPMOF/blob/master/doc/Figures/Fig6.PNG" width="300">

where _N<sub>p</sub>_ and _N<sub>a</sub>_ are number of atoms in active and passive MOFs; _i_ and _j_ are atom indices for passive and active MOFs; _V<sub>LJ</sub>_ is interatomic potential (Lennard-Jones) energy, and _V<sub>cell</sub>_ is unit cell volume in Å<sup>3</sup>.

There are additional algorithm parameters available for user to manipulate how the structures are generated, how the data is reported etc. however these parameters do not affect the interpenetration. More information on algorithm parameters can be found in algorithm documentation.

Simulation parameters can be changed by modifying _~/settings/sim_par.sample.yaml_ file.
The file name must be changed to _sim_par.yaml_ after modification in order to be read by the algorithm.
Without the _sim_par.yaml_ file default simulation parameters are read from _~/ipmof/parameters.py_.

### Simulation directories
Default directories to read input files such as energy maps, MOF files and output result files
are read from _~/ipmof/parameters.py_.

If you wish to use different directories export sim_dir.yaml file into _~/settings_ folder with
the directories you want to use. You can use methods in _~/ipmof/parameters.py_ to export this file.

### Documents
To get a better understanding of the functions used in ipmof libraries and/or analyze your
results you can go through jupyter notebook files in _~/doc/Notebooks_ directory.

Notebooks list:
- MOF-5 Example Run (Energymap generation and interpenetration test)
- Hetero-interpenetration Test (LEHXUT + XAMDUM02)
- Supercell Generation (LEHXUT + XAMDUM02)
- [CoRE Database][CORE-ref] Usage

The _~/doc_ directory also contains force field parameters (UFF and DRE) and supplementary information for
MOFs in [CoRE database][CORE-ref].

#### Candidate Structures
Candidate interpenetrated structure files can be accessed here -> [IPMOF-candidates](https://github.com/kbsezginel/IPMOF-candidates)

### Contact
For any questions, ideas or feedback please contact me!

mail: kbs37@pitt.edu

http://wilmerlab.com/

### Algorithm Flowsheet
![alt text][Fig3]

[Fig0]: https://github.com/kbsezginel/IPMOF/blob/master/doc/Figures/Fig0.PNG "CandidateStructures"
[Fig1]: https://github.com/kbsezginel/IPMOF/blob/master/doc/Figures/Fig1.PNG "Energymap"
[Fig2]: https://github.com/kbsezginel/IPMOF/blob/master/doc/Figures/Fig2.PNG "Interpenetration"
[Fig3]: https://github.com/kbsezginel/IPMOF/blob/master/doc/Figures/Fig3.PNG "AlgorithmFlowsheet"
[Fig5]: https://github.com/kbsezginel/IPMOF/blob/master/doc/Figures/Fig5.PNG "SimulationParameters"
[Animation1]:  https://github.com/kbsezginel/IPMOF/blob/master/doc/Figures/UQOFOX_VEHJUP_rotation.gif "Rotation"

[ASE-formats]: https://wiki.fysik.dtu.dk/ase/ase/io/io.html
[LJ-wikipedia]: https://en.wikipedia.org/wiki/Lennard-Jones_potential
[UFF-ref]: http://pubs.acs.org/doi/abs/10.1021/ja00051a040
[DRE-ref]: http://pubs.acs.org/doi/abs/10.1021/j100389a010
[CORE-ref]: http://pubs.acs.org/doi/abs/10.1021/cm502594j
