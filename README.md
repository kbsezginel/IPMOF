# IPMOF
### Discovering interpenetration in MOFs
![alt text][Fig0]

## Installation
IPMOF uses Python 3.5.1 with required libraries listed in requirements.txt file.

To install dependencies just type:

`pip install -r requirements.txt`

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
- 'full': Full atom list in the force field parameters database. (103 atoms)
- 'uniq': Unique atoms for a given list of MOFs
- 'dummy': Single dummy atom with predefined force field parameters
- 'qnd': Simplified atom list consisting of 1 dummy atom and 10 most frequent atoms in CoRE MOF database

The runtime of the energy map generation and interpenetration tests depend on the selection of atom list in the following
fashion: 'full' > 'uniq' > 'qnd' > 'dummy'. All list types except for 'uniq' allows checking interpenetration of the passive MOF with any other MOF. The 'uniq' energy map uses only the atoms in the given MOF set, therefore any other MOF that has other types of atoms requires new energy map. 'qnd' and 'dummy' enery map types were designed to work with large scale screening calculations to have one type of energy map for large amount of MOF combinations. 'full' energy map can also be used in the same application to increase accuracy however that would result in slower runtime and bigger file sizes for energy maps.

**To generate energy map type following in a command-line window:**

`python ipmof_energymap.py`

By default this will create energy maps for each MOF file in ~/mof directory. The atom list and energy map are stored as a numpy array. This can be changed to a human readable format (yaml) by changing the 'energy_map_type' simulation parameter to 'yaml'.

## Interpenetration
![alt text][Fig2]
After generating the energy map for the _passive_ MOF, an _active_ MOF is selected to test interpenetration. The _active_ unit cell is first aligned with the _passive_ unit cell and then rotated around the global x, y, and z axes by increments defined by the user. Then the _active_ unit cell is translated in x, y, and z directions, also by increments defined by the user, within the _passive_ unit cell. After testing all orientations the most favorable ones (if any) are reported.

Energy map is required only for the _passive_ MOF (should be in ~/energymap) and structure files (such as 'cif') are required for both MOFs (should be in ~/mof).

**To test intepenetration type following in a command-line window:**

`python ipmof_interpenetration.py`

Simulation summary, simulation parameters, and information on discovered structures will be exported to ~/results/*Structure1_Structure2*/results.yaml. Methods to analyze results are included in ~ipmof/analysis.py library. The structures discovered will be exported to ~/results/*Structure1_Structure2*.

### Simulation Parameters
Default simulation parameters are read from ~/ipmof/parameters.py.

To change the parameters either modify the parameters in that file or create a 'yaml'
file using the functions in ipmof.parameters library. Necessary functions for creating
simulation parameters input files and reading them are given in this library with instructions.
Also, a sample simulation parameters file (sim_par_sample.yaml) is given in ~/settings.
This file can also be used by modifying its name to sim_par.yaml in the same directory.

### Simulation directories
Default directories to read input files such as energy maps, MOF files and output result files
are read from ~/ipmof/parameters.py.

If you wish to use different directories export sim_dir.yaml file into ~/setting folder with
the directories you want to use. You can use methods in ~/ipmof/parameters.py to export this file.

### Documents
To get a better understanding of the functions used in ipmof libraries and/or analyze your
results you can go through jupyter notebook files in '/doc' directory.

This directory also contains force field parameters (UFF and DRE) and supplementary information for
MOFs in [CoRE database][CORE-ref].

### Contact
For any questions, ideas or feedback please contact me!

mail: kbs37@pitt.edu

http://wilmerlab.com/

### Algorithm Flowsheet
![alt text][Fig3]

[Fig0]: https://github.com/kbsezginel/IPMOF/blob/documentation/doc/Figures/Fig0.PNG "CandidateStructures"
[Fig1]: https://github.com/kbsezginel/IPMOF/blob/documentation/doc/Figures/Fig1.PNG "Energymap"
[Fig2]: https://github.com/kbsezginel/IPMOF/blob/documentation/doc/Figures/Fig2.PNG "Interpenetration"
[Fig3]: https://github.com/kbsezginel/IPMOF/blob/documentation/doc/Figures/Fig3.PNG "AlgorithmFlowsheet"
[Animation1]:  https://github.com/kbsezginel/IPMOF/blob/documentation/doc/Figures/UQOFOX_VEHJUP_rotation.gif "Rotation"

[ASE-formats]: https://wiki.fysik.dtu.dk/ase/ase/io/io.html
[LJ-wikipedia]: https://en.wikipedia.org/wiki/Lennard-Jones_potential
[UFF-ref]: http://pubs.acs.org/doi/abs/10.1021/ja00051a040
[DRE-ref]: http://pubs.acs.org/doi/abs/10.1021/j100389a010
[CORE-ref]: http://pubs.acs.org/doi/abs/10.1021/cm502594j
