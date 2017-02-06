# Create simulation input files for IPMOF simulations
import os

import yaml

# Simulation Parameters Data:
sim_par_data = {'structure_energy_limit': 1E8,   # Maximum allowed potential energy for structure
                'atom_energy_limit': 1E8,        # Maximum allowed potential energy for atom
                'energy_density_limit': 1.0,     # Maximum allowed potential energy for atom
                'rotation_limit': 20,            # Total number of rotations for each point
                'rotation_freedom': 30,          # Increments of rotation (degrees)
                'try_all_rotations': False,      # Try all possible rotations for given angle
                'summary_percent': 5,            # Percentage increment to acquire summary data
                'cut_off': 12,                   # Cut-off radius for interpenetration (Angstrom)
                'ext_cut_off': 50,               # Cut-off radius for checking extension (Angstrom)
                'check_extension': True,         # Check extended unit cells for collisions
                'grid_size': 1,                  # Grid size for potential energy map (Angstrom)
                'force_field': 'uff',            # Force field selection for LJ ('uff' or 'dre')
                'core_database': False,          # Use CoRE database information or not
                'core_limit': 1000,              # Number of MOF combinations to select from CoRE database
                'energy_map_atom_list': 'uniq',  # Atom list for energy map ('full', 'uniq', 'dummy', 'qnd')
                'energy_map_type': 'numpy',      # Energy map file format ('numpy' or 'yaml')
                'self_interpenetration': True,   # Test for homo-interpenetration or not
                'report_structures': 10,         # Number of min. energy structures to report in results
                'export_structures': 5,          # Number of min. energy structures to export
                'export_format': 'cif',          # Export structure file format
                'export_pbc': True,              # Export coordinates after applying PBC
                'export_single': True,           # Export structures with original atom names
                'export_single_color': True,     # Export structures with each I.P. layer colored
                'export_packed': True,           # Export structures with packed unit cell
                'export_packed_color': True,     # Export structures with packed unit cell and color
                'export_summary': True,          # Export summary information for interpenetration
                'directory_separation': False    # Separate results directories alphabetically
                }

# Working Directories:
main_dir = os.getcwd()
python_lib_dir = os.path.join(main_dir, 'ipmof')
force_field_path = os.path.join(main_dir, 'doc', 'FF_Parameters.xlsx')
core_path = os.path.join(main_dir, 'doc', 'CoRE.xlsx')
vf_list_path = os.path.join(main_dir, 'doc', 'core_mof_vf_list.yaml')
core_mof_dir = r'CoRE MOFs database directory'
mof_dir = os.path.join(main_dir, 'mof')
energy_map_dir = os.path.join(main_dir, 'energymap')
export_dir = os.path.join(main_dir, 'results')
settings_dir = os.path.join(main_dir, 'settings')
if not os.path.isdir(export_dir):
    os.mkdir(export_dir)
sim_par_path = os.path.join(settings_dir, 'sim_par.yaml')
sim_dir_path = os.path.join(settings_dir, 'sim_dir.yaml')

# Simulation Directories Data:
sim_dir_data = {'main_dir': main_dir,
                'python_lib_dir': python_lib_dir,
                'force_field_path': force_field_path,
                'core_path': core_path,
                'core_mof_dir': core_mof_dir,
                'mof_dir': mof_dir,
                'energy_map_dir': energy_map_dir,
                'export_dir': export_dir,
                'settings_dir': settings_dir,
                'vf_list_path': vf_list_path
                }


def read_parameters(sim_par_path=sim_par_path, sim_dir_path=sim_dir_path):
    """
    Read simulation parameters and directories from yaml files if they exist.
    Otherwise they are read from ~/ipmof/parameters.py.
    """
    if os.path.exists(sim_par_path):
        sim_par = yaml.load(open(sim_par_path, 'r'))
    else:
        sim_par = sim_par_data
    if os.path.exists(sim_dir_path):
        sim_dir = yaml.load(open(sim_dir_path, 'r'))
    else:
        sim_dir = sim_dir_data

    return sim_par, sim_dir


def export_sim_par(exp_dir=settings_dir):
    """
    Export simulation parameters input file to given directory.
    If no directory is given the file will be exported to default location defined in this library.
    """
    sim_par_path = os.path.join(exp_dir, 'sim_par.yaml')
    with open(sim_par_path, 'w') as sim_par_file:
        yaml.dump(sim_par_data, sim_par_file, default_flow_style=False)


def export_sim_dir(exp_dir=settings_dir):
    """
    Export simulation directories to given directory.
    If no directory is given the file will be exported to default location defined in this library.
    """
    sim_dir_path = os.path.join(exp_dir, 'sim_dir.yaml')
    with open(sim_dir_path, 'w') as sim_dir_file:
        yaml.dump(sim_dir_data, sim_dir_file, default_flow_style=False)


def export_interpenetration_results(sim_par, structure_info, summary, export_dir):
    """
    Export interpenetration results in yaml format.
    """
    sim_par_dict = {'simulation_parameters': sim_par}
    structure_dict = {'structure_info': structure_info}
    summary_dict = {'summary': summary}

    results_path = os.path.join(export_dir, 'results.yaml')
    with open(results_path, 'w') as yaml_file:
        yaml.dump(sim_par_dict, yaml_file, default_flow_style=False, indent=4)
        yaml.dump(structure_dict, yaml_file, explicit_start=True, indent=4)
        yaml.dump(summary_dict, yaml_file, explicit_start=True, indent=4)
