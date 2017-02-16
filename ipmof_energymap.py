import os
import sys

# Load IPMOF python libraries
from ipmof.forcefield import read_ff_parameters
from ipmof.energymap import energy_map, get_mof_list, energy_map_atom_list
from ipmof.parameters import read_parameters

# Read simulation parameters and directories
sim_par, sim_dir = read_parameters()

# Read excel file containing force field information
force_field = read_ff_parameters(sim_dir['force_field_path'], sim_par['force_field'])

# Read MOF list from CoRE or from a given directory
mof_path_list = get_mof_list(sim_par, sim_dir)

# Calculate atom list according to 'energy_map_atom_list' simulation parameter
atom_list = energy_map_atom_list(sim_par, force_field, mof_path_list)

# Export initialization file containing MOF names and simulation parameters
print('Starting energy map calculation for', len(mof_path_list), 'MOF(s) (grid size:',
      sim_par['grid_size'], '| cut-off radius:', str(sim_par['cut_off']) + ')')
print('Atom list ->', atom_list['atom'])
print('Energy map(s) will be exported in', sim_par['energy_map_type'], 'format')

# Main Loop (Energy Map)
for mof_index, mof_path in enumerate(mof_path_list):

    print('-' * 80)
    print(mof_index + 1, 'Calculating energy map for ->', os.path.basename(mof_path))
    # Submit jobs here
    if sys.argv[-1] == 'q':
        # Load job server libraries
        from rq import Queue
        from redis import Redis
        import sjs

        # Load job queue
        sjs.load(os.path.join("settings", "sjs.yaml"))
        job_queue = sjs.get_job_queue()

        # Calculate energy map
        job_queue.enqueue(energy_map, sim_par, mof_path, atom_list, force_field)
    else:
        emap = energy_map(sim_par, mof_path, atom_list, force_field)
