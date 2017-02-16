# IPMOF CoRE Database Functions
# Date: June 2016
# Author: Kutay B. Sezginel
import os

import xlrd
import yaml
import math


def core_mof_properties(core_file_path):
    """
    Reads properties of MOFs in the CoRE database from 'CoRE.xlsx' file in 'doc' folder.
    Properties are exported in dictonary format with following inputs:
    - mof_names - metal_type - unmofied_mofs
    - density (g/cm3) - pore_volume (cm3/g) - void_fraction
    - pld (pore limiting diamater)
    - lcd (largest cavity diamater)
    - vsa (volumetric surface area m2/cm3)
    - gsa (gravimetric surface area m2/g)
    CoRE Database: 10.1021/cm502594j
    """
    # Read Excel File
    core_data = xlrd.open_workbook(core_file_path)
    # Read columns to acquire force field parameters
    mof_names = core_data.sheets()[0].col_values(1)[9:]
    metal_type = core_data.sheets()[0].col_values(2)[9:]
    density = core_data.sheets()[0].col_values(3)[9:]
    pld = core_data.sheets()[0].col_values(4)[9:]
    lcd = core_data.sheets()[0].col_values(5)[9:]
    vsa = core_data.sheets()[0].col_values(6)[9:]
    gsa = core_data.sheets()[0].col_values(7)[9:]
    void_fraction = core_data.sheets()[0].col_values(8)[9:]
    pore_volume = core_data.sheets()[0].col_values(9)[9:]
    unmofied_mofs = core_data.sheets()[9].col_values(1)[3:]

    mof_properties = {'mof_name': mof_names,
                      'metal_type': metal_type,
                      'density': density,
                      'pld': pld,
                      'lcd': lcd,
                      'vsa': vsa,
                      'gsa': gsa,
                      'void_fraction': void_fraction,
                      'pore_volume': pore_volume,
                      'unmofied_mofs': unmofied_mofs
                      }

    return mof_properties


def core_mof_sort(mof_properties, sort='void_fraction', limit='0.9'):
    """
    Selects MOFs with given property higher than given limit.
    Default property is void fraction with a limit of 0.9.
    Available properties are:
        - density (g/cm3) - pore_volume (cm3/g) - void_fraction
        - pld (pore limiting diamater) - lcd (largest cavity diamater)
        - vsa (volumetric surface area m2/cm3)
        - gsa (gravimetric surface area m2/g)
    """
    sorted_mofs = {sort: [], 'name': []}
    sort_prop = mof_properties[sort]

    for prop_index, prop in enumerate(sort_prop):
        if prop > limit:
            sorted_mofs[sort].append(prop)
            sorted_mofs['name'].append(mof_properties['mof_name'][prop_index])

    print('Gathered a total of', len(sorted_mofs['name']), 'MOFs')
    print('With', sort, '>', limit)

    return sorted_mofs


def core_mof_dir(sorted_mofs, core_mof_dir):
    """
    Collects directories for the given MOF list (sorted_mofs).
    The directories are taken from the database directory (core_mof_dir) provided.
    Collected directories are returned in list format.
    """
    all_mofs = os.listdir(core_mof_dir)
    mof_dirs = []
    missing_mofs = []

    for mof in sorted_mofs['name']:
        mof_found = False
        for mof_file_name in all_mofs:
            mof_name = mof_file_name.split('.')[0].split('_')[0]
            if mof == mof_name:
                mof_dir = os.path.join(core_mof_dir, mof_file_name)
                mof_dirs.append(mof_dir)
                mof_found = True
        if not mof_found:
            missing_mofs.append(mof)
            print(mof)

    print(len(missing_mofs), 'mofs are missing |^|')
    print(len(mof_dirs), 'total mofs found')

    return mof_dirs


def core_mof_vf_list(target_vf, vf_list_path, limit=math.inf):
    """ Generate MOF combination list from void fraction values """
    core_vf = yaml.load(open(vf_list_path, 'r'))
    vf_mofs = []
    count = 0
    vf_count = 0
    list_complete = False
    print('Reading CoRE void fraction list...')
    for mof1, vf1 in enumerate(core_vf):
        if not list_complete:
            for mof2, vf2 in enumerate(core_vf):
                if mof2 >= mof1:
                    count += 1
                    if vf1[1] + vf2[1] >= target_vf:
                        if vf_count < limit:
                            vf_mofs.append([vf1[0], vf2[0], vf1[1] + vf2[1]])
                            vf_count += 1
                        else:
                            break
                            list_complete = True
                        # print("\rCounting MOF combinations with total Vf > %d | %d / %d selected" % (target_vf, vf_count, count), end="")
    print("MOF combinations with Vf > %d : %d / %d selected. Limit: %.0f" % (target_vf, vf_count, count, limit), end="")

    return vf_mofs


def core_interpenetration_list(sim_dir, limit=math.inf, target_vf=1.0):
    """ Generate interpenetration list from a given MOF combination list """
    vf_list_path = os.path.join(sim_dir['main_dir'], 'doc', 'core_mof_vf_list.yaml')
    mof_list = core_mof_vf_list(target_vf, vf_list_path, limit=limit)
    emap_dir = sim_dir['energy_map_dir']
    mof_dir = sim_dir['mof_dir']
    interpenetration_list = []
    for mof in mof_list:
        emap_path = os.path.join(emap_dir, "%s_emap.npy" % mof[0])
        emap_mof_path = os.path.join(mof_dir, "%s.cif" % mof[0])
        ip_mof_path = os.path.join(mof_dir, "%s.cif" % mof[1])
        interpenetration_list.append((emap_path, emap_mof_path, ip_mof_path))
    return interpenetration_list
