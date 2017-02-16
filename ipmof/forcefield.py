# IPMOF Force Field Functions
# Date: June 2016
# Author: Kutay B. Sezginel
import os
import math

import xlrd
import numpy as np


def read_ff_parameters(excel_file_path, ff_selection):
    """
    Read force field parameters from an excel file according to force field selection
    """
    # Read Excel File
    force_field_data = xlrd.open_workbook(excel_file_path)
    # Read columns to acquire force field parameters
    atom_names = force_field_data.sheets()[0].col_values(0)[2:]
    uff_sigma = force_field_data.sheets()[0].col_values(1)[2:]
    uff_epsilon = force_field_data.sheets()[0].col_values(2)[2:]
    dre_sigma = force_field_data.sheets()[0].col_values(3)[2:]
    dre_epsilon = force_field_data.sheets()[0].col_values(4)[2:]

    uff = {'atom': atom_names, 'sigma': uff_sigma, 'epsilon': uff_epsilon}
    dre = {'atom': atom_names, 'sigma': dre_sigma, 'epsilon': dre_epsilon}

    if ff_selection == 'uff':
        return uff
    if ff_selection == 'dre':
        return dre
    else:
        print('No such force field')


def get_ff_parameters(atom_names, ff_parameters):
    """
    Get force field parameters of the atom list and force field parameters you provide
    """
    atom_ff_parameters = []
    for atom in atom_names:
        for ff_index, ff_atom in enumerate(ff_parameters['atom']):
            if atom == ff_atom:
                atom_name = atom
                sigma = ff_parameters['sigma'][ff_index]
                epsilon = ff_parameters['epsilon'][ff_index]
                atom_ff_parameters.append([atom_name, sigma, epsilon])
    return atom_ff_parameters


def lorentz_berthelot_mix(sigmaList1, sigmaList2, epsilonList1, epsilonList2):
    """
    Lorentz-Berthelot mixing rules for given lists of sigma1, sigma2, epsilon1, and epsilon2
    """
    sig = np.zeros([len(sigmaList1), len(sigmaList2)])
    eps = np.zeros([len(epsilonList1), len(epsilonList2)])

    for index1 in range(len(sigmaList1)):
        for index2 in range(len(sigmaList2)):
            sig[index1][index2] = (sigmaList1[index1] + sigmaList2[index2]) / 2
            eps[index1][index2] = math.sqrt(epsilonList1[index1] * epsilonList2[index2])

    return sig, eps


def lennard_jones(r, sig, eps):
    """
    Calculate Lennard Jones potential for given distance, sigma, and epsilon values.
    Energy unit: (kB)
    """
    return 4 * eps * ((sig / r)**12 - (sig / r)**6)
