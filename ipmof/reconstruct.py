# IPMOF Functions for reconstructing unit cell and finding common unit cell
# Date: January 2017
# Author: Kutay B. Sezginel
import math
from ipmof.geometry import car2frac, frac2car, sub3, add3, xyz_rotation


def lcm(x, tolerance=1):
    """
    Calculates least common multiplier for 1 and a given number with a tolerance given in percentage
    """
    tolerance = 100 / tolerance
    multiplier = max(int(x), 1)
    multiplier_found = False
    while(not multiplier_found):
        ratio = round(multiplier / x)
        lower = ratio - ratio / tolerance
        upper = ratio + ratio / tolerance
        if(multiplier / x <= upper and multiplier / x >= lower):
            multiplier_found = True
            break
        else:
            multiplier += 1
    return multiplier


def common_cell_parameters(mof1, mof2, tolerance=1):
    """
    Calculates common cell parameters for two given MOF objects.
    Scaling factor and distortion for the second unit cell is returned as well.
    """
    if not hasattr(mof2, 'edge_points'):
        mof2.extend_unit_cell(pack=[1, 1, 1])
    v111 = mof2.edge_points[-1]                         # 111 vector for UC2
    v111 = car2frac(v111, mof1.to_frac)                 # Convert to fractional for UC1
    lcm1 = [lcm(i, tolerance=tolerance) for i in v111]  # Calculate lcm for each dimension
    lcm2 = [lcm1[i] / v for i, v in enumerate(v111)]    # Divide lcm to vector for lcm2 -> actual lcm2
    lcm2_round = [int(round(i)) for i in lcm2]          # Calculate lcm round -> lcm2
    scaling_factor = [round(i) / i for i in lcm2]       # Scaling factor for UC2
    distortion = sum([abs((i - round(i))) / round(i) * 100 for i in lcm2])  # Cell distortion percent
    return dict(lcm1=lcm1, lcm2=lcm2_round, lcm2_actual=lcm2, sf=scaling_factor, dist=distortion)


def reshape(mof, rotation, initial_coordinate):
    """
    Apply rotation and translation operations to given MOF
        - rotation = [x_angle, y_angle, z_angle] in degrees
        - initial_coordinate = [x, y, z] in cartesian coordinates
    """
    first_point = initial_coordinate
    x_angle, y_angle, z_angle = [math.radians(a) for a in rotation]
    structure = {'atom_names': [], 'atom_coors': [], 'name': mof.name}
    rot_coor = mof.packed_coors[0][0]
    translation_vector = sub3(first_point, rot_coor)
    for unit_cell in mof.packed_coors:
        for atom_coor in unit_cell:
            rot_coor = atom_coor
            rot_coor = xyz_rotation(rot_coor, [x_angle, y_angle, z_angle])
            new_coor = add3(rot_coor, translation_vector)
            structure['atom_coors'].append(new_coor)
    structure['atom_names'] = len(mof.packed_coors) * mof.atom_names
    return structure


def rescale(mof, scaling_factor):
    """
    Rescale unit cell of a MOF object in a, b, c vectors by given amount in list format
    """
    structure = {'atom_names': [], 'atom_coors': [], 'name': mof.name}
    if not hasattr(mof, 'packing_factor'):
        mof.extend_unit_cell(pack=[1, 1, 1])
    to_frac = mof.to_frac
    to_car = mof.to_car
    for unit_cell in mof.packed_coors:
        for atom_coor in unit_cell:
            # Convert to fractional
            frac_coor = car2frac(atom_coor, to_frac)
            # Scale
            scaled_coor = [i * j for i, j in zip(scaling_factor, frac_coor)]
            # Apply pbc (using common unit cell parameters)
            pbc_coor = [i - math.floor(i) for i in scaled_coor]
            # Convert to cartesian
            new_coor = frac2car(pbc_coor, to_car)
            structure['atom_coors'].append(new_coor)
    structure['atom_names'] = len(mof.packed_coors) * mof.atom_names
    return structure


def get_common_cell(mof1, mof2, rotation, initial_coordinate, export_dir=os.getcwd(), colorify=False):
    """
    Export common cell for two given MOFs with rotation ans translation for active MOF.
        - rotation = [x_angle, y_angle, z_angle] -> Rotate atoms in degrees around each axis.
        - initial_coordinate = [x, y, z] -> Translate atoms in cartesian coordinates.
        - colorify=True/False -> Export unit cells with assigned atom name for each cell.
    """
    # 1. Get common cell parameters
    common_cell = common_cell_parameters(mof1, mof2)
    # 2. Extend both MOFs according to new cell parameters
    extended1 = mof1.extend_unit_cell(pack=common_cell['lcm1'])
    extended2 = mof2.extend_unit_cell(pack=common_cell['lcm2'])
    # 3. Apply rotation and translation to second MOF
    extended2 = reshape(mof2, rotation, initial_coordinate)
    # 4. Scale second MOF to fit rounded least common multiplier
    common_uc_size = [i * j for i, j in zip(mof1.uc_size, common_cell['lcm1'])]
    common_uc_angle = mof1.uc_angle
    new_mof2 = MOF(extended2, file_format='dict')
    new_mof2.uc_size = common_uc_size
    new_mof2.uc_angle = common_uc_angle
    new_mof2.unit_cell_volume()
    new_mof2.pbc_parameters()
    extended2 = rescale(new_mof2, common_cell['sf'])
    new_mof2 = MOF(extended2, file_format='dict')
    # 5. Join extended MOFs into single common cell and export
    new_mof1 = mof1.clone()
    if len(mof1.packed_coors) > 1:
        for cell in mof1.packed_coors[1:]:
            for coor, name in zip(cell, mof1.atom_names):
                atom = Atom(symbol=name, position=coor)
                mof1.ase_atoms.append(atom)
    a, b, c = common_uc_size
    alpha, beta, gamma = common_uc_angle
    if colorify:
        joined_mof_color = new_mof1.join(new_mof2, colorify=True)
        joined_mof_color.name = '%s_%sJC' % (mof1.name, mof2.name)
        joined_mof_color.ase_atoms.set_cell([a, b, c, alpha, beta, gamma])
        joined_mof_color.export(export_dir, file_format='cif')

    else:
        joined_mof = new_mof1.join(new_mof2, colorify=False)
        joined_mof.name = '%s_%sJ' % (mof1.name, mof2.name)
        joined_mof.ase_atoms.set_cell([a, b, c, alpha, beta, gamma])
        joined_mof.export(export_dir, file_format='cif')
