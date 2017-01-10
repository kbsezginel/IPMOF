# IPMOF Functions for reconstructing unit cell and finding common unit cell
# Date: January 2017
# Author: Kutay B. Sezginel


def lcm(x, tolerance=1):
    """ Calculates least common multiplier for 1 and a given number
        with a tolerance given in percentage """
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
    v111 = frac3(v111, mof1.to_frac)                    # Convert to fractional for UC1
    lcm1 = [lcm(i, tolerance=tolerance) for i in v111]  # Calculate lcm for each dimension
    lcm2 = [lcm1[i] / v for i, v in enumerate(v111)]    # Divide lcm to vector for lcm2 -> actual lcm2
    lcm2_round = [int(round(i)) for i in lcm2]          # Calculate lcm round -> lcm2
    scaling_factor = [round(i) / i for i in lcm2]       # Scaling factor for UC2
    distortion = sum([abs((i - round(i))) / round(i) * 100 for i in lcm2])  # Cell distortion percent
    return dict(lcm1=lcm1, lcm2=lcm2_round, lcm2_actual=lcm2, sf=scaling_factor, dist=distortion)


def reshape(mof, rotation, initial_coordinate):
    """ Apply rotation and translation operations to given MOF """
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
