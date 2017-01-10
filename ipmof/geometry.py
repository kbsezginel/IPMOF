# Python library for 3 dimensional geometric operations
# Date: October 2016
# Author: Kutay B. Sezginel
import math
from random import random

import mathutils


def rotation(p, a1, a2, angle):
    """
    Quaternion rotation using mathutils library.
    Rotation is performed for a given point and angle around an axis defined by two given points.
     >>> new_point = rotation(rotation_point, axis_point1, axis_point2, rotation_angle)
            p: point of rotation
            a1: first point for rotation axis
            a2: second point for rotation axis
            angle: angle of rotaion in radians
    """
    Q_point = mathutils.Quaternion((0, p[0] - a2[0], p[1] - a2[1], p[2] - a2[2]))
    Q_rot = mathutils.Quaternion((a2[0] - a1[0], a2[1] - a1[1], a2[2] - a1[2]), angle).normalized()
    Quat = (Q_rot * Q_point) * Q_rot.inverted()

    return [Quat.x + a2[0], Quat.y + a2[1], Quat.z + a2[2]]


def xyz_rotation(p, angle):
    """
    Quaternion rotation using mathutils library.
    Rotation is performed for a given point and angle list around x, y, and z axes respectively.
     >>> new_point = rotation(rotation_point, rotation_angle_list)
            p: point of rotation
            rotation_angle_list: angle of rotations for x, y, z axes in radians.
    """

    Q_point = mathutils.Quaternion((0, p[0], p[1], p[2]))
    Qrx = mathutils.Quaternion((1, 0, 0), angle[0])
    Qry = mathutils.Quaternion((0, 1, 0), angle[1])
    Qrz = mathutils.Quaternion((0, 0, 1), angle[2])
    Quat = (Qrz * Qry * Qrx * Q_point) * Qrx.inverted() * Qry.inverted() * Qrz.inverted()

    return [Quat.x, Quat.y, Quat.z]


def add3(p1, p2):
    """ 3D vector addition. """
    return [p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]]


def sub3(p1, p2):
    """ 3D vector subtraction. """
    return [p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]]


def pbc3(pbc_coor, to_frac, to_car):
    """ 3D vector periodic boundary conditions. """
    x, y, z = pbc_coor
    x_frac = to_frac[0] * x + to_frac[1] * y + to_frac[2] * z
    x_frac -= math.floor(x_frac)
    y_frac = to_frac[3] * y + to_frac[4] * z
    y_frac -= math.floor(y_frac)
    z_frac = to_frac[5] * z
    z_frac -= math.floor(z_frac)
    x_frac = to_car[0] * x_frac + to_car[1] * y_frac + to_car[2] * z_frac
    y_frac = to_car[3] * y_frac + to_car[4] * z_frac
    z_frac = to_car[5] * z_frac

    return [x_frac, y_frac, z_frac]


def car2frac(car_coor, to_frac):
    """ Convert cartesian coordinates to fractional coordinates.
        Requires 'to_frac' constants which is calculated for MOF objects. """
    x, y, z = car_coor
    x_frac = to_frac[0] * x + to_frac[1] * y + to_frac[2] * z
    y_frac = to_frac[3] * y + to_frac[4] * z
    z_frac = to_frac[5] * z
    return [x_frac, y_frac, z_frac]


def frac2car(frac_coor, to_car):
    """ Convert fractional coordinates to cartesian coordinates.
        Requires 'to_car' constants which is calculated for MOF objects. """
    x_frac, y_frac, z_frac = frac_coor
    x = to_car[0] * x_frac + to_car[1] * y_frac + to_car[2] * z_frac
    y = to_car[3] * y_frac + to_car[4] * z_frac
    z = to_car[5] * z_frac
    return [x, y, z]


def possible_rotations(rot_degree, num_of_points=10):
    """ Calculate possible degrees of rotation about x, y, z axes for given rotation freedom. """
    # Calculate possible rotation angles for the given rotation freedom
    degrees = []
    increments = round(360 / rot_degree)
    for i in range(increments):
        degrees.append(2 * math.pi * i / increments)

    # Rotate an arbitrary point around x, y, and z axes with combinations of possible rot angles
    p = [-34.34, -3.1231, 31.123]
    points = []
    degree_combinations = []
    for dx in degrees:
        for dy in degrees:
            for dz in degrees:
                degree_combinations.append([dx, dy, dz])
                points.append(xyz_rotation(p, [dx, dy, dz]))

    # Look for points that are equal to each other after rotations were performed
    eq_list = []
    for p1_index, p1 in enumerate(points):
        for p2_index, p2 in enumerate(points):
            if p2_index > p1_index:
                p1r = [round(i, 3) for i in p1]
                p2r = [round(i, 3) for i in p2]
                if p1r == p2r and [p2_index, p1_index] not in eq_list:
                    eq_list.append([p1_index, p2_index])

    # Identify all identical points (in the above step more than one equalities can be returned)
    eq1 = [i[0] for i in eq_list]
    eq2 = [i[1] for i in eq_list]
    new_eq = []
    for eq in eq1:
        if eq not in eq2 and eq not in new_eq:
            new_eq.append(eq)

    # Gather rot angle combinations that give unique results
    possible_rot_degrees = []
    for d_index, d in enumerate(degree_combinations):
        if d_index in new_eq:
            possible_rot_degrees.append(d)

    return possible_rot_degrees
