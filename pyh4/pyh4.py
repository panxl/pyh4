from pathlib import Path
from ctypes import Structure, POINTER, CDLL, c_long, c_double

import numpy as np


class Coord(Structure):
    """Cartesian coordinates"""
    _fields_ = [
        ('x', c_double),
        ('y', c_double),
        ('z', c_double),
    ]


class Atom(Structure):
    """Atom: coordinates + proton number"""
    _fields_ = [
        ('x', c_double),
        ('y', c_double),
        ('z', c_double),
        ('e', c_long),
    ]


class H4Parameters(Structure):
    """H4 Parameters"""
    _fields_ = [
        ('para_oh_o', c_double),
        ('para_oh_n', c_double),
        ('para_nh_o', c_double),
        ('para_nh_n', c_double),
        ('multiplier_wh_o', c_double),
        ('multiplier_nh4', c_double),
        ('multiplier_coo', c_double),
        ('hh_rep_k', c_double),
        ('hh_rep_e', c_double),
        ('hh_rep_r0', c_double),
    ]


class H4Library(object):
    def __init__(self):
        self.library = CDLL(Path(__file__).parent.joinpath('libh4.so'))

        self.library.energy_corr_h4.argtypes = [
            c_long,
            POINTER(Atom),
            POINTER(Coord),
            H4Parameters,
        ]

        self.library.energy_corr_h4.restype = c_double

        self.library.energy_corr_hh_rep.argtypes = [
            c_long,
            POINTER(Atom),
            POINTER(Coord),
            H4Parameters,
        ]

        self.library.energy_corr_hh_rep.restype = c_double

    def H4Correction(self, natoms, positions, numbers, parameters):

        atom = np.zeros(natoms, dtype=[('x', c_double), ('y', c_double), ('z', c_double), ('e', c_long)])
        atom['x'] = np.asarray(positions)[:, 0]
        atom['y'] = np.asarray(positions)[:, 1]
        atom['z'] = np.asarray(positions)[:, 2]
        atom['e'] = np.asarray(numbers)
        gradient = np.zeros((natoms, 3), dtype=c_double)

        args = [
            c_long(natoms),
            atom.ctypes.data_as(POINTER(Atom)),
            gradient.ctypes.data_as(POINTER(Coord)),
            H4Parameters(**parameters),
        ]

        energy = self.library.energy_corr_h4(*args)

        return energy, gradient

    def HHRepulsion(self, natoms, positions, numbers, parameters):
        atom = np.zeros(natoms, dtype=[('x', c_double), ('y', c_double), ('z', c_double), ('e', c_long)])
        atom['x'] = np.asarray(positions)[:, 0]
        atom['y'] = np.asarray(positions)[:, 1]
        atom['z'] = np.asarray(positions)[:, 2]
        atom['e'] = np.asarray(numbers)
        gradient = np.zeros((natoms, 3), dtype=c_double)

        args = [
            c_long(natoms),
            atom.ctypes.data_as(POINTER(Atom)),
            gradient.ctypes.data_as(POINTER(Coord)),
            H4Parameters(**parameters),
        ]

        energy = self.library.energy_corr_hh_rep(*args)

        return energy, gradient

    def H4Calculation(self, natoms, positions, numbers, parameters):
        energy, gradient = self.H4Correction(natoms, positions, numbers, parameters)
        energy2, gradient2 = self.HHRepulsion(natoms, positions, numbers, parameters)

        return energy + energy2, gradient + gradient2
