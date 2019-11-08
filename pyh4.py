from ctypes import Structure, POINTER, c_long, c_double, cdll

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


class H4Library(object):
    def __init__(self):
        self.library = cdll.LoadLibrary('./libh4.so')

        self.library.energy_corr_h4.argtypes = {
            c_long,
            POINTER(Atom),
            POINTER(Coord),
        }

        self.library.energy_corr_h4.restype = c_double

        self.library.energy_corr_hh_rep.argtypes = {
            c_long,
            POINTER(Atom),
            POINTER(Coord),
        }

        self.library.energy_corr_hh_rep.restype = c_double

    def H4Correction(self, natoms, positions, numbers):

        atom = np.zeros(natoms, dtype=[('x', c_double), ('y', c_double), ('z', c_double), ('e', c_long)])
        atom['x'] = positions[:, 0]
        atom['y'] = positions[:, 1]
        atom['z'] = positions[:, 2]
        atom['e'] = numbers
        gradient = np.zeros((natoms, 3), dtype=c_double)

        args = [
            c_long(natoms),
            atom.ctypes.data_as(POINTER(Atom)),
            gradient.ctypes.data_as(POINTER(Coord)),
        ]

        energy = self.library.energy_corr_h4(*args)

        return energy, gradient

    def HHRepulsion(self, natoms, positions, numbers):
        atom = np.zeros(natoms, dtype=[('x', c_double), ('y', c_double), ('z', c_double), ('e', c_long)])
        atom['x'] = positions[:, 0]
        atom['y'] = positions[:, 1]
        atom['z'] = positions[:, 2]
        atom['e'] = numbers
        gradient = np.zeros((natoms, 3), dtype=c_double)

        args = [
            c_long(natoms),
            atom.ctypes.data_as(POINTER(Atom)),
            gradient.ctypes.data_as(POINTER(Coord)),
        ]

        energy = self.library.energy_corr_hh_rep(*args)

        return energy, gradient

    def H4Calculation(self, natoms, positions, numbers):
        energy, gradient = self.H4Correction(natoms, positions, numbers)
        energy2, gradient2 = self.HHRepulsion(natoms, positions, numbers)

        return energy + energy2, gradient + gradient2
