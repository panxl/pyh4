import numpy as np

from pyh4 import H4Library


Elements = ["None", 'H', 'He',
            'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
            'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
            'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
            'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
            'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb',
            'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi']


# PM6-D3H4 parameters
parameters = {
    'para_oh_o': 2.32,
    'para_oh_n': 3.10,
    'para_nh_o': 1.07,
    'para_nh_n': 2.01,
    'multiplier_wh_o': 0.42,
    'multiplier_nh4': 3.61,
    'multiplier_coo': 1.41,
    'hh_rep_k': 0.4,
    'hh_rep_e': 12.7,
    'hh_rep_r0': 2.3,
}


def test_h4_correction(xyz, out):
    natoms = np.loadtxt(xyz, dtype=int, max_rows=1)
    positions = np.loadtxt(xyz, usecols=(1, 2, 3), skiprows=2, max_rows=natoms)

    elements = np.loadtxt(xyz, dtype="U2", usecols=0, skiprows=2, max_rows=natoms)
    numbers = np.zeros(natoms, dtype=int)
    for i, e in enumerate(elements):
        numbers[i] = Elements.index(e)

    lib = H4Library()

    ene, grad = lib.H4Correction(natoms, positions, numbers, parameters)
    ene_ref = np.loadtxt("TIP4P-2.out", usecols=2, max_rows=1)
    grad_ref = np.loadtxt("TIP4P-2.out", skiprows=5, max_rows=natoms)

    np.testing.assert_allclose(ene, ene_ref, atol=1e-6)
    np.testing.assert_allclose(grad, grad_ref, atol=1e-6)


def test_hh_repulsion(xyz, out):
    natoms = np.loadtxt(xyz, dtype=int, max_rows=1)
    positions = np.loadtxt(xyz, usecols=(1, 2, 3), skiprows=2, max_rows=natoms)

    elements = np.loadtxt(xyz, dtype="U2", usecols=0, skiprows=2, max_rows=natoms)
    numbers = np.zeros(natoms, dtype=int)
    for i, e in enumerate(elements):
        numbers[i] = Elements.index(e)

    lib = H4Library()

    ene, grad = lib.HHRepulsion(natoms, positions, numbers, parameters)

    ene_ref = np.loadtxt("TIP4P-2.out", usecols=2, skiprows=1, max_rows=1)
    grad_ref = np.loadtxt("TIP4P-2.out", skiprows=(natoms + 7), max_rows=natoms)

    np.testing.assert_allclose(ene, ene_ref, atol=1e-6)
    np.testing.assert_allclose(grad, grad_ref, atol=1e-6)


if __name__ == "__main__":
    test_h4_correction("TIP4P-2.xyz", "TIP4P-2.out")
    test_hh_repulsion("TIP4P-2.xyz", "TIP4P-2.out")
