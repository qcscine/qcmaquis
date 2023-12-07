

import pytest
from pyscf import gto, mcscf, scf

from scine_qcmaquis.pyscf_interface.pyscf_interface import QcMaquis


def test_dmrgci():
    mol = gto.Mole()
    mol.build(
        atom=[['H', (0., 0., i)] for i in range(8)],
        basis={'H': '6-31g'},
        symmetry=False, spin=2, verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASCI(m, 4, 4)
    mc.fcisolver = QcMaquis(mol)
    mc.fcisolver.file_path = None
    emc_0 = mc.casci()[0]

    mol = gto.Mole()
    mol.build(
        atom=[['H', (0., 0., i)] for i in range(8)],
        basis={'H': '6-31g'},
        symmetry=False, spin=2, verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()
    mc = mcscf.CASCI(m, 4, 4)
    emc_0ref = mc.casci()[0]
    # print('Maquis-CI  = %.15g CASCI  = %.15g' % (emc_0, emc_0ref))
    # print('Diff = %.15g' % (emc_0 - emc_0ref))
    assert abs(emc_0ref - emc_0) < 1e-9


def test_dmrgscf():
    mol = gto.Mole()
    mol.build(
        atom=[['H', (0., 0., i)] for i in range(8)],
        basis={'H': '6-31g'},
        symmetry=False, spin=2, verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASSCF(m, 4, 4)
    mc.max_cycle_macro = 20
    mc.fcisolver = QcMaquis(mol)
    mc.fcisolver.file_path = None
    emc_1 = mc.mc2step()[0]

    mol = gto.Mole()
    mol.build(
        atom=[['H', (0., 0., i)] for i in range(8)],
        basis={'H': '6-31g'},
        symmetry=False, spin=2, verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()
    mc = mcscf.CASSCF(m, 4, 4)
    mc.max_cycle_macro = 20
    emc_1ref = mc.mc2step()[0]
    # print('Maquis-SCF = %.15g CASSCF = %.15g' % (emc_1, emc_1ref))
    # print('Diff = %.15g' % (emc_1 - emc_1ref))
    assert abs(emc_1ref - emc_1) < 1e-9
