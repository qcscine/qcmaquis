import os
import shutil

import pytest
from pyscf import gto, mcscf, scf

from scine_qcmaquis.pyscf_interface.pyscf_interface import QcMaquis


def remove_checkpoint_and_results_file():
    path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(path)
    shutil.rmtree("checkpoint_DMRGSCF.h5", ignore_errors=True)
    shutil.rmtree("checkpoint_gs.h5", ignore_errors=True)
    os.remove("results_DMRGSCF.h5")
    os.remove("results_file.h5")


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
    mc.fcisolver.parameters.set_orbital_optimization()
    # mc.fcisolver.verbose = 1
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
    print('Maquis-CI  = %.15g CASCI  = %.15g' % (emc_0, emc_0ref))
    print('Diff = %.15g' % (emc_0 - emc_0ref))
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
    mc.fcisolver.parameters.set_orbital_optimization()
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
    print('Maquis-SCF = %.15g CASSCF = %.15g' % (emc_1, emc_1ref))
    print('Diff = %.15g' % (emc_1 - emc_1ref))
    assert abs(emc_1ref - emc_1) < 1e-9


def test_dmrgscf_with_checkpoint():
    mol = gto.Mole()
    mol.build(
        atom=[['N', (0., 0., i)] for i in range(8)],
        basis={'N': '6-31g'},
        symmetry=False, spin=2, verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASSCF(m, 6, 6)
    mc.max_cycle_macro = 20
    mc.fcisolver = QcMaquis(mol)
    mc.fcisolver.verbose = 1

    mc.fcisolver.file_path = os.path.dirname(os.path.abspath(__file__))
    mc.fcisolver.parameters.set_orbital_optimization()
    mc.fcisolver.parameters.set_entropies()
    emc_1 = mc.mc2step()[0]

    mol = gto.Mole()
    mol.build(
        atom=[['N', (0., 0., i)] for i in range(8)],
        basis={'N': '6-31g'},
        symmetry=False, spin=2, verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()
    mc = mcscf.CASSCF(m, 6, 6)
    mc.max_cycle_macro = 20
    emc_1ref = mc.mc2step()[0]
    print('Maquis-SCF = %.15g CASSCF = %.15g' % (emc_1, emc_1ref))
    print('Diff = %.15g' % (emc_1 - emc_1ref))
    assert abs(emc_1ref - emc_1) < 5e-9


def test_dmrgscf_with_existing_checkpoint():
    mol = gto.Mole()
    mol.build(
        atom=[['N', (0., 0., i)] for i in range(8)],
        basis={'N': '6-31g'},
        symmetry=False, spin=2, verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASSCF(m, 6, 6)
    mc.max_cycle_macro = 20
    mc.fcisolver = QcMaquis(mol)
    mc.fcisolver.verbose = 1

    mc.fcisolver.file_path = os.path.dirname(os.path.abspath(__file__))
    mc.fcisolver.parameters.set_orbital_optimization()
    mc.fcisolver.parameters.set_entropies()
    emc_1 = mc.mc2step()[0]

    mol = gto.Mole()
    mol.build(
        atom=[['N', (0., 0., i)] for i in range(8)],
        basis={'N': '6-31g'},
        symmetry=False, spin=2, verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()
    mc = mcscf.CASSCF(m, 6, 6)
    mc.max_cycle_macro = 20
    emc_1ref = mc.mc2step()[0]
    print('Maquis-SCF = %.15g CASSCF = %.15g' % (emc_1, emc_1ref))
    print('Diff = %.15g' % (emc_1 - emc_1ref))
    assert abs(emc_1ref - emc_1) < 5e-9


if __name__ == "__main__":
    test_dmrgci()
    test_dmrgscf()
    test_dmrgscf_with_checkpoint()
    test_dmrgscf_with_existing_checkpoint()
    remove_checkpoint_and_results_file()
