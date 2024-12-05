import os
import shutil

import pytest
from pyscf import gto, mcscf, scf, fci
import numpy as np

from scine_qcmaquis.pyscf_interface.pyscf_interface import DMRGSolver


def test_dmrgci():
    mol = gto.Mole()
    mol.build(
        atom=[["H", (0.0, 0.0, i)] for i in range(8)],
        basis={"H": "6-31g"},
        symmetry=False,
        spin=2,
        verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASCI(m, 4, 4)
    mc.fcisolver = DMRGSolver(mol)
    mc.fcisolver.file_path = None
    mc.fcisolver.parameters.set_orbital_optimization()
    # mc.fcisolver.verbose = 1
    emc_0 = mc.casci()[0]

    mol = gto.Mole()
    mol.build(
        atom=[["H", (0.0, 0.0, i)] for i in range(8)],
        basis={"H": "6-31g"},
        symmetry=False,
        spin=2,
        verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()
    mc = mcscf.CASCI(m, 4, 4)
    emc_0ref = mc.casci()[0]
    print("Maquis-CI  = %.15g CASCI  = %.15g" % (emc_0, emc_0ref))
    print("Diff = %.15g" % (emc_0 - emc_0ref))
    assert abs(emc_0ref - emc_0) < 1e-9
    remove_checkpoint_and_results_file()


def test_dmrgscf():
    mol = gto.Mole()
    mol.build(
        atom=[["H", (0.0, 0.0, i)] for i in range(8)],
        basis={"H": "6-31g"},
        symmetry=False,
        spin=2,
        verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASSCF(m, 4, 4)
    mc.max_cycle_macro = 20
    mc.fcisolver = DMRGSolver(mol)
    mc.fcisolver.file_path = None
    mc.fcisolver.parameters.set_orbital_optimization()
    emc_1 = mc.mc2step()[0]

    mol = gto.Mole()
    mol.build(
        atom=[["H", (0.0, 0.0, i)] for i in range(8)],
        basis={"H": "6-31g"},
        symmetry=False,
        spin=2,
        verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()
    mc = mcscf.CASSCF(m, 4, 4)
    mc.max_cycle_macro = 20
    emc_1ref = mc.mc2step()[0]
    print("Maquis-SCF = %.15g CASSCF = %.15g" % (emc_1, emc_1ref))
    print("Diff = %.15g" % (emc_1 - emc_1ref))
    assert abs(emc_1ref - emc_1) < 1e-9
    remove_checkpoint_and_results_file()


def test_dmrgscf_with_checkpoint():
    mol = gto.Mole()
    mol.build(
        atom=[["N", (0.0, 0.0, i)] for i in range(8)],
        basis={"N": "6-31g"},
        symmetry=False,
        spin=2,
        verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASSCF(m, 6, 6)
    mc.max_cycle_macro = 20
    mc.fcisolver = DMRGSolver(mol)
    mc.fcisolver.verbose = 1

    mc.fcisolver.file_path = os.path.dirname(os.path.abspath(__file__))
    mc.fcisolver.parameters.set_orbital_optimization()
    mc.fcisolver.parameters.set_entropies()
    emc_1 = mc.mc2step()[0]

    mol = gto.Mole()
    mol.build(
        atom=[["N", (0.0, 0.0, i)] for i in range(8)],
        basis={"N": "6-31g"},
        symmetry=False,
        spin=2,
        verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()
    mc = mcscf.CASSCF(m, 6, 6)
    mc.max_cycle_macro = 20
    emc_1ref = mc.mc2step()[0]
    print("Maquis-SCF = %.15g CASSCF = %.15g" % (emc_1, emc_1ref))
    print("Diff = %.15g" % (emc_1 - emc_1ref))
    assert abs(emc_1ref - emc_1) < 5e-9
    remove_checkpoint_and_results_file()


def test_dmrgscf_with_existing_checkpoint():
    mol = gto.Mole()
    mol.build(
        atom=[["N", (0.0, 0.0, i)] for i in range(8)],
        basis={"N": "6-31g"},
        symmetry=False,
        spin=2,
        verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASSCF(m, 6, 6)
    mc.max_cycle_macro = 20
    mc.fcisolver = DMRGSolver(mol)
    mc.fcisolver.verbose = 1

    mc.fcisolver.file_path = os.path.dirname(os.path.abspath(__file__))
    mc.fcisolver.parameters.set_orbital_optimization()
    mc.fcisolver.parameters.set_entropies()
    emc_1 = mc.mc2step()[0]

    mol = gto.Mole()
    mol.build(
        atom=[["N", (0.0, 0.0, i)] for i in range(8)],
        basis={"N": "6-31g"},
        symmetry=False,
        spin=2,
        verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()
    mc = mcscf.CASSCF(m, 6, 6)
    mc.max_cycle_macro = 20
    emc_1ref = mc.mc2step()[0]
    print("Maquis-SCF = %.15g CASSCF = %.15g" % (emc_1, emc_1ref))
    print("Diff = %.15g" % (emc_1 - emc_1ref))
    assert abs(emc_1ref - emc_1) < 5e-9
    remove_checkpoint_and_results_file()


def test_rdms():
    mol = gto.M(atom="Li 0 0 0; H 0 0 1.6", basis="sto-3g")
    mf = scf.RHF(mol).run()

    # full ci calculation
    norb, nelec = mol.nao_nr(), mol.nelec
    mc = fci.FCI(mf)
    e_fci, fcivec = mc.kernel()
    dm1_ref, dm2_ref = mc.make_rdm12(fcivec, norb, nelec)

    mc = mcscf.CASCI(mf, norb, nelec)
    mc.fcisolver = DMRGSolver(mol)
    mc.fcisolver.file_path = None
    emc_0 = mc.casci()[0]

    dm1_qc = mc.fcisolver._get_rdm1(norb)
    dm2_qc = mc.fcisolver._get_rdm2(norb)

    # Is this small enough
    assert (np.abs(dm1_ref - dm1_qc) < 1e-6).all()
    assert (np.abs(dm2_ref - dm2_qc) < 1e-6).all()


def remove_checkpoint_and_results_file():
    path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(path)
    try:
        shutil.rmtree("checkpoint_DMRGSCF.h5", ignore_errors=True)
    except:
        pass
    try:
        shutil.rmtree("checkpoint_gs.h5", ignore_errors=True)
    except:
        pass
    try:
        os.remove("results_DMRGSCF.h5")
    except:
        pass
    try:
        os.remove("results_file.h5")
    except:
        pass


if __name__ == "__main__":
    test_dmrgci()
    test_dmrgscf()
    test_dmrgscf_with_checkpoint()
    test_dmrgscf_with_existing_checkpoint()
    remove_checkpoint_and_results_file()
