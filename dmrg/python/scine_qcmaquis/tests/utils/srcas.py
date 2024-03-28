import os
import shutil

import pytest
from pyscf import gto, mcscf, scf

from scine_qcmaquis import MaquisDmrg
from scine_qcmaquis.pyscf_interface.pyscf_interface import QcMaquis
from scine_qcmaquis.utils.srcas import Srcas

# TODO: python bindings for SRCAS not testable atm, because results have to be modified in cpp code


def test_srcas_with_interface():
    mol = gto.Mole()
    mol.build(
        atom="""N 0 0 0
                N 0 0 2""",
        basis='6-31g', symmetry=False, spin=2, verbose=1,
    )
    m = scf.RHF(mol)
    m.scf()
    mc = mcscf.CASCI(m, 6, 6)
    mc.fcisolver = QcMaquis(mol)
    mc.fcisolver.file_path = os.path.dirname(os.path.abspath(__file__))
    mc.fcisolver.parameters.set_orbital_optimization()
    mc.casci()
    dmrg = mc.fcisolver.dmrg

    srcas = Srcas(dmrg)
    srcas.run()
    srcas.print()
    # results = srcas.results()
    # for i in results:
    #     print(i)

    path = os.path.dirname(os.path.abspath(__file__))
    os.chdir(path)
    shutil.rmtree("checkpoint_gs.h5", ignore_errors=True)
    os.remove("results_file.h5")


def test_srcas():

    dmrg = MaquisDmrg()
    dmrg.set_parameter("symmetry", "2u1pg")
    dmrg.init_dmrg("checkpoint_n2_triplet.2.2.h5", 6, 6, 2)

    srcas = Srcas(dmrg)
    # srcas.run(dmrg, "4,4,3,3,1,1")
    srcas.run()
    srcas.print()
    # results = srcas.results()
    # for i in results:
    #     print(i, results[i])


if __name__ == "__main__":
    test_srcas()
    test_srcas_with_interface()
    # setup()
    # try:
    # test_srcas(dmrg)
    # except:
    #     pass
    # remove_checkpoint_and_results_file()
