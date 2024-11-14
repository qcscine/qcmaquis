import os

import numpy as np
import pytest

from scine_qcmaquis.parameters_wrapper import ParametersWrapper


def test_parameters_wrapper_hf_occupation():
    parameters = ParametersWrapper()
    norbs = 10
    nelec = 9
    spin = 1
    parameters._make_hf_occupation(norbs, nelec, spin)
    hf_occ = "4,4,4,4,3,1,1,1,1,1"
    assert hf_occ == parameters._parameter_dict["hf_occ"]

    parameters = ParametersWrapper()
    norbs = 10
    nelec = 9
    spin = 3
    parameters._make_hf_occupation(norbs, nelec, spin)
    hf_occ = "4,4,4,3,3,3,1,1,1,1"
    assert hf_occ == parameters._parameter_dict["hf_occ"]


def test_parameters_wrapper_set_system():
    parameters = ParametersWrapper()
    norbs = 10
    nelec = 9
    spin = 3
    parameters.set_system(norbs, nelec, spin)

    assert parameters._parameter_dict["u1_total_charge1"] == 6
    assert parameters._parameter_dict["u1_total_charge2"] == 3
    assert parameters._parameter_dict["spin"] == spin
    assert parameters._parameter_dict["nelec"] == nelec
    assert parameters._parameter_dict["L"] == norbs
    assert parameters._parameter_dict["hf_occ"] == "4,4,4,3,3,3,1,1,1,1"
    assert parameters._parameter_dict["site_types"] == "0,0,0,0,0,0,0,0,0,0"


def test_parameters_wrapper_convenience_functions():
    parameters = ParametersWrapper()
    assert parameters._parameter_dict["MODEL"] == "quantum_chemistry"
    assert parameters._parameter_dict["init_type"] == "default"
    assert parameters._parameter_dict["irrep"] == 0
    assert parameters._parameter_dict["nsweeps"] == 100
    assert parameters._parameter_dict["max_bond_dimension"] == 250
    assert parameters._parameter_dict["optimization"] == "twosite"
    assert parameters._parameter_dict["conv_thresh"] == 1e-6
    assert parameters._parameter_dict["symmetry"] == "su2u1pg"
    assert parameters._parameter_dict["CONSERVED_QUANTUMNUMBERS"] == "Nup,Ndown"
    assert parameters._parameter_dict["lattice_library"] == "coded"
    assert parameters._parameter_dict["model_library"] == "coded"
    assert parameters._parameter_dict["LATTICE"] == "orbitals"

    parameters.set_integral_file("abc_fcidump_yzx")
    assert parameters._parameter_dict["integral_file"] == "abc_fcidump_yzx"
    assert "integrals" not in parameters._parameter_dict

    parameters.set_orbital_optimization()
    assert parameters._parameter_dict["MEASURE[1rdm]"] is True
    assert parameters._parameter_dict["MEASURE[2rdm]"] is True

    parameters.set_entropies()
    assert parameters._parameter_dict["MEASURE[ChemEntropy]"] is True

    parameters.set_result_path("/abc/cde/results.h5")
    assert parameters._parameter_dict["resultfile"] == "/abc/cde/results.h5"
    assert parameters._results_path == "/abc/cde/results.h5"

    parameters.set_checkpoint_path("/abc/cde/checkpoint")
    assert parameters._parameter_dict["chkpfile"] == "/abc/cde/checkpoint.h5"
    assert parameters._checkpoint_path == "/abc/cde/checkpoint.h5"


def test_parameters_wrapper_maquis_dump_files():
    parameters = ParametersWrapper()

    parameters.set_result_path("/abc/cde/results.h5")
    assert parameters._parameter_dict["resultfile"] == "/abc/cde/results.h5"
    assert parameters._results_path == "/abc/cde/results.h5"

    parameters.set_checkpoint_path("/abc/cde/checkpoint")
    assert parameters._parameter_dict["chkpfile"] == "/abc/cde/checkpoint.h5"
    assert parameters._checkpoint_path == "/abc/cde/checkpoint.h5"

    parameters.set_storage_dir("/abc/cde/storagedir")
    assert parameters._parameter_dict["storagedir"] == "/abc/cde/storagedir"
    assert parameters._storage_dir == "/abc/cde/storagedir"


if __name__ == "__main__":
    test_parameters_wrapper_hf_occupation()
    test_parameters_wrapper_set_system()
    test_parameters_wrapper_convenience_functions()
    test_parameters_wrapper_maquis_dump_files()
