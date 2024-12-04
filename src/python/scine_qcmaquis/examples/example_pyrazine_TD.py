#---------------------------------------------------------
# Exmaple TD-DMRG calculation with a 4-mode vibronic
# Hamiltonian for pyrazine. This script calls the analyze
# function at the end, generating plots of desired
# measurements.
#
# Before running, compile with "QD=ON"
#---------------------------------------------------------


import os
import shutil

import numpy as np
import pytest

import importlib
from scine_qcmaquis.entropy_builder import EntropyBuilder
from scine_qcmaquis.maquis_dmrg import MaquisDmrg

def _write_fcidump():
    tmp_fcidump = open("Pyrazine_VibHam_RedDim.txt", "w")
    tmp_fcidump.write("EL_ST 0 0\n")
    tmp_fcidump.write("-4114.23 0 0\n")
    tmp_fcidump.write("325.84799 1 1\n")
    tmp_fcidump.write("-325.84799 -1 -1\n")
    tmp_fcidump.write("559.34550 2 2\n")
    tmp_fcidump.write("-559.34550 -2 -2\n")
    tmp_fcidump.write("674.68277 3 3\n")
    tmp_fcidump.write("-674.68277 -3 -3\n")
    tmp_fcidump.write("518.21124 4 4\n")
    tmp_fcidump.write("-518.21124 -4 -4\n")
    tmp_fcidump.write("EL_ST 1 1\n")
    tmp_fcidump.write("4114.23 0 0\n")
    tmp_fcidump.write("325.84799 1 1\n")
    tmp_fcidump.write("-325.84799 -1 -1\n")
    tmp_fcidump.write("559.34550 2 2\n")
    tmp_fcidump.write("-559.34550 -2 -2\n")
    tmp_fcidump.write("674.68277 3 3\n")
    tmp_fcidump.write("-674.68277 -3 -3\n")
    tmp_fcidump.write("518.21124 4 4\n")
    tmp_fcidump.write("-518.21124 -4 -4\n")
    tmp_fcidump.write("EL_ST 0 0\n")
    tmp_fcidump.write("-791.22989 1 0\n")
    tmp_fcidump.write("-405.69687 2 0\n")
    tmp_fcidump.write("1171.11703 3 0\n")
    tmp_fcidump.write("EL_ST 1 1\n")
    tmp_fcidump.write("1092.881251 1 0\n")
    tmp_fcidump.write("-1379.877165 2 0\n")
    tmp_fcidump.write("302.457910 3 0\n")
    tmp_fcidump.write("EL_ST 0 0\n")
    tmp_fcidump.write("0.16131088 1 1\n")
    tmp_fcidump.write("8.71078783 1 2\n")
    tmp_fcidump.write("-16.45371035 1 3\n")
    tmp_fcidump.write("8.71078783 2 1\n")
    tmp_fcidump.write("-65.33090875 2 2\n")
    tmp_fcidump.write("38.23067993 2 3\n")
    tmp_fcidump.write("-16.45371035 3 1\n")
    tmp_fcidump.write("38.23067993 3 2\n")
    tmp_fcidump.write("-9.35603137 3 3\n")
    tmp_fcidump.write("-93.47965832 4 4\n")
    tmp_fcidump.write("EL_ST 1 1\n")
    tmp_fcidump.write("-73.96104114 1 1\n")
    tmp_fcidump.write("-24.03532198 1 2\n")
    tmp_fcidump.write("-15.24387871 1 3\n")
    tmp_fcidump.write("-24.03532198 2 1\n")
    tmp_fcidump.write("39.35985614 2 2\n")
    tmp_fcidump.write("9.27537593 2 3\n")
    tmp_fcidump.write("-15.24387871 3 1\n")
    tmp_fcidump.write("9.27537593 3 2\n")
    tmp_fcidump.write("1.77441974 3 3\n")
    tmp_fcidump.write("-93.47965832 4 4\n")
    tmp_fcidump.write("EL_ST 0 1\n")
    tmp_fcidump.write("1677.63321200 4 0\n")
    tmp_fcidump.write("50.65161814 4 1\n")
    tmp_fcidump.write("-44.44114904 4 2\n")
    tmp_fcidump.write("10.24324125 4 3\n")
    tmp_fcidump.write("50.65161814 1 4\n")
    tmp_fcidump.write("-44.44114904 2 4\n")
    tmp_fcidump.write("10.24324125 3 4\n")
    tmp_fcidump.write("EL_ST 1 0\n")
    tmp_fcidump.write("1677.63321200 4 0\n")
    tmp_fcidump.write("50.65161814 4 1\n")
    tmp_fcidump.write("-44.44114904 4 2\n")
    tmp_fcidump.write("10.24324125 4 3\n")
    tmp_fcidump.write("50.65161814 1 4\n")
    tmp_fcidump.write("-44.44114904 2 4\n")
    tmp_fcidump.write("10.24324125 3 4\n")
    tmp_fcidump.close()
_write_fcidump()

_write_fcidump()

dmrg = MaquisDmrg()
dmrg.set_parameter("L", 6)
dmrg.set_parameter("symmetry", "u1")
dmrg.set_parameter("LATTICE", "vibronic lattice")
dmrg.set_parameter("MODEL", "vibronic")
dmrg.set_parameter("Nmax", 6)
dmrg.set_parameter("vibronic_num_elestates", 2)
dmrg.set_parameter("vibronic_num_vibmodes", 4)
dmrg.set_parameter("init_type", "coherent")
dmrg.set_parameter("init_coeffs", "0.8,0.2")
dmrg.set_parameter("init_basis_state", "0,1,0,0,0,0|1,0,0,0,0,0")
dmrg.set_parameter("nsweeps", 40)
dmrg.set_parameter("max_bond_dimension", 20)
dmrg.set_parameter("time_step", 1)
dmrg.set_parameter("time_units", "fs")
dmrg.set_parameter("propagator_maxiter", 40)
dmrg.set_parameter("propagator_max_accuracy", 1.0E-10)
dmrg.set_parameter("max_bond_dimension", 20)
dmrg.set_parameter("hamiltonian_units", "cm-1")
dmrg.set_parameter("resultfile", "res.h5")
dmrg.set_parameter("MEASURE[Autocorrelation]", 1)
dmrg.set_parameter("MEASURE[Population]", 1)
dmrg.set_parameter("ALWAYS_MEASURE", "Autocorrelation,PopulationState0,PopulationState1")
dmrg.set_fcidump("Pyrazine_VibHam_RedDim.txt")
dmrg.evolve()
print(dmrg._dmrg._run_option)
print("time evolution finished")
measurements = ["autocorrelation", "spectrum", "population"]
dmrg.analyze(measurements)
