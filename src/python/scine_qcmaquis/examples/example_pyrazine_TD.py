#---------------------------------------------------------
# Example TD-DMRG calculation with a 4-mode vibronic
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

dmrg = MaquisDmrg(model="vibronic", elec_states=2, vib_modes=4)
dmrg.set_fcidump("Pyrazine_VibHam_RedDim.txt")
dmrg.init_mps(init_type="coherent", init_string="0,1,0,0,0,0|1,0,0,0,0,0", init_coeffs="0.8,0.2")
dmrg.measure_population()
dmrg.evolve(t_step=1, n_steps=40, t_units="fs")
print(dmrg._dmrg._run_option)
print("time evolution finished")
measurements = ["autocorrelation", "spectrum", "population"]
dmrg.analyze_results(measurements)
