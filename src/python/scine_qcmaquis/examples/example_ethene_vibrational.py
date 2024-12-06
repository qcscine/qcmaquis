#------------------------------------------------------
# Example of a vibrational DMRG calculation of 
# a one-body vibrational Hamiltonian of ethene
# with the python interface
#
# Before running, compile with "VIB=ON" and "QN=12"
#-----------------------------------------------------

import os
import shutil

import numpy as np
import pytest

import importlib
from scine_qcmaquis.maquis_dmrg import MaquisDmrg

dmrg = MaquisDmrg(model = "vibrational", vib_modes=12, num_basis="6,6,6,6,6,6,6,6,6,6,6,6")
dmrg.set_fcidump("FCIDUMP_ethene")
dmrg.init_mps(init_type="basis_state_generic", init_string="1,0,0,0,0,0,0,0,0,0,0,0")
dmrg.run_vibrational(bond_dimension=20)
