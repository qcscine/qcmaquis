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

dmrg = MaquisDmrg()
dmrg.set_parameter("nsweeps", 10)
dmrg.set_parameter("ngrowsweeps", 2)
dmrg.set_parameter("nmainsweeps", 3)
dmrg.set_parameter("max_bond_dimension", 25)
dmrg.set_parameter("alpha_initial", 1.0E-9)
dmrg.set_parameter("alpha_main", 1.0E-10)
dmrg.set_parameter("alpha_final", 0)
dmrg.set_parameter("truncation_initial", 1.0E-8)
dmrg.set_parameter("truncation_final", 1.0E-10)
dmrg.set_parameter("eigensolver", "IETL_JCD")
dmrg.set_parameter("conv_thresh", 1.0E-3)
dmrg.set_parameter("optimization", "singlesite")
dmrg.set_parameter("symmetry", "nu1")
dmrg.set_parameter("MODEL", "nmode")
dmrg.set_parameter("model_library", "coded")
dmrg.set_parameter("lattice library", "coded")
dmrg.set_parameter("LATTICE", "nmode lattice")
dmrg.set_parameter("L", 72)
dmrg.set_parameter("nmode_num_modes", 12)
dmrg.set_parameter("nmode_num_basis", "6,6,6,6,6,6,6,6,6,6,6,6")
dmrg.set_parameter("init_type", "basis_state_generic")
dmrg.set_parameter("init_basis_state", "0,0,0,0,0,0,0,0,0,0,0,0")
dmrg.set_parameter("integral_cutoff", 1.0E-8)
dmrg.set_fcidump("FCIDUMP_ethene")
dmrg.run_vibrational()
