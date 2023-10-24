"""Blub blub."""

from enum import Enum
from typing import List

# pylint: disable=import-error
from _dmrg import DmrgComplex, DmrgReal

from dmrg.python.integral_wrapper import ComplexTCIntegralMap, IntegralMapWrapper, TCIntegralMap
from dmrg.python.parameters_wrapper import ParametersWrapper

# pylint: enable=import-error


class RunOptions(Enum):
    """Set run option for qcmaquis."""
    OPTIMIZE = "optimize"
    """For conventional DMRG."""
    EVOLVE = "evolve"
    """For time dependent DMRG."""
    FEAST = "feast"
    """For FEAST calculations."""


class DmrgWrapper:
    """Wrapper for Dmrg Python Interface."""

    def __init__(self):
        """Construct Wrapper."""
        self._dmrg = None
        """The DMRG interface."""
        self._run_flag = False
        """Flag indicating if dmrg run."""
        self._measure_flag = False
        """Flag indicating if dmrg did measurements."""
        self._run_option = RunOptions.OPTIMIZE
        """Determine the dmrg algorithm."""
        self._feast_states = 0
        """Number of feast states."""

    # def _get_dmrg(self):
    #     """Get dmrg object."""
    #     return self._dmrg

    def set_parameters(self, parameters: ParametersWrapper):
        """Initialize Dmrg object with DmrgParameters."""
        # in case new measurements appear in parameters
        self._measure_flag = False
        # pylint: disable=protected-access
        if "transcorrelated_hamiltonian" in parameters._parameter_dict:
            self._run_option = RunOptions.EVOLVE
            self._dmrg = DmrgComplex(parameters._get_parameters())
        elif "feast_num_states" in parameters._parameter_dict:
            self._run_option = RunOptions.FEAST
            self._feast_states = parameters._parameter_dict["feast_num_states"]
            self._dmrg = DmrgComplex(parameters._get_parameters())
        else:
            self._dmrg = DmrgReal(parameters._get_parameters())
        # pylint: enable=protected-access

    def get_fiedler(self, hf_occupations: List[List[int]] = None, n_states: int = None):
        if self._dmrg is None:
            raise ValueError("Set parameters before running dmrg!")

        if n_states is None or n_states == 0:
            n_states = 1
        if hf_occupations is None:
            hf_occupations = []

        # fiedler_calculator = self._dmrg
        zero_based_fiedler_string = self._dmrg.fiedler_order(n_states, hf_occupations, "fiedler")
        print("huihuhuh")
        fiedler_string = ""
        for i in zero_based_fiedler_string.split(","):
            fiedler_string += str(int(i) + 1) + ","

        fiedler_string = fiedler_string[:-1]
        print(fiedler_string)

        self._dmrg = None
        return fiedler_string

    def set_integrals(self, integral_map: IntegralMapWrapper):
        """Set integrals."""
        if self._dmrg is None:
            raise ValueError("Set parameters before running dmrg!")
        if type(integral_map.get()) is ComplexTCIntegralMap:
            print("hihi")
            self._dmrg.update_tc_integrals(integral_map.get())
        else:
            print("haha")
            self._dmrg.update_integrals(integral_map.get())

    def run(self):
        """Run Dmrg Calculation"""
        print("run")
        if self._dmrg is None:
            raise ValueError("Set parameters before running dmrg!")
        if self._run_option == RunOptions.OPTIMIZE:
            print("optimize")
            self._dmrg.optimize()
        elif self._run_option == RunOptions.EVOLVE:
            print("evolve")
            self._dmrg.evolve()
        elif self._run_option == RunOptions.FEAST:
            # self._dmrg.optimize()
            self._dmrg.runFEAST()
        self._run_flag = True

    def get_energy(self):
        """Get the energy from last calculation."""
        if self._run_flag is False:
            raise ValueError("Run DMRG before asking for energies")
        # Feast gives you all energies at once
        if self._run_option == RunOptions.FEAST:
            energies = []
            for i in range(self._feast_states):
                try:
                    energies.append(self._dmrg.energyFEAST(i))
                # there are more states requested by feast than valid
                except RuntimeError:
                    pass
            return energies
        return self._dmrg.energy()

    def measure(self):
        """Measure set measurements."""
        if self._run_flag is False:
            raise ValueError("Run DMRG before asking for energies")
        if self._measure_flag is True:
            return
        self._measure_flag = True
        self._dmrg.measure()

    def get_ci_coefficient(self, det_string: str) -> float:
        # print(det_string)
        return self._dmrg.getCICoefficient(det_string)

    def onerdm(self):
        """Get 1rdm."""
        self.measure()
        return self._dmrg.onerdm()

    def twordm(self):
        """Get 1rdm."""
        self.measure()
        return self._dmrg.twordm()
