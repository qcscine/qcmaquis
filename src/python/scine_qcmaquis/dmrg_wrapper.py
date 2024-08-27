from enum import Enum
from typing import Any, List, Optional

# pylint: disable=import-error
from _dmrg import DmrgComplex, DmrgReal

from .integral_wrapper import ComplexTCIntegralMap, IntegralMapWrapper, TCIntegralMap
from .parameters_wrapper import ParametersWrapper

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
    """Wrapper for Dmrg Python Interface.

    Attributes
    ----------
    _dmrg : Union[DmrgReal, DmrgComplex]
        The DMRG interface.
    _run_flag : bool, default = False
        Flag indicating if dmrg run.
    _measure_flag : bool, default = False
        Flag indicating if dmrg did measurements.
    _run_option : RunOptions, default = RunOptions.OPTIMIZE
        Determine the dmrg algorithm.
    _feast_states : int
        Number of feast states.
    """

    __slots__ = ("_dmrg", "_run_flag", "_measure_flag", "_run_option", "_feast_states")

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

    def set_parameters(self, parameters: ParametersWrapper):
        """Initialize Dmrg object with DmrgParameters.

        Parameters
        ----------
        parameters : ParametersWrapper
            Wrapper around QCMaquis parameters
        """
        # in case new measurements appear in parameters
        self._measure_flag = False
        if "transcorrelated_hamiltonian" in parameters.get_parameters_dict():
            self._run_option = RunOptions.EVOLVE
            self._dmrg = DmrgReal(parameters.get_parameters())
        elif "feast_num_states" in parameters.get_parameters_dict():
            self._run_option = RunOptions.FEAST
            self._feast_states = parameters.get_parameters_dict()["feast_num_states"]
            self._dmrg = DmrgComplex(parameters.get_parameters())
        else:
            self._dmrg = DmrgReal(parameters.get_parameters())

    def get_fiedler(self, hf_occupations: Optional[List[List[int]]] = None, n_states: Optional[int] = None) -> str:
        """Evaluate Fiedler ordering.

        Parameters
        ----------
        hf_occupations : List[List[int]], optional
            The mean field occupation for each state
        n_states : int, optional
            number of states

        Return
        ------
        fiedler_string : str
            Orbital order from fiedler
        """
        if self._dmrg is None:
            raise ValueError("Set parameters before running dmrg!")

        if n_states is None or n_states == 0:
            n_states = 1
        if hf_occupations is None:
            hf_occupations = []

        # fiedler_calculator = self._dmrg
        zero_based_fiedler_string = self._dmrg.fiedler_order(n_states, hf_occupations, "fiedler")
        fiedler_string = ""
        for i in zero_based_fiedler_string.split(","):
            fiedler_string += str(int(i) + 1) + ","

        fiedler_string = fiedler_string[:-1]

        self._dmrg = None
        return fiedler_string

    def set_integrals(self, integral_map: IntegralMapWrapper):
        """Set integrals.

        Parameters
        ----------
        integral_map : IntegralMapWrapper
            wrapper around qcmaquis integral map
        """
        if self._dmrg is None:
            raise ValueError("Set parameters before running dmrg!")

        if isinstance(integral_map.get(), (ComplexTCIntegralMap, TCIntegralMap)):
            self._dmrg.update_tc_integrals(integral_map.get())
        else:
            self._dmrg.update_integrals(integral_map.get())

    def get_dmrg(self):
        """Get base dmrg object."""
        if self._dmrg is None:
            raise ValueError("Set parameters before running dmrg!")
        return self._dmrg

    def run(self):
        """Run Dmrg Calculation."""
        if self._dmrg is None:
            raise ValueError("Set parameters before running dmrg!")

        if self._run_option == RunOptions.OPTIMIZE:
            self._dmrg.optimize()
        elif self._run_option == RunOptions.EVOLVE:
            self._dmrg.evolve()
        elif self._run_option == RunOptions.FEAST:
            self._dmrg.runFEAST()

        self._run_flag = True

    # TODO: Make this for feast
    # def get_energy(self) -> Union[List[float], float]:
    def get_energy(self) -> float:
        """Get the energy from last calculation.

        Return
        ------
        energy : Union[List[float], float]
            one energy per state, if only one state is evaluated it's one float
        """
        if self._run_flag is False:
            raise ValueError("Run DMRG before asking for energies")
        # Feast gives you all energies at once
        # if self._run_option == RunOptions.FEAST:
        #     energies = []
        #     for i in range(self._feast_states):
        #         try:
        #             energies.append(self._dmrg.energyFEAST(i))
        #         # there are more states requested by feast than valid
        #         except RuntimeError:
        #             pass
        #     return energies
        return self._dmrg.energy()

    def measure(self) -> Any:
        """Measure set measurements."""
        if self._run_flag is False:
            raise ValueError("Run DMRG before asking for energies")
        if self._measure_flag is True:
            return
        self._measure_flag = True
        self._dmrg.measure()

    def get_ci_coefficient(self, det_string: str) -> float:
        """Getter for ci coeffs.

        Parameters
        ----------
        det_string : str
            QCMaquis compatible determinant string

        Return
        ------
        ci coeff : float
            The correponding CI coefficient.
        """
        return self._dmrg.getCICoefficient(det_string)

    def onerdm(self) -> Any:
        """Get 1rdm."""
        self.measure()
        return self._dmrg.onerdm()

    def twordm(self) -> Any:
        """Get 2rdm."""
        self.measure()
        return self._dmrg.twordm()
