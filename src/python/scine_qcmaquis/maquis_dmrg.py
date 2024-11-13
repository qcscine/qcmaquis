from typing import Any, List, Optional, Tuple, Union

import numpy as np

# pylint: disable=import-error
from .dmrg_wrapper import DmrgWrapper
from .entropy_builder import EntropyBuilder
from .integral_wrapper import ComplexTCIntegralMap, IntegralMap, IntegralMapWrapper, IntegralType, TCIntegralMap
from .parameters_wrapper import ExcitedStates, ParametersWrapper

# pylint: enable=import-error


class MaquisDmrg:
    """Python Interface to QCMaquis.

    Attributes
    ----------
    _dmrg : DmrgWrapper
        handler for executing dmrg calculation
    _parameters : ParametersWrapper
        handler for parameters
    _integral_map : IntegralMapWrapper
        handler for integrals
    _transcorrelated : bool, default = False
        flag to enable transcorrelated calculations
    _orbital_optimization : bool, default = False
        flag to enable measurements of the 1 and 2 rdm
    _excited_states : bool, default = False
        flag to enable excited state calculations
    _energy : Union[float, List[float]]
        the energy of one or more states
    """
    __slots__ = (
        "_dmrg",
        "_parameters",
        "_integral_map",
        "_transcorrelated",
        "_orbital_optimization",
        "_excited_states",
        "_energy",
        "_entropy_builder",
    )

    def __init__(self) -> None:
        """Construct Wrapper.

        Note
        ----
        All attributes should be modified by the corresponding functions, to ensure the expected behavior.
        """
        self._dmrg = DmrgWrapper()
        """Handler for calculations."""
        self._parameters = ParametersWrapper()
        """Handler for parameters."""
        self._integral_map = IntegralMapWrapper()
        """Handler for integrals."""
        self._transcorrelated = False
        """Flag for transcorrelation."""
        self._orbital_optimization = False
        """Flag for orbital optimization."""
        self._excited_states = False
        """Flag for excited states."""
        self._energy: Union[float, List[float]] = 0.0
        """Final energy of the system."""
        self._entropy_builder: Optional[EntropyBuilder] = None
        """Assembly s1, s2 and mut inf from qcmaquis"""

    def replace_parameters(self, parameters_wrapper: ParametersWrapper):
        """Replace existing parameters wrapper with new parameters.

        Parameters
        ----------
        parameters_wrapper : ParametersWrapper
            parameter object
        """
        self._parameters = parameters_wrapper

    def set_parameter(self, parameter_name: str, parameter_value: Any):
        """Set any parameter in DmrgParameters.

        Parameters
        ----------
        parameter_name : str
            name of the parameter to set
        parameter_value : Any
            the value of the parameter to set

        Note
        ----
        The parameter names and value types are the same as defined by QcMaquis.
        """
        self._parameters.set(parameter_name, parameter_value)

    def set_bond_dimension(self, bond_dimension: int):
        """Set bond dimension.

        Parameters
        ----------
        bond_dimension : int
            the bond dimension
        """
        self._parameters.set("max_bond_dimension", bond_dimension)

    def set_excited_states(
        self,
        n_excited_states: int,
        method: ExcitedStates = ExcitedStates.ORTHO,
        feast_window: Optional[List[float]] = None
    ):
        """Enable Excited States.

        This function acts as general interface to excited states, independent of the requested method.

        Parameters
        ----------
        method : ExcitedStates
            The enum value of the requested method
        n_excited_states : int
            number of excited states
        feast_window: List[float]
            energy window for feast calculations (only used for FEAST)
        """
        self._excited_states = True
        if method == ExcitedStates.ORTHO and n_excited_states is None:
            raise ValueError("Orthogonal excited states require 'n_excited_states'")
        if method == ExcitedStates.FEAST and feast_window is None:
            raise ValueError("Feast excited states require to set 'feast_window'")

        if method == ExcitedStates.ORTHO:
            self._parameters.set_excited_states_ortho(n_excited_states)
        elif method == ExcitedStates.FEAST:
            raise ValueError("Feast is not implemented yet")

    def set_entropies(self):
        """Enable entropy measurements."""
        self._parameters.set_entropies()

    def set_orbital_optimization(self):
        """Enable orbital optimization.

        Set flags and parameters required for orbital optimization
        with QcMaquis as FCI solver.
        """
        self._orbital_optimization = True
        self._parameters.set_orbital_optimization()

    def set_transcorrelation(self):
        """Enable transcorrelation.

        Set flags and parameters required for transcorrelated DMRG
        calculations with QcMaquis.
        """
        self._transcorrelated = True
        self._parameters.set_transcorrelation_values()
        self._integral_map.set_type(IntegralType.TRANSCORRLEATED)

    def get_entropies(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Getter for entropies.

        Return
        ------
        s1 : np.ndarray
            The s1 entropies
        s2 : np.ndarray
            The s2 entropies
        ipq : np.ndarray
            The mutual information
        """
        self._dmrg.measure()
        if self._entropy_builder is None:
            raise AssertionError("Initialize entropy builder before extracting entropies")
        self._entropy_builder.make_diagnostics(self._dmrg.get_dmrg())
        return (
            self._entropy_builder.s1_entropy,
            self._entropy_builder.s2_entropy,
            self._entropy_builder.mutual_information
        )

    def get_reduced_density_matrices(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get 1 and 2 RDM.

        Run measurements after a succesful DMRG calculation.
        This is required for DMRGSCF calculations.

        Return
        ------
        onerdm : np.ndarray
            the one particle reduced density matrix
        twordm : np.ndarray
            the two particle reduced density matrix
        """
        self._dmrg.measure()
        onerdm = self._dmrg.onerdm()
        twordm = self._dmrg.twordm()
        return onerdm, twordm

    def get_energy(self) -> Union[float, List[float]]:
        """Get energy.

        The type of energy depends on the type of calculation.
        If multiple states were requested, e.g. excited states,
        it will return a list with these energies.

        Return
        ------
        energy : Union[float, List[float]]
            the final energies
        """
        # try:
        self._dmrg._run_flag = True
        self._energy = self._dmrg.get_energy()
        # except:

        if self._energy == 0.0:
            raise RuntimeWarning("Run DMRG before requesting energies")

        return self._energy

    def set_feast(self, window: Tuple[float, float], n_states: int):
        """Enable FEAST calculations.

        Parameters
        ----------
        window : Tuple[float, float]
            the energy window with [e_min, e_max]
        n_states : int
            number of states in dmrg feast
        """
        self._parameters.set_excited_states_feast(window, n_states)

    def set_fcidump(self, fcidump: str):
        """Read integrals from fcidump.

        Parameters
        ----------
        fcidump : str
            Path to the fcidump
        """
        self._parameters.set_integral_file(fcidump)

    def run(
        self,
        n_orbitals: int,
        n_electrons: int,
        spin: int = 0,
        n_states: Optional[int] = None,
        fiedler: bool = False
    ):
        """Run Dmrg.

        Note
        ----
        This function is responsible for ALL dmrg calculations, independent
        of the type of calculation, e.g. transcorrelation, ortho, feast,
        groundstate, ...

        Parameteters
        ------------
        n_orbitals : int
            the number of orbitals in the active space
        n_electrons : int
            the number of electrons in the active space
        spin : int, default = 0
            the total spin in the active space, e.g. 2S
        n_states : int, default = None
            number of orthogonal states
        fiedler : int, default = False
            Flag to enable Fiedler ordering
        """

        # excited states
        if n_states is not None:
            # start with ground state
            self._parameters.set_excited_states_ortho(0)

        self._parameters.set_system(n_orbitals, n_electrons, spin)

        if fiedler is True and "orbital_order" not in self._parameters.get_parameters_dict():
            # don't dump anything for fiedler
            try:
                tmp_chkpfile = self._parameters.get_parameters_dict()["chkpfile"]
            except KeyError:
                tmp_chkpfile = ""
            try:
                tmp_result_file = self._parameters.get_parameters_dict()["resultfile"]
            except KeyError:
                tmp_result_file = ""

            self._parameters.erase("chkpfile", verbose=False)
            self._parameters.erase("resultfile", verbose=False)

            fiedler_orderer = DmrgWrapper()
            fiedler_orderer.set_parameters(self._parameters)
            if "integral_file" not in self._parameters.get_parameters_dict():
                fiedler_orderer.set_integrals(self._integral_map)
            orbital_order = fiedler_orderer.get_fiedler()
            self._parameters.set("orbital_order", orbital_order)
            self._entropy_builder = EntropyBuilder(n_orbitals, orbital_order)

            # but enable dumping for real calc
            if tmp_chkpfile:
                self._parameters.set_checkpoint_path(tmp_chkpfile)
            if tmp_result_file:
                self._parameters.set_result_path(tmp_result_file)

        elif "orbital_order" in self._parameters.get_parameters_dict():
            orbital_order = self._parameters.get("orbital_order")
            self._entropy_builder = EntropyBuilder(n_orbitals, orbital_order)

        else:
            self._entropy_builder = EntropyBuilder(n_orbitals)

        self._dmrg.set_parameters(self._parameters)

        if "integral_file" not in self._parameters.get_parameters_dict():
            self._dmrg.set_integrals(self._integral_map)

        for i in self._parameters.get_parameters_dict():
            print(i, self._parameters.get_parameters_dict()[i])

        self._dmrg.run()
        # excited states
        if n_states is not None:
            self._energy = []
            self._energy.append(self._dmrg.get_energy())
            for state in range(1, n_states):
                self._parameters.set_system(n_orbitals, n_electrons, spin)
                self._parameters.set_excited_states_ortho(state)
                self._dmrg.set_parameters(self._parameters)
                if "integral_file" not in self._parameters.get_parameters_dict():
                    self._dmrg.set_integrals(self._integral_map)
                self._dmrg.run()
                self._energy.append(self._dmrg.get_energy())
        else:
            self._energy = self._dmrg.get_energy()

    def update_integrals(self, integral_map: Union[IntegralMap, TCIntegralMap, ComplexTCIntegralMap]):
        """Update integrals.

        Parameters
        ----------
        integral_map : Union[IntegralMap, TCIntegralMap, ComplexTCIntegralMap]
            the integral map in corresponding symmetry, e.g.
            8-fold for conventional (4 indices)
            2-fold for transcorrelated (6 indices)
        """
        self._integral_map.set(integral_map)

    def set_integrals(self, core_value: float, one_body: np.ndarray, two_body: np.ndarray, norb: int):
        """Set integrals from PySCF.

        Parameters
        ----------
        core_value : float
            the core integral
        one_body : np.ndarray
            the one body integrals
        two_body : np.ndarray
            the two body integrals
        norb : int
            number of orbitals, e.g. the shape of one and two body integrals
        """
        self._integral_map.fill_from_pyscf(core_value, one_body, two_body, norb)

    def init_dmrg(self, checkpoint: str, norb: int, nelec: int, spin: int):
        """Initialize DMRG object.

        Parameters
        ----------
        checkpoint : str
            Path to the checkpoint file
        norb : int
            number of orbitals
        nelec : int
            number of electrons
        spin : int
            Total spin of the system (2S)
        """
        self._parameters.set_system(norb, nelec, spin)
        print(checkpoint)
        self._parameters.set_checkpoint_path(checkpoint)
        self._dmrg.set_parameters(self._parameters)
        self._dmrg.set_integrals(self._integral_map)
        # self._dmrg.run()

    def get_ci_coefficient(self, determinant_string: str) -> float:
        """Get CI coefficient from determinant string.

        Parameters
        ----------
        determinant_string : str
            a qcmaquis compatible determinant string

        Returns
        -------
        ci_coeff : float
            The corresponding ci coefficient
        """
        return self._dmrg.get_ci_coefficient(determinant_string)
