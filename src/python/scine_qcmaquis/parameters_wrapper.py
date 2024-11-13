from enum import Enum
from typing import Any, Dict, Tuple

# pylint: disable=import-error
from _dmrg import DmrgParameters

# pylint: enable=import-error


class ExcitedStates(Enum):
    """Store excited state methods.

    Attributes
    ----------
    NO
    ORTHO
    FEAST
    """
    NO = "no"
    """Disable excited states"""
    ORTHO = "ortho"
    """Optimize excited states with ortho"""
    FEAST = "feast"
    """Optimize excited states with FEAST"""


class ParametersWrapper:
    """Wrapper for DmrgParameters Python Interface.

    Attributes
    ----------
    _parameters : DmrgParameters
        dmrg parameters map
    _parameter_dict: Dict[str, Any]
        for reference, same as dmrg parameters but native python
    _checkpoint_path : str, default = "checkpoint_gs"
        name and path to store the checkpoint file.
    _results_path : str, default = "results_file.h5"
        name and path to store the results file.
    _excited_state_name : str
        name for excited states checkpoint files.
    """

    __slots__ = (
        "_parameters",
        "_parameter_dict",
        "_checkpoint_path",
        "_results_path",
        "_excited_state_name",
        "_storage_dir",
    )

    def __init__(self, set_defaults: bool = True) -> None:
        """Constructor.

        Parameters
        ----------
        set_defaults : bool, default = True
            disable default parameters and use empty parameters
        """
        self._parameters = DmrgParameters()
        """dmrg parameters map"""
        self._parameter_dict: Dict[str, Any] = {}
        """for reference, same as dmrg parameters but native python"""
        self._checkpoint_path: str = "checkpoint_gs.h5"
        """name and path to store the checkpoint file."""
        self._results_path: str = "results_file.h5"
        """name and path to store the results file."""
        self._storage_dir: str = "qcmaquis_storage_dir"
        """name and path to store the results file."""
        self._excited_state_name: str = self._checkpoint_path[:-3] + "ex0"
        """name for excited states checkpoint files."""

        # set default parameters
        if set_defaults:
            self._set_defaults()

    def set_storage_dir(self, path: str):
        """Set path and name of storagedir.

        Parameters
        ----------
        path : str
            the path
        """
        self._storage_dir = path
        self.set("storagedir", path)

    def get_result_path(self) -> str:
        """Return current result path.

        Returns
        -------
        self.result_path : str
            Path to the result file
        """
        return self._results_path

    def set_result_path(self, path: str):
        """Set path and name to checkpoint file.

        Parameters
        ----------
        path : str
            the path
        """
        if not path.endswith(".h5"):
            path += ".h5"
            # raise ValueError("result_path has to point to <.h5> file")
        self._results_path = path
        self.set("resultfile", path, verbose=False)

    def get_checkpoint_path(self) -> str:
        """Return current checkpoint path.

        Returns
        -------
        self.checkpoint_path : str
            Path to the checkpoint file
        """
        return self._checkpoint_path

    def set_checkpoint_path(self, path: str):
        """Set path and name to checkpoint file.

        Parameters
        ----------
        path : str
            the path
        """
        if not path.endswith(".h5"):
            path += ".h5"
        self._checkpoint_path = path
        self.set("chkpfile", self._checkpoint_path, verbose=False)

    def _set_defaults(self):
        """Set default parameters.

        Note
        ----
        "MODEL" = "quantum_chemistry"
        "init_type" = "default"
        "irrep" = 0
        "nsweeps" = 100
        "max_bond_dimension" = 250
        "optimization" = "twosite"
        "conv_thresh" = 1e-6
        "symmetry" = "su2u1pg"
        "CONSERVED_QUANTUMNUMBERS" = "Nup,Ndown"
        "lattice_library" = "coded"
        "model_library" = "coded"
        "LATTICE" = "orbitals"

        "integrals" are set one arbitrary entry, since QcMaquis requires this keyword to be set.
        However it's not used anywhere, since we update the integral map later anyways.
        """
        self.set("MODEL", "quantum_chemistry")
        self.set("init_type", "default")
        self.set("irrep", 0)
        self.set("nsweeps", 100)
        self.set("max_bond_dimension", 250)
        self.set("optimization", "twosite")
        self.set("conv_thresh", 1e-6)
        self.set("symmetry", "su2u1pg")
        self.set("CONSERVED_QUANTUMNUMBERS", "Nup,Ndown")
        self.set("lattice_library", "coded")
        self.set("model_library", "coded")
        self.set("LATTICE", "orbitals")

        # dummy parameters for qcmaquis
        # dmrg requires integrals, integrals_binary or integrals_file to be set
        # already in constructor
        # Here we update the integrals later anyways with a new integral map
        self.set("integrals", "   0.00000000000              1     1     1     1")

    def _make_site_types(self, n_orbitals: int):
        """Generate the string for site types.

        Parameters
        ----------
        n_orbitals : int
            the number of orbitals

        Note
        ----
        Only implemented for C1 symmetry yet.
        """
        if self._parameter_dict["irrep"] != 0:
            raise ValueError("irrep was changed, and this function is only implemented for irrep = 0")

        site_types = ""
        for _ in range(n_orbitals):
            site_types += "0,"
        self.set("site_types", site_types[:-1])

    def set_transcorrelation_values(self):
        """Defaults to transcorrelated Hamiltonian.

        Note
        ----
        "simulation_type" = "TD"
        "transcorrelated_hamiltonian" = "yes"
        "propagator_accuracy" = 1.0E-10
        "propagator_maxiter" = 10
        "hamiltonian_units" = "Hartree"
        "time_units" = "fs"
        "imaginary_time" = "yes"
        "TD_backpropagation" = "no"
        "symmetry" = "2u1"
        "nsweeps" = 10
        "time_step" = 0.2
        """

        self.set("simulation_type", "evolve")
        self.set("transcorrelated_hamiltonian", "yes")
        self.set("propagator_accuracy", 1.0E-10)
        self.set("propagator_maxiter", 10)
        self.set("hamiltonian_units", "Hartree")
        self.set("time_units", "fs")
        self.set("imaginary_time", "yes")
        self.set("TD_backpropagation", "no")
        self.set("symmetry", "2u1")
        self.set("nsweeps", 10)
        self.set("time_step", 0.2)
        # self.set("delta_t", 0.02)
        # self.set("COMPLEX", True)

    def set_excited_states_feast(self, energy_window: Tuple[float, float], n_states: int = 2):
        """Set FEAST parameters.

        Parameters
        ----------
        energy_window: Tuple[float, float]
            the energy window
        n_states: int, default = 2
            the number of states within this window

        Note
        ----
        The number of states should always be larger than the required number of states.
        "feast_emin" = energy_window[0]
        "feast_emax" = energy_window[1]
        "feast_num_states" = n_states
        "feast_num_points" = 8
        """
        if energy_window[0] > energy_window[1]:
            raise ValueError("The energy window must be set as (<e_min>, <e_max>)")

        self.set("COMPLEX", True)
        self.set("feast_emin", energy_window[0])
        self.set("feast_emax", energy_window[1])
        self.set("feast_num_states", n_states)
        self.set("feast_max_iter", 0)
        self.set("feast_num_points", 8)

    def set_excited_states_ortho(self, n_excited_states: int):
        """Set Ortho States.

        Parameters
        ----------
        n_excited_states : int
            the number of states to be evaluated with ortho

        Note
        ----
        this funciton update the _checkpoint_path to the checkpoint file of the current state.
        "chkpfile" = self._checkpoint_path
        "chkpfile" = self._excited_state_name[:-1] + n_excited_states
        "n_ortho_states" = n_excited_states
        "ortho_states" = self._checkpoint_path
        """
        if n_excited_states == 0:
            self.set("chkpfile", self._checkpoint_path)
        else:
            self.set("chkpfile", self._excited_state_name[:-1] + str(n_excited_states))
            self.set("n_ortho_states", n_excited_states, verbose=False)
            self.set("ortho_states", self._checkpoint_path, verbose=False)
            self._checkpoint_path =\
                self._checkpoint_path + "," + self._excited_state_name[:-1] + str(n_excited_states)

    def set_orbital_optimization(self):
        """Set Orbital optimization parameters.

        Note
        ----
        Even though the parameters are set to measure the 1 and 2rdm,
        qcmaquis still has to perform the measurement.
        "MEASURE[1rdm]" = True
        "MEASURE[2rdm]" = True
        """
        self.set("MEASURE[1rdm]", True)
        self.set("MEASURE[2rdm]", True)

    def set_entropies(self):
        """Set Entropy parameters.

        Note
        ----
        Even though the parameters are set to measure the entropies
        qcmaquis still has to perform the measurement.
        "MEASURE[ChemEntropy]" = True
        """
        self.set("MEASURE[ChemEntropy]", True)

    def set_system(self, n_orbitals: int, n_electrons: int, spin: int = 0):
        """Set system specifics.

        Parameters
        ----------
        n_orbitals : int
            number of orbitals
        n_electrons : int
            number of electrons
        spin : int, default = 0
            spin of the system
        """
        # number of "3"s in det string
        self.set("u1_total_charge1", int((n_electrons - spin) / 2) + int(spin))
        # number of "2"s in det string
        self.set("u1_total_charge2", int((n_electrons - spin) / 2))
        self.set("spin", spin)
        self.set("nelec", n_electrons)
        self.set("L", n_orbitals)

        if n_orbitals > 50 and self._storage_dir:
            print("More than 50 orbitals detected")
            print(f"Dumping boundaries to {self._storage_dir}")
            print("If you want to disable automatic dumping, call <dmrg_object>.parameter.disable_storage_dir()")
            self.set("storagedir", self._storage_dir, verbose=False)

        self._make_hf_occupation(n_orbitals, n_electrons, spin)
        self._make_site_types(n_orbitals)

    def set_integral_file(self, integral_file: str):
        """Set integral_file.

        Parameters
        ----------
        integral_file : str
            path to fcidump
        """
        self.erase("integrals", verbose=False)
        self.set("integral_file", integral_file)

    def _make_hf_occupation(self, n_orbitals: int, n_electrons: int, spin: int = 0):
        """Make Hf occupation.

        Parameters
        ----------
        n_orbitals : int
            number of orbitals
        n_electrons : int
            number of electrons
        spin : int
            total spin (2S)
        """
        doubly_occupied_electrons = n_electrons - spin
        singly_occupied_electrons = spin
        assert doubly_occupied_electrons % 2 == 0

        occupation = ""
        for _ in range(n_orbitals):
            if doubly_occupied_electrons > 0:
                occupation += "4,"
                doubly_occupied_electrons -= 2
            elif singly_occupied_electrons > 0:
                occupation += "3,"
                singly_occupied_electrons -= 1
            else:
                occupation += "1,"
        self.set("hf_occ", occupation[:-1])

    def get(self, parameter_name: str) -> Any:
        """Set any parameter in DmrgParameters.

        Parameters
        ----------
        parameter_name : str
            name of the parameter
        verbose : bool, default = True
            verbosity option

        Returns
        -------
        prameter_val : Any
            value of the parameter
        """
        if parameter_name in self._parameter_dict:
            return self._parameter_dict[parameter_name]
        return None

    def set(self, parameter_name: str, parameter_value: Any, verbose: bool = True):
        """Set any parameter in DmrgParameters.

        Parameters
        ----------
        parameter_name : str
            name of the parameter
        parameter_value : Any
            value of the parameter
        verbose : bool, default = True
            verbosity option
        """
        if parameter_name in self._parameter_dict and verbose is True:
            if self._parameter_dict[parameter_name] != parameter_value:
                message = f"{parameter_name} already set with value"
                message += f" {self._parameter_dict[parameter_name]}; new value {parameter_value}"
                print(message)

        self._parameter_dict[parameter_name] = parameter_value
        self._parameters.set(parameter_name, parameter_value)

    def erase(self, parameter_name: str, verbose: bool = True):
        """Erase any parameter from DmrgParameters.

        Parameters
        ----------
        parameter_name : str
            name of the parameter
        verbose : bool, default = True
            verbosity option
        """
        if verbose is True:
            message = f"{parameter_name} will be removed"
            print(message)

        if parameter_name in self._parameter_dict:
            del self._parameter_dict[parameter_name]

        self._parameters.erase(parameter_name)

    def get_parameters(self) -> DmrgParameters:
        """Get DmrgParameters  object.

        Return
        ------
        prameters : DmrgParameters
            the qcmaquis parameters binding
        """
        return self._parameters

    def get_parameters_dict(self) -> Dict[str, Any]:
        """Get dict with parameters.

        Return
        ------
        prameters : Dict[str, Any]
            all set parameters
        """
        return self._parameter_dict
