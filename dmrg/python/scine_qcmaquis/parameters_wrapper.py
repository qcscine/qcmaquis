"""Blub."""

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
    """Wrapper for DmrgParameters Python Interface."""

    def __init__(self):
        """Constructor."""
        self._parameters = DmrgParameters()
        """dmrg parameters map"""
        self._parameter_dict: Dict[str, Any] = {}
        """for reference, same as dmrg parameters but native python"""
        self._checkpoint_path = "checkpoint_gs"
        """name and path to store the checkpoint file."""
        self._results_path = "results_file.h5"
        """name and path to store the results file."""
        self._excited_state_name = self._checkpoint_path[:-3] + "ex0"
        """name for excited states checkpoint files."""
        # set default parameters
        self._set_defaults()

    def set_result_path(self, path: str):
        """Set path and name to checkpoint file.

        Parameters
        ----------
        path : str
            the path
        """
        self.set("resultfile", path)

    def set_checkpoint_path(self, path: str):
        """Set path and name to checkpoint file.

        Parameters
        ----------
        path : str
            the path
        """
        self._checkpoint_path = path
        self.set("chkpfile", self._checkpoint_path)

    def _set_defaults(self):
        """Set default parameters.

        Note
        ----
        "MODEL" = "quantum_chemistry"
        "init_type" = "default"
        "irrep" = 0
        "nsweeps" = 20
        "max_bond_dimension" = 2000
        "optimization" = "twosite"
        "symmetry" = "su2u1pg"
        "CONSERVED_QUANTUMNUMBERS" = "Nup,Ndown"
        "lattice_library" = "coded"
        "model_library" = "coded"
        "LATTICE" = "orbitals"
        "COMPLEX" = True

        "integrals" are set one arbitrary entry, since QcMaquis requires this keyword to be set.
        However it's not used anywhere, since we update the integral map later anyways.
        """
        self.set("MODEL", "quantum_chemistry")
        self.set("init_type", "default")
        self.set("irrep", 0)
        self.set("nsweeps", 100)
        # should be okay, due to truncation
        self.set("max_bond_dimension", 3000)
        self.set("optimization", "twosite")
        self.set("truncation_initial", 1e-6)
        self.set("truncation_final", 1e-10)
        self.set("conv_thresh", 1e-6)
        # TODO: this alpha is noise ???
        self.set("alpha_main", 1e-6)
        self.set("alpha_final", 1e-16)
        # self.set("symmetry", "2u1pg")
        self.set("symmetry", "su2u1pg")
        self.set("CONSERVED_QUANTUMNUMBERS", "Nup,Ndown")
        self.set("lattice_library", "coded")
        self.set("model_library", "coded")
        self.set("LATTICE", "orbitals")
        # self.set("COMPLEX", False)

        # dummy parameters for qcmaquis
        # dmrg requires integrals, integrals_binary or integrals_file to be set
        # already in constructor
        # Here we update the integrals later anyways with a new integral map
        # TODO: change this behavior in qcmaquis
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
        "transcorrelated_hamiltonian" = "yes"
        "time_step" = 10.
        "propagator_maxiter" = 10
        "imaginary_time" = "yes"
        "TD_backpropagation" = "no"
        "simulation_type" = "evolve"
        "COMPLEX" = True
        "time_units" = "fs"
        """
        # self.set("time_step", 1.)
        # self.set("propagator_maxiter", 10)
        # self.set("imaginary_time", "yes")
        # self.set("TD_backpropagation", "no")
        # self.set("simulation_type", "evolve")
        # self.set("propagator_accuracy", 1.0E-10)
        # self.set("symmetry", "2u1")
        # self.set("COMPLEX", True)
        # self.set("time_units", "fs")

        self.set("simulation_type", "evolve")
        self.set("transcorrelated_hamiltonian", "yes")
        # self.set("transcorrelated_hamiltonian", True)
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
            # self.set("resultfile", self._checkpoint_path + ".h5")
            self.set("chkpfile", self._checkpoint_path)
            # self.set("n_ortho_states", n_excited_states, verbose=False)
            # self.set("ortho_states", self._checkpoint_path, verbose=False)
        else:
            # self.set("resultfile", self._excited_state_name[:-1] + str(n_excited_states) + ".h5")
            self.set("chkpfile", self._excited_state_name[:-1] + str(n_excited_states))
            self.set("n_ortho_states", n_excited_states, verbose=False)
            self.set("ortho_states", self._checkpoint_path, verbose=False)
            self._checkpoint_path =\
                self._checkpoint_path + "," + self._excited_state_name[:-1] + str(n_excited_states)
            # self.set(
            #    "ortho_states",
            #    self._checkpoint_path,
            #    # verbose=False
            # )

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
        # if spin != 0:
        #     raise ValueError("Only implemented for singlet.")
        self.set("u1_total_charge1", int(n_electrons / 2) + int(spin / 2))
        self.set("u1_total_charge2", int(n_electrons / 2) - int(spin / 2))
        self.set("spin", spin)
        self.set("nelec", n_electrons)
        self.set("L", n_orbitals)

        self._make_hf_occupation(n_orbitals, n_electrons, spin)
        self._make_site_types(n_orbitals)
        # integral_file = "fcidump in correct symmetry"

    def set_integral_file(self, integral_file: str):
        """Set integral_file."""
        self.set("integral_file", integral_file)

    def _make_hf_occupation(self, n_orbitals: int, n_electrons: int, spin: int = 0):
        """Make Hf occupation."""

        occupation = ""
        for _ in range(n_orbitals):
            if n_electrons == 1:
                occupation += "3,"
                n_electrons -= 1
            elif n_electrons > 0:
                occupation += "4,"
                n_electrons -= 2
            else:
                occupation += "1,"
        self.set("hf_occ", occupation[:-1])

    def set(self, parameter_name: str, parameter_value: Any, verbose: bool = True):
        """Set any parameter in DmrgParameters."""
        if parameter_name in self._parameter_dict and verbose is True:
            message = f"{parameter_name} already set with value"
            message += f" {self._parameter_dict[parameter_name]}; new value {parameter_value}"
            print(message)

        self._parameter_dict[parameter_name] = parameter_value
        self._parameters.set(parameter_name, parameter_value)

    def erase(self, parameter_name: str, verbose: bool = True):
        """Erase any parameter from DmrgParameters."""
        if parameter_name in self._parameter_dict and verbose is True:
            message = f"{parameter_name} will be removed"
            print(message)
            del self._parameter_dict[parameter_name]

        self._parameters.erase(parameter_name)

    def get_parameters(self) -> DmrgParameters:
        """Get DmrgParameters  object."""
        return self._parameters

    def get_parameters_dict(self) -> Dict[str, Any]:
        """Get dict with parameters."""
        return self._parameter_dict
