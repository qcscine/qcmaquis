"""This module provides the public QCMaquis DMRG interface and related functions."""
from pathlib import Path
from typing import Any, List, Tuple, Union

import numpy as np

# pylint: disable=import-error
from .dmrg_wrapper import DmrgWrapper
from .integral_wrapper import ComplexTCIntegralMap, IntegralMap, IntegralMapWrapper, IntegralType, TCIntegralMap
from .parameters_wrapper import ExcitedStates, ParametersWrapper
from .utils.ci_coeffs import (make_doubles_aa, make_doubles_ab, make_doubles_bb, make_ref, make_singles_aa,
                              make_singles_bb)

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
    # TODO: add __slots__

    def __init__(self):
        """Construct Wrapper."""
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
        """Set bond dimension."""
        self._parameters.set("max_bond_dimension", bond_dimension)

    def set_excited_states(
        self,
        method: ExcitedStates = ExcitedStates.ORTHO,
        n_excited_states: int = None,
        feast_window: List[float] = None
    ):
        """Enable Excited States."""
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

    def get_entropies(self):
        self._dmrg.measure()
        return self._dmrg.entropies()

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
        # self._energy = self._dmrg.get_dmrg().energy()
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

    def run(self, n_orbitals: int, n_electrons: int, spin: int = 0, n_states: int = None, fiedler: bool = False):
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
        """
        for i in self._parameters._parameter_dict:
            print(i, self._parameters._parameter_dict[i])

        # excited states
        # TODO: If should not be required here ...
        if n_states is not None:
            self._energy = []
            self._parameters.set_system(n_orbitals, n_electrons, spin)
            self._parameters.set_excited_states_ortho(0)

            # fiedler needs to be implemented in the interface
            if fiedler is True:
                raise NotImplementedError
                # fiedler_orderer = DmrgWrapper()
                # self._parameters.set("MEASURE[ChemEntropy]", True)
                # fiedler_orderer.set_parameters(self._parameters)
                # if "integral_file" not in self._parameters._parameter_dict:
                #     fiedler_orderer.set_integrals(self._integral_map)
                # # orbital_order = fiedler_orderer.get_fiedler(n_states=n_states)
                # orbital_order = fiedler_orderer.get_fiedler()
                # self._parameters.erase("MEASURE[ChemEntropy]")
                # self._parameters.set("orbital_order", orbital_order)

                # orbital_order = self._dmrg.get_fiedler(n_states=n_states)
                # self._parameters.set("orbital_order", orbital_order)

            self._dmrg.set_parameters(self._parameters)

            if "integral_file" not in self._parameters._parameter_dict:
                self._dmrg.set_integrals(self._integral_map)

            self._dmrg.run()
            self._energy.append(self._dmrg.get_energy())

            for i in range(1, n_states):
                self._parameters.set_system(n_orbitals, n_electrons, spin)
                self._parameters.set_excited_states_ortho(i)
                self._dmrg.set_parameters(self._parameters)
                if "integral_file" not in self._parameters._parameter_dict:
                    self._dmrg.set_integrals(self._integral_map)
                self._dmrg.run()
                self._energy.append(self._dmrg.get_energy())

        # ground state only
        else:

            self._parameters.set_system(n_orbitals, n_electrons, spin)

            # fiedler needs to be implemented in the interface
            if fiedler is True:
                raise NotImplementedError
                # fiedler_orderer = DmrgWrapper()
                # self._parameters.set("MEASURE[ChemEntropy]", True)
                # fiedler_orderer.set_parameters(self._parameters)
                # if "integral_file" not in self._parameters._parameter_dict:
                #     fiedler_orderer.set_integrals(self._integral_map)
                # orbital_order = fiedler_orderer.get_fiedler()
                # self._parameters.erase("MEASURE[ChemEntropy]")
                # self._parameters.set("orbital_order", orbital_order)

            # self._parameters.set("orbital_order", "1,2")
            self._dmrg.set_parameters(self._parameters)
            if "integral_file" not in self._parameters._parameter_dict:
                self._dmrg.set_integrals(self._integral_map)
            self._dmrg.run()
            self._energy = self._dmrg.get_energy()

    def update_integrals(self, integral_map: IntegralMap):
        """Update integrals.

        Parameters
        ----------
        integral_map : IntegralMap
            the integral map in corresponding symmetry, e.g.
            8-fold for conventional
            2-fold for transcorrelated
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

    # unused
    def dummy_run_excited_states(self, n_excited_states: int):
        """Dummy function."""
        energies = []
        self._parameters.set_system(2, 2, 0)
        self._dmrg.set_parameters(self._parameters.get())
        self._dmrg.get_dmrg().update_integrals(self._integral_map.get())
        for i in range(n_excited_states):
            self._parameters.set_excited_states_ortho(i)
            self._dmrg.set_parameters(self._parameters.get())
            self._dmrg.get_dmrg().optimize()
            energies.append(self.get_energy())
        print(energies)

    def init_dmrg(self, checkpoint: str, norb: int, nelec: int, spin: int):
        self._parameters.set_system(norb, nelec, spin)
        self._parameters.set_checkpoint_path(checkpoint)
        self._dmrg.set_parameters(self._parameters)

    def get_ci_coefficient(self, determinant_string: str):
        return self._dmrg.get_ci_coefficient(determinant_string)

    def get_singles_and_doubles(self, nocc: int, norb: int):
        nvir = norb - nocc
        hf_string, sign = make_ref(nocc, norb)
        coeff_hf = self._dmrg.get_ci_coefficient(hf_string).real
        print(hf_string, coeff_hf)
        singles = np.ndarray((nocc * 2, nvir * 2))
        singles.fill(0)
        for i in range(nocc):
            for a in range(nvir):
                singles_aa_string, sign_aa = make_singles_aa(nocc, norb, i, a)
                coeff_aa = self._dmrg.get_ci_coefficient(singles_aa_string)
                singles_bb_string, sign_bb = make_singles_bb(nocc, norb, i, a)
                coeff_bb = self._dmrg.get_ci_coefficient(singles_bb_string)
                singles[i * 2, a * 2] = coeff_aa.real
                singles[i * 2 + 1, a * 2 + 1] = coeff_bb.real

        doubles = np.ndarray((nocc * 2, nocc * 2, nvir * 2, nvir * 2))
        doubles.fill(0)
        for i in range(nocc):
            for j in range(nocc):
                for a in range(nvir):
                    for b in range(nvir):
                        if nocc > 1 and nvir > 1 and i != j and a != b:
                            doubles_aa_string, sign_aa = make_doubles_aa(nocc, norb, i, j, a, b)
                            coeff_aa = self._dmrg.get_ci_coefficient(doubles_aa_string)
                            doubles_bb_string, sign_bb = make_doubles_bb(nocc, norb, i, j, a, b)
                            coeff_bb = self._dmrg.get_ci_coefficient(doubles_bb_string)
                            doubles[i * 2, j * 2, a * 2, b * 2] = coeff_aa.real
                            doubles[i * 2 + 1, j * 2 + 1, a * 2 + 1, b * 2 + 1] = coeff_bb.real
                        doubles_ab_string, sign_ab = make_doubles_ab(nocc, norb, i, j, a, b)
                        coeff_ab = self._dmrg.get_ci_coefficient(doubles_ab_string)
                        doubles[i * 2, j * 2 + 1, a * 2, b * 2 + 1] = coeff_ab.real
                        doubles[i * 2 + 1, j * 2, a * 2 + 1, b * 2] = coeff_ab.real

                        # if nocc > 1 and nvir > 1:
                        #     print(doubles_aa_string, coeff_aa)
                        #     print(doubles_bb_string, coeff_bb)
                        # print(doubles_ab_string, coeff_ab)
        return coeff_hf, singles, doubles


if __name__ == "__main__":

    test_current_path = str(Path.cwd()) + "/python_tests"

    integrals = ComplexTCIntegralMap()
    integrals.set((1, 1, 0, 0, 0, 0), -1.9410228773342559e+00)
    integrals.set((1, 2, 0, 0, 0, 0), -3.1641663736652514e-01)
    integrals.set((2, 1, 0, 0, 0, 0), -3.1641663736652437e-01)
    integrals.set((2, 2, 0, 0, 0, 0), -9.0227670561564222e-02)
    integrals.set((3, 3, 0, 0, 0, 0), 7.8499729043522848e-01)
    integrals.set((4, 4, 0, 0, 0, 0), 7.8499729043522848e-01)
    integrals.set((5, 5, 0, 0, 0, 0), 7.8499729043522892e-01)
    integrals.set((1, 1, 1, 1, 0, 0), 1.0183556583001547e+00)
    integrals.set((1, 1, 2, 1, 0, 0), 3.2114288110690942e-01)
    integrals.set((1, 1, 2, 2, 0, 0), 8.5366263969801692e-01)
    integrals.set((1, 1, 3, 3, 0, 0), 9.5116233984037879e-01)
    integrals.set((1, 1, 4, 4, 0, 0), 9.5116233984037879e-01)
    integrals.set((1, 1, 5, 5, 0, 0), 9.5116233984037879e-01)
    integrals.set((1, 2, 1, 1, 0, 0), 3.0476675403355274e-01)
    integrals.set((1, 2, 1, 2, 0, 0), 2.5174364834784085e-01)
    integrals.set((1, 2, 2, 1, 0, 0), 2.2385496280331504e-01)
    integrals.set((1, 2, 2, 2, 0, 0), 2.4714874264182368e-01)
    integrals.set((1, 2, 3, 3, 0, 0), 2.1082535928181162e-01)
    integrals.set((1, 2, 4, 4, 0, 0), 2.1082535928181162e-01)
    integrals.set((1, 2, 5, 5, 0, 0), 2.1082535928181179e-01)
    integrals.set((1, 3, 1, 3, 0, 0), 1.8330678004452328e-01)
    integrals.set((1, 3, 2, 3, 0, 0), 2.5134760449013423e-02)
    integrals.set((1, 3, 3, 1, 0, 0), 1.8098179143498977e-01)
    integrals.set((1, 3, 3, 2, 0, 0), 4.4024363522944421e-02)
    integrals.set((1, 4, 1, 4, 0, 0), 1.8330678004452328e-01)
    integrals.set((1, 4, 2, 4, 0, 0), 2.5134760449013416e-02)
    integrals.set((1, 4, 4, 1, 0, 0), 1.8098179143498977e-01)
    integrals.set((1, 4, 4, 2, 0, 0), 4.4024363522944407e-02)
    integrals.set((1, 5, 1, 5, 0, 0), 1.8330678004452333e-01)
    integrals.set((1, 5, 2, 5, 0, 0), 2.5134760449013447e-02)
    integrals.set((1, 5, 5, 1, 0, 0), 1.8098179143498982e-01)
    integrals.set((1, 5, 5, 2, 0, 0), 4.4024363522944449e-02)
    integrals.set((2, 1, 2, 1, 0, 0), 1.9596627725878948e-01)
    integrals.set((2, 1, 3, 3, 0, 0), 2.4775473575494533e-01)
    integrals.set((2, 1, 4, 4, 0, 0), 2.4775473575494533e-01)
    integrals.set((2, 1, 5, 5, 0, 0), 2.4775473575494544e-01)
    integrals.set((2, 2, 2, 1, 0, 0), 2.5546841384748981e-01)
    integrals.set((2, 2, 2, 2, 0, 0), 7.4940475756065872e-01)
    integrals.set((2, 2, 3, 3, 0, 0), 7.8286400008585233e-01)
    integrals.set((2, 2, 4, 4, 0, 0), 7.8286400008585233e-01)
    integrals.set((2, 2, 5, 5, 0, 0), 7.8286400008585233e-01)
    integrals.set((2, 3, 2, 3, 0, 0), 3.6676840894247775e-02)
    integrals.set((2, 3, 3, 1, 0, 0), 4.7410516610267971e-03)
    integrals.set((2, 3, 3, 2, 0, 0), 2.5466242990294468e-02)
    integrals.set((2, 4, 2, 4, 0, 0), 3.6676840894247768e-02)
    integrals.set((2, 4, 4, 1, 0, 0), 4.7410516610267971e-03)
    integrals.set((2, 4, 4, 2, 0, 0), 2.5466242990294475e-02)
    integrals.set((2, 5, 2, 5, 0, 0), 3.6676840894247782e-02)
    integrals.set((2, 5, 5, 1, 0, 0), 4.7410516610267902e-03)
    integrals.set((2, 5, 5, 2, 0, 0), 2.5466242990294482e-02)
    integrals.set((3, 1, 3, 1, 0, 0), 1.7865680282545626e-01)
    integrals.set((3, 2, 3, 1, 0, 0), 2.3630654734957762e-02)
    integrals.set((3, 2, 3, 2, 0, 0), 1.4255645086341169e-02)
    integrals.set((3, 3, 3, 3, 0, 0), 9.3158111092842422e-01)
    integrals.set((3, 3, 4, 4, 0, 0), 9.3158111092842422e-01)
    integrals.set((3, 3, 5, 5, 0, 0), 9.3158111092842444e-01)
    integrals.set((4, 1, 4, 1, 0, 0), 1.7865680282545626e-01)
    integrals.set((4, 2, 4, 1, 0, 0), 2.3630654734957766e-02)
    integrals.set((4, 2, 4, 2, 0, 0), 1.4255645086341169e-02)
    integrals.set((4, 4, 4, 4, 0, 0), 9.3158111092842422e-01)
    integrals.set((4, 4, 5, 5, 0, 0), 9.3158111092842444e-01)
    integrals.set((5, 1, 5, 1, 0, 0), 1.7865680282545632e-01)
    integrals.set((5, 2, 5, 1, 0, 0), 2.3630654734957790e-02)
    integrals.set((5, 2, 5, 2, 0, 0), 1.4255645086341176e-02)
    integrals.set((5, 5, 5, 5, 0, 0), 9.3158111092842388e-01)
    # integrals.set((1, 1, 1, 1, 1, 1), 1.9115381141742809e-03)
    # integrals.set((1, 1, 1, 1, 1, 2), -2.5509485850552923e-04)
    # integrals.set((1, 1, 1, 1, 2, 2), 2.7878203188368704e-03)
    # integrals.set((1, 1, 1, 1, 3, 3), 1.3963811163335059e-03)
    # integrals.set((1, 1, 1, 1, 4, 4), 1.3963811163335059e-03)
    # integrals.set((1, 1, 1, 1, 5, 5), 1.3963811163335066e-03)
    # integrals.set((1, 1, 1, 2, 1, 2), 5.5178136572181643e-04)
    # integrals.set((1, 1, 1, 2, 2, 2), 2.7868471658280803e-04)
    # integrals.set((1, 1, 1, 2, 3, 3), -9.4939835529712555e-04)
    # integrals.set((1, 1, 1, 2, 4, 4), -9.4939835529712566e-04)
    # integrals.set((1, 1, 1, 2, 5, 5), -9.4939835529712501e-04)
    # integrals.set((1, 1, 1, 3, 1, 3), -7.3106461867098875e-05)
    # integrals.set((1, 1, 1, 3, 2, 3), 1.9549001030377992e-04)
    # integrals.set((1, 1, 1, 4, 1, 4), -7.3106461867098740e-05)
    # integrals.set((1, 1, 1, 4, 2, 4), 1.9549001030377946e-04)
    # integrals.set((1, 1, 1, 5, 1, 5), -7.3106461867098631e-05)
    # integrals.set((1, 1, 1, 5, 2, 5), 1.9549001030377946e-04)
    # integrals.set((1, 1, 2, 2, 2, 2), 3.3096558215331652e-03)
    # integrals.set((1, 1, 2, 2, 3, 3), 2.2226209281468959e-03)
    # integrals.set((1, 1, 2, 2, 4, 4), 2.2226209281468959e-03)
    # integrals.set((1, 1, 2, 2, 5, 5), 2.2226209281468973e-03)
    # integrals.set((1, 1, 2, 3, 2, 3), 1.2063223168393608e-04)
    # integrals.set((1, 1, 2, 4, 2, 4), 1.2063223168393637e-04)
    # integrals.set((1, 1, 2, 5, 2, 5), 1.2063223168393584e-04)
    # integrals.set((1, 1, 3, 3, 3, 3), 1.3115851469828498e-03)
    # integrals.set((1, 1, 3, 3, 4, 4), 1.3115851469828498e-03)
    # integrals.set((1, 1, 3, 3, 5, 5), 1.3115851469828502e-03)
    # integrals.set((1, 1, 4, 4, 4, 4), 1.3115851469828502e-03)
    # integrals.set((1, 1, 4, 4, 5, 5), 1.3115851469828502e-03)
    # integrals.set((1, 1, 5, 5, 5, 5), 1.3115851469828492e-03)
    # integrals.set((1, 2, 1, 2, 1, 2), 1.0902359255668213e-03)
    # integrals.set((1, 2, 1, 2, 2, 2), 7.2073724685171848e-04)
    # integrals.set((1, 2, 1, 2, 3, 3), -1.8874709287954700e-04)
    # integrals.set((1, 2, 1, 2, 4, 4), -1.8874709287954776e-04)
    # integrals.set((1, 2, 1, 2, 5, 5), -1.8874709287954603e-04)
    # integrals.set((1, 2, 1, 3, 1, 3), 2.9814078664790850e-04)
    # integrals.set((1, 2, 1, 3, 2, 3), 1.8134560599690198e-04)
    # integrals.set((1, 2, 1, 4, 1, 4), 2.9814078664790818e-04)
    # integrals.set((1, 2, 1, 4, 2, 4), 1.8134560599690182e-04)
    # integrals.set((1, 2, 1, 5, 1, 5), 2.9814078664790867e-04)
    # integrals.set((1, 2, 1, 5, 2, 5), 1.8134560599690163e-04)
    # integrals.set((1, 2, 2, 2, 2, 2), 6.4345903846134921e-04)
    # integrals.set((1, 2, 2, 2, 3, 3), -3.9185982211398463e-04)
    # integrals.set((1, 2, 2, 2, 4, 4), -3.9185982211398555e-04)
    # integrals.set((1, 2, 2, 2, 5, 5), -3.9185982211398360e-04)
    # integrals.set((1, 2, 2, 3, 2, 3), 1.0680638237659020e-04)
    # integrals.set((1, 2, 2, 4, 2, 4), 1.0680638237659016e-04)
    # integrals.set((1, 2, 2, 5, 2, 5), 1.0680638237659009e-04)
    # integrals.set((1, 2, 3, 3, 3, 3), -1.2390840425682437e-03)
    # integrals.set((1, 2, 3, 3, 4, 4), -1.2390840425682435e-03)
    # integrals.set((1, 2, 3, 3, 5, 5), -1.2390840425682446e-03)
    # integrals.set((1, 2, 4, 4, 4, 4), -1.2390840425682433e-03)
    # integrals.set((1, 2, 4, 4, 5, 5), -1.2390840425682442e-03)
    # integrals.set((1, 2, 5, 5, 5, 5), -1.2390840425682450e-03)
    # integrals.set((1, 3, 1, 3, 2, 2), -1.1571371462836206e-04)
    # integrals.set((1, 3, 1, 3, 3, 3), 2.3764355447230526e-05)
    # integrals.set((1, 3, 1, 3, 4, 4), -3.2866895450366616e-04)
    # integrals.set((1, 3, 1, 3, 5, 5), -3.2866895450366605e-04)
    # integrals.set((1, 3, 1, 4, 3, 4), 1.7621665497544853e-04)
    # integrals.set((1, 3, 1, 5, 3, 5), 1.7621665497544859e-04)
    # integrals.set((1, 3, 2, 2, 2, 3), 1.7274800027940331e-04)
    # integrals.set((1, 3, 2, 3, 3, 3), 2.0044058066931485e-04)
    # integrals.set((1, 3, 2, 3, 4, 4), 5.7447590338582306e-05)
    # integrals.set((1, 3, 2, 3, 5, 5), 5.7447590338582225e-05)
    # integrals.set((1, 3, 2, 4, 3, 4), 7.1496495165366387e-05)
    # integrals.set((1, 3, 2, 5, 3, 5), 7.1496495165366319e-05)
    # integrals.set((1, 4, 1, 4, 2, 2), -1.1571371462836227e-04)
    # integrals.set((1, 4, 1, 4, 3, 3), -3.2866895450366616e-04)
    # integrals.set((1, 4, 1, 4, 4, 4), 2.3764355447231122e-05)
    # integrals.set((1, 4, 1, 4, 5, 5), -3.2866895450366600e-04)
    # integrals.set((1, 4, 1, 5, 4, 5), 1.7621665497544859e-04)
    # integrals.set((1, 4, 2, 2, 2, 4), 1.7274800027940299e-04)
    # integrals.set((1, 4, 2, 3, 3, 4), 7.1496495165366387e-05)
    # integrals.set((1, 4, 2, 4, 3, 3), 5.7447590338581601e-05)
    # integrals.set((1, 4, 2, 4, 4, 4), 2.0044058066931458e-04)
    # integrals.set((1, 4, 2, 4, 5, 5), 5.7447590338581954e-05)
    # integrals.set((1, 4, 2, 5, 4, 5), 7.1496495165366332e-05)
    # integrals.set((1, 5, 1, 5, 2, 2), -1.1571371462836157e-04)
    # integrals.set((1, 5, 1, 5, 3, 3), -3.2866895450366573e-04)
    # integrals.set((1, 5, 1, 5, 4, 4), -3.2866895450366567e-04)
    # integrals.set((1, 5, 1, 5, 5, 5), 2.3764355447231285e-05)
    # integrals.set((1, 5, 2, 2, 2, 5), 1.7274800027940299e-04)
    # integrals.set((1, 5, 2, 3, 3, 5), 7.1496495165366305e-05)
    # integrals.set((1, 5, 2, 4, 4, 5), 7.1496495165366373e-05)
    # integrals.set((1, 5, 2, 5, 3, 3), 5.7447590338581683e-05)
    # integrals.set((1, 5, 2, 5, 4, 4), 5.7447590338581872e-05)
    # integrals.set((1, 5, 2, 5, 5, 5), 2.0044058066931436e-04)
    # integrals.set((2, 2, 2, 2, 2, 2), 3.5462315265775328e-03)
    # integrals.set((2, 2, 2, 2, 3, 3), 2.7438361990735384e-03)
    # integrals.set((2, 2, 2, 2, 4, 4), 2.7438361990735375e-03)
    # integrals.set((2, 2, 2, 2, 5, 5), 2.7438361990735388e-03)
    # integrals.set((2, 2, 2, 3, 2, 3), 1.1536186423993642e-04)
    # integrals.set((2, 2, 2, 4, 2, 4), 1.1536186423993647e-04)
    # integrals.set((2, 2, 2, 5, 2, 5), 1.1536186423993645e-04)
    # integrals.set((2, 2, 3, 3, 3, 3), 2.0381744505766700e-03)
    # integrals.set((2, 2, 3, 3, 4, 4), 2.0381744505766700e-03)
    # integrals.set((2, 2, 3, 3, 5, 5), 2.0381744505766682e-03)
    # integrals.set((2, 2, 4, 4, 4, 4), 2.0381744505766700e-03)
    # integrals.set((2, 2, 4, 4, 5, 5), 2.0381744505766682e-03)
    # integrals.set((2, 2, 5, 5, 5, 5), 2.0381744505766687e-03)
    # integrals.set((2, 3, 2, 3, 3, 3), 1.1958381178876656e-04)
    # integrals.set((2, 3, 2, 3, 4, 4), 3.4062460360147071e-05)
    # integrals.set((2, 3, 2, 3, 5, 5), 3.4062460360147139e-05)
    # integrals.set((2, 3, 2, 4, 3, 4), 4.2760675714310064e-05)
    # integrals.set((2, 3, 2, 5, 3, 5), 4.2760675714309895e-05)
    # integrals.set((2, 4, 2, 4, 3, 3), 3.4062460360147179e-05)
    # integrals.set((2, 4, 2, 4, 4, 4), 1.1958381178876691e-04)
    # integrals.set((2, 4, 2, 4, 5, 5), 3.4062460360147220e-05)
    # integrals.set((2, 4, 2, 5, 4, 5), 4.2760675714309969e-05)
    # integrals.set((2, 5, 2, 5, 3, 3), 3.4062460360146949e-05)
    # integrals.set((2, 5, 2, 5, 4, 4), 3.4062460360146922e-05)
    # integrals.set((2, 5, 2, 5, 5, 5), 1.1958381178876591e-04)
    # integrals.set((3, 3, 3, 3, 3, 3), 1.4624785379663617e-03)
    # integrals.set((3, 3, 3, 3, 4, 4), 1.4624785379663613e-03)
    # integrals.set((3, 3, 3, 3, 5, 5), 1.4624785379663634e-03)
    # integrals.set((3, 3, 4, 4, 4, 4), 1.4624785379663608e-03)
    # integrals.set((3, 3, 4, 4, 5, 5), 1.4624785379663634e-03)
    # integrals.set((3, 3, 5, 5, 5, 5), 1.4624785379663639e-03)
    # integrals.set((4, 4, 4, 4, 4, 4), 1.4624785379663600e-03)
    # integrals.set((4, 4, 4, 4, 5, 5), 1.4624785379663626e-03)
    # integrals.set((4, 4, 5, 5, 5, 5), 1.4624785379663639e-03)
    # integrals.set((5, 5, 5, 5, 5, 5), 1.4624785379663643e-03)

    # integrals = IntegralMap()
    # integrals.set((1, 1, 1, 1), 0.354237848011)
    # integrals.set((1, 1, 2, 1), -0.821703816101E-13)
    # integrals.set((2, 1, 2, 1), 0.185125251547)
    # integrals.set((2, 2, 2, 1), 0.782984788117E-13)
    # integrals.set((1, 1, 2, 2), 0.361001163519)
    # integrals.set((2, 2, 2, 2), 0.371320200119)
    # integrals.set((1, 1, 0, 0), -0.678487901790)
    # integrals.set((2, 1, 0, 0), -0.539801158857E-14)
    # integrals.set((2, 2, 0, 0), -0.653221638776)
    # integrals.set((0, 0, 0, 0), 0.176392403557)

    blub = []
    this_dmrg = MaquisDmrg()
    print("---------")
    # this_dmrg.set_transcorrelation()
    print("---------")
    # this_dmrg.update_integrals(integrals)
    print("---------")
    # this_dmrg._parameters.erase("integrals")
    # this_dmrg.
    # _parameters.set("integral_file", "/home/max/Programs/coupled_wick_scf/maquis-dmrg_python/dmrg/IntegralFile_H2_Transcorrelated")
    this_dmrg._parameters.erase("integrals")
    this_dmrg._parameters.set("nsweeps", 1)
    this_dmrg._parameters.set("max_bond_dimension", 1000)
    this_dmrg._parameters.set("integral_file", "/home/max/Programs/coupled_wick_scf/scripts/test/cc-pvdz/trans/0.0/He_cc-pvdz.FCIDUMP")
    this_dmrg._parameters.set("optimization", "singlesite")
    this_dmrg._parameters.set("simulation_type", "time_dep")
    this_dmrg._parameters.set("propagator_accuracy", 1.0E-10)
    this_dmrg._parameters.set("propagator_maxiter", 10)
    this_dmrg._parameters.set("time_step", "0.2")
    this_dmrg._parameters.set("hamiltonian_units", "Hartree")
    this_dmrg._parameters.set("time_units", "fs")
    this_dmrg._parameters.set("imaginary_time", "yes")
    this_dmrg._parameters.set("TD_backpropagation", "no")
    this_dmrg._parameters.set("transcorrelated_hamiltonian", "yes")

    # this_dmrg._parameters.set("chh", 1000)
    # this_dmrg._parameters.set("chkpfile", "/home/max/Programs/coupled_wick_scf/maquis-dmrg_python/dmrg/python/checkpoint")
    # this_dmrg.run(28, 14, fiedler=False)
    # this_dmrg.run(2, 2, fiedler=True)
    # blub.append(this_dmrg.get_energy())
    # this_dmrg.run(2, 2, fiedler=False)
    # this_dmrg.run(28, 14, fiedler=False, n_states=2)
    # blub.append(this_dmrg.get_energy())
    # print("------------ FIEDLER --------------------")
    # this_dmrg.run(2, 2, fiedler=True, n_states=2)
    # this_dmrg._parameters._checkpoint_path = "blub_gs"
    #
    # this_dmrg.run(28, 14, fiedler=True, n_states=4)
    # [-108.8661510547466, -108.7528692473976, -108.74233988457854, -108.71384766020647]
    # this_dmrg.set_feast((-108.75, -108.74,), 4)
    # this_dmrg.run(28, 14, fiedler=True)
    # this_dmrg.set_feast((-0.7, -0.5,), 8)
    this_dmrg.run(5, 2, fiedler=False)
    # print(this_dmrg._dmrg.get_ci_coefficients("3,2,1,1,1"))
    # print("--------------------")
    this_dmrg.get_singles_and_doubles(1, 2)
    print(this_dmrg._dmrg._dmrg.getCICoefficients(2))
    # print("0000000000")
    blub.append(this_dmrg.get_energy())
    # -3.882045755
    """
    this_dmrg_2 = MaquisDmrg()
    this_dmrg_2.set_feast((-0.7, -0.5,), 8)
    this_dmrg_2.update_integrals(integrals)
    this_dmrg_2.run(2, 2)
    this_dmrg_3 = MaquisDmrg()
    this_dmrg_3.update_integrals(integrals)
    this_dmrg_3.run(2, 2, n_states=2, fiedler=False)
    print(this_dmrg.get_energy())
    print(this_dmrg_2.get_energy())
    print(this_dmrg_3.get_energy())
    """
    for i in blub:
        print(i)
        print(i)
