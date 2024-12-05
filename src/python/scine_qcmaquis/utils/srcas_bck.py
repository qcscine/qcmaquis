from typing import Dict, List, Optional, Tuple

import numpy as np

from scine_qcmaquis import QCMaquis


class DeterminantQueen:

    def __init__(self, determinant_string: str):
        pass


class Srcas:

    def __init__(self, excitation_speed: float = 0.1) -> None:
        self.seed: int = 42
        self.overlap_threshold: float = 1e-5
        self.target_completeness: float = 0.95
        self.max_samples: int = 1000000

        self._current_completeness: float = 0.0
        self._n_samples: int = 0
        self._sampled_dets: Dict[str, float] = {}
        self._queen_det: str
        self._queen_coeff: float

        self._np_generator = np.random.default_rng(self.seed)

        self._excitation_speed: float = 1 - excitation_speed

        # TODO: check if this gets to big
        self._too_small_dets: Dict[str, float] = {}

    def set_seed(self, seed: int):
        """Set seed.

        Parameters
        ----------
        seed : int
            seed for random number generation
        """
        self.seed = seed
        self._np_generator = np.random.default_rng(self.seed)

    def _check_symmetry(self, det: List[int], norb: int, spin: int, nelec: int) -> bool:
        new_norb, new_spin, new_nelec = self._get_det_props(det)
        sanity_check = True
        for i in det:
            if i > 4 or i < 1:
                sanity_check = False
                break
        return norb == new_norb and new_spin == spin and sanity_check and nelec == new_nelec

    def _get_det_props(self, det: List[int]) -> Tuple[int, int, int]:
        spin = 0
        norb = 0
        nelec = 0
        for i in det:
            norb += 1
            if i == 1:
                pass
            elif i == 2:
                nelec += 1
                spin -= 1
            elif i == 3:
                nelec += 1
                spin += 1
            elif i == 4:
                nelec += 2
        return norb, spin, nelec

    # this is stupied, because this has to be generated for every excitation, even if nothing changed
    # it would be smarter to precompute and copy determinant information

    def _update_determinant(self, det: List[int]) -> List[int]:
        """Generate a new determinant

        Parameters
        ----------
        det : List[int]
            for example [4, 4, 4, 4, 1, 1, 1]
        """

        # print("update")
        norb, spin, nelec = self._get_det_props(det)

        n_excited_electrons = self._np_generator.geometric(self._excitation_speed)
        occ_orbs = []
        vir_orbs = []
        # positive indices indicate alpha occupation
        # negative indices indicate beta occupation
        for i, orbital_occupation in enumerate(det):
            i += 1
            if orbital_occupation == 4:
                occ_orbs.append(i)
                occ_orbs.append(-i)
            elif orbital_occupation == 3:
                occ_orbs.append(i)
                vir_orbs.append(-i)
            elif orbital_occupation == 2:
                occ_orbs.append(-i)
                vir_orbs.append(i)
            elif orbital_occupation == 1:
                vir_orbs.append(i)
                vir_orbs.append(-i)

        symmetry_check = False
        while not symmetry_check:
            new_det = det[:]
            for i in range(n_excited_electrons):
                # index of randomly choosen orbital (note index can be negative)
                # TODO: use s1 for weighted choice
                choosen_occ_orb_index = self._np_generator.choice(occ_orbs)
                choosen_vir_orb_index = self._np_generator.choice(vir_orbs)
                if choosen_occ_orb_index < 0:
                    new_det[abs(choosen_occ_orb_index) - 1] -= 1
                else:
                    new_det[abs(choosen_occ_orb_index) - 1] -= 2

                if choosen_vir_orb_index < 0:
                    new_det[abs(choosen_vir_orb_index) - 1] += 2
                else:
                    new_det[abs(choosen_vir_orb_index) - 1] += 1

            symmetry_check = self._check_symmetry(new_det, norb, spin, nelec)
        return new_det

    def _generate_new_determinant(self) -> str:
        """Generate a new determinant from current Queen

        Returns
        -------
        new_det : str
            new determinant
        """
        # print("generate")
        new_det_list = [int(i) for i in self._queen_det.replace(",", "")]
        new_det_list = self._update_determinant(new_det_list)
        new_det = ",".join([str(i) for i in new_det_list])
        # if new_det in self._sampled_dets or new_det in self._too_small_dets:
        #     # TODO: check if new det is new queen
        #     return self._generate_new_determinant()

        return new_det

    def run(self, dmrg: QCMaquis, initial_det: str, additional_dets: Optional[List[str]] = None):
        """Run SRCAS

        Parameters
        ----------
        dmrg : QCMaquis
            dmrg object
        initial_det : str
            initial determinant (for example HF, e.g. 44444111)
        additional_dets : Optional[List[str]]
            list with determinants to be measured
        """
        overlap = dmrg.get_ci_coefficient(initial_det)
        self._current_completeness += overlap*overlap
        self._sampled_dets[initial_det] = overlap
        self._n_samples += 1
        self._queen_det = initial_det
        print(self._n_samples, overlap, initial_det, self._current_completeness)

        if additional_dets:
            for det in additional_dets:
                overlap = dmrg.get_ci_coefficient(det)
                self._current_completeness += overlap*overlap
                self._sampled_dets[det] = overlap
                self._n_samples += 1

        self._loop(dmrg, initial_det)

    def _loop(self, dmrg: QCMaquis, det: str):

        while self._n_samples < self.max_samples:
            det = self._generate_new_determinant()
            if det not in self._sampled_dets and det not in self._too_small_dets:
                # TODO: check if new det is new queen
                # return self._generate_new_determinant()
                overlap = dmrg.get_ci_coefficient(det)

                if overlap >= self.overlap_threshold:
                    self._sampled_dets[det] = overlap
                    self._current_completeness += overlap*overlap
                    self._n_samples += 1
                    print(self._n_samples, overlap, det, self._current_completeness)
                else:
                    self._too_small_dets[det] = overlap
                    self._current_completeness += overlap*overlap
                    # print(self._n_samples, overlap, det, self._current_completeness)
                    # self._n_samples += 1

            if self._current_completeness >= self.target_completeness:
                break


"""





















"""
