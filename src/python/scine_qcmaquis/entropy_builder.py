from dataclasses import dataclass
from typing import List, Tuple, Union

import numpy as np


@dataclass
class Datasets:
    """
    Class to store the measurements from the QCMaquis for entropy calculations

    Note
    ----
    All class member names correspond are similar to the names in the QCMaquis

    Attributes
    ----------
    self.nup : np.ndarray
        nup measurement for orbital RDM
    self.ndown : np.ndarray
        ndown measurement for orbital RDM
    self.nupdown : np.ndarray
        nupdown measurement for orbital RDM
    self.dmUpSym : np.ndarray
        dmUpSym measurement for orbital RDM
    self.dmDownSym : np.ndarray
        dmDownSym measurement for orbital RDM
    self.nupNupSym : np.ndarray
        nupNupSym measurement for orbital RDM
    self.nupNdownSym : np.ndarray
        nupNdownSym measurement for orbital RDM
    self.ndownNupSym : np.ndarray
        ndownNupSym measurement for orbital RDM
    self.ndownNdownSym : np.ndarray
        ndownNdownSym measurement for orbital RDM
    self.doccDoccSym : np.ndarray
        doccDoccSym measurement for orbital RDM
    self.transferUpWhileDownSym : np.ndarray
        transferUpWhileDownSym measurement for orbital RDM
    self.transferDownWhileUpSym : np.ndarray
        transferDownWhileUpSym measurement for orbital RDM
    self.transferPairSym : np.ndarray
        transferPairSym measurement for orbital RDM
    self.spinflipSym : np.ndarray
        spinflipSym measurement for orbital RDM
    self.transUpDown2Asym : np.ndarray
        transUpDown2Asym measurement for orbital RDM
    self.transUpDown1Asym : np.ndarray
        transUpDown1Asym measurement for orbital RDM
    self.transDownUp2Asym : np.ndarray
        transDownUp2Asym measurement for orbital RDM
    self.transDownUp1Asym : np.ndarray
        transDownUp1Asym measurement for orbital RDM
    self.nupDoccAsym : np.ndarray
        nupDoccAsym measurement for orbital RDM
    self.doccNupAsym : np.ndarray
        doccNupAsym measurement for orbital RDM
    self.ndownDoccAsym : np.ndarray
        ndownDoccAsym measurement for orbital RDM
    self.doccNdownAsym : np.ndarray
        doccNdownAsym measurement for orbital RDM
    """

    __slots__ = [
        "nup", "ndown", "nupdown", "dmUpSym", "dmDownSym", "nupNupSym", "nupNdownSym", "ndownNupSym",
        "ndownNdownSym", "doccDoccSym", "transferUpWhileDownSym", "transferDownWhileUpSym", "transferPairSym",
        "spinflipSym", "transUpDown2Asym", "transUpDown1Asym", "transDownUp2Asym", "transDownUp1Asym",
        "nupDoccAsym", "doccNupAsym", "ndownDoccAsym", "doccNdownAsym", "occ_num_vector", "spinorCdagCSymMat",
        "spinorDoccDoccSymMat",
    ]

    def __init__(self):
        """Init."""
        self.nup: np.ndarray
        """nup measurement for orbital RDM"""
        self.ndown: np.ndarray
        """ndown measurement for orbital RDM"""
        self.nupdown: np.ndarray
        """nupdown measurement for orbital RDM"""
        self.dmUpSym: np.ndarray
        """dmUpSym measurement for orbital RDM"""
        self.dmDownSym: np.ndarray
        """dmDownSym measurement for orbital RDM"""
        self.nupNupSym: np.ndarray
        """nupNupSym measurement for orbital RDM"""
        self.nupNdownSym: np.ndarray
        """nupNdownSym measurement for orbital RDM"""
        self.ndownNupSym: np.ndarray
        """ndownNupSym measurement for orbital RDM"""
        self.ndownNdownSym: np.ndarray
        """ndownNdownSym measurement for orbital RDM"""
        self.doccDoccSym: np.ndarray
        """doccDoccSym measurement for orbital RDM"""
        self.transferUpWhileDownSym: np.ndarray
        """transferUpWhileDownSym measurement for orbital RDM"""
        self.transferDownWhileUpSym: np.ndarray
        """transferDownWhileUpSym measurement for orbital RDM"""
        self.transferPairSym: np.ndarray
        """transferPairSym measurement for orbital RDM"""
        self.spinflipSym: np.ndarray
        """spinflipSym measurement for orbital RDM"""
        self.transUpDown2Asym: np.ndarray
        """transUpDown2Asym measurement for orbital RDM"""
        self.transUpDown1Asym: np.ndarray
        """transUpDown1Asym measurement for orbital RDM"""
        self.transDownUp2Asym: np.ndarray
        """transDownUp2Asym measurement for orbital RDM"""
        self.transDownUp1Asym: np.ndarray
        """transDownUp1Asym measurement for orbital RDM"""
        self.nupDoccAsym: np.ndarray
        """nupDoccAsym measurement for orbital RDM"""
        self.doccNupAsym: np.ndarray
        """doccNupAsym measurement for orbital RDM"""
        self.ndownDoccAsym: np.ndarray
        """ndownDoccAsym measurement for orbital RDM"""
        self.doccNdownAsym: np.ndarray
        """doccNdownAsym measurement for orbital RDM"""


class EntropyBuilder:
    """Build entropies from base measurements.

    Attributes
    ----------
    self.L : int
        number of orbitals
    self.orbital_order : List[int]
        Orbital order used in calculation
    """

    __slots__ = ("L", "orbital_order", "s1_entropy", "s2_entropy", "mutual_information")

    def __init__(self, lattice_size: int, orbital_order: Union[List[int], str] = None):
        """Set orbital order to "canonical" if none is set.

        Parameters
        ----------
        lattice_size : int
            n_orbitals, L used in DMRG calculation.
        orbital_order : List[int], optional
            Orbital order, default will be canonical
        """
        self.L: int = lattice_size
        """size of the lattice"""
        if not orbital_order:
            orbital_order = [i for i in range(lattice_size)]
        elif type(orbital_order) is str:
            orbital_order = self._get_list_from_fiedler_string(orbital_order)
        self.orbital_order: np.ndarray = np.array(orbital_order)
        """Order of the orbitals"""
        self.s1_entropy: np.ndarray
        """sinlge orbital entropy"""
        self.s2_entropy: np.ndarray
        """two orbital entropy"""
        self.mutual_information: np.ndarray
        """mutual information"""

    def _get_list_from_fiedler_string(self, fiedler_string: str) -> List[int]:
        return [int(i) for i in fiedler_string.split(",")]

    def __mat_measurement(self, label_vector: List[Tuple[int, int]], value_vector: np.ndarray) -> np.ndarray:
        """Build a matrix from a label_vector and a value_vector.

        Parameters
        ----------
        label_vector : List[Tuple[int]]
            contains the converted labels from a measurement
        value_vector : np.ndarray
            contains the values from a measurement

        Returns
        -------
        matrix : np.ndarray
            the full, ordered matrix of the corresponding measurement
        """
        if len(label_vector) != len(value_vector):
            raise AssertionError("problem with vector length in 'matMeasurement")
        matrix = np.zeros((self.L, self.L))
        for i, label in enumerate(label_vector):
            matrix[label[0], label[1]] = value_vector[i]
            matrix[label[1], label[0]] = value_vector[i]
        return matrix

    def __mat_merge_transpose(
            self, label_vector_1: List[Tuple[int, int]], value_vector_1: np.ndarray,
            label_vector_2: List[Tuple[int, int]], value_vector_2: np.ndarray
    ) -> np.ndarray:
        """Build merged matrix from two measurements and corresponding label
        vectors.

        Parameters
        ----------
        label_vector_1 : List[Tuple[int]]
            converted labels for the first measurement
        value_vector_1 : np.ndarray
            contains the values from the first measurement
        label_vector_2 : List[Tuple[int]]
            converted labels for the second measurement
        value_vector_2 : np.ndarray
            contains the values from the second measurement

        Returns
        -------
        matrix : np.ndarray
            the full, ordered, merged matrix of the measurements
        """
        if (len(label_vector_1) != len(value_vector_1)) or (
            len(label_vector_2) != len(value_vector_2)
        ):
            raise AssertionError("problem with vector length in 'matMergeTranspose")
        matrix = np.zeros((self.L, self.L))
        for i, label in enumerate(label_vector_1):
            matrix[label[0], label[1]] = value_vector_1[i]
        for i, label in enumerate(label_vector_2):
            matrix[label[1], label[0]] = value_vector_2[i]
        return matrix

    def get_measurements(self, dmrg_obj) -> Datasets:
        """Read measurements form the QCMaquis HDF5 file.

        Parameters
        ----------
        dmrg_obj :
            The base qcmaquis object

        Returns
        -------
        data : Datasets
            contains all measurements as ordered vectors/matrices
        """
        data = Datasets()

        data.nup = dmrg_obj.getMeasurement("Nup")[1]
        data.ndown = dmrg_obj.getMeasurement("Ndown")[1]
        data.nupdown = dmrg_obj.getMeasurement("Nupdown")[1]

        dm_up = dmrg_obj.getMeasurement("dm_up")
        dm_down = dmrg_obj.getMeasurement("dm_down")
        nupnup = dmrg_obj.getMeasurement("nupnup")
        nupndown = dmrg_obj.getMeasurement("nupndown")
        ndownnup = dmrg_obj.getMeasurement("ndownnup")
        ndownndown = dmrg_obj.getMeasurement("ndownndown")
        doccdocc = dmrg_obj.getMeasurement("doccdocc")
        transfer_up_while_down = dmrg_obj.getMeasurement("transfer_up_while_down")
        transfer_down_while_up = dmrg_obj.getMeasurement("transfer_down_while_up")
        transfer_up_while_down_at_2 = dmrg_obj.getMeasurement("transfer_up_while_down_at_2")
        transfer_up_while_down_at_1 = dmrg_obj.getMeasurement("transfer_up_while_down_at_1")
        transfer_down_while_up_at_2 = dmrg_obj.getMeasurement("transfer_down_while_up_at_2")
        transfer_down_while_up_at_1 = dmrg_obj.getMeasurement("transfer_down_while_up_at_1")
        transfer_pair = dmrg_obj.getMeasurement("transfer_pair")
        spinflip = dmrg_obj.getMeasurement("spinflip")
        nupdocc = dmrg_obj.getMeasurement("nupdocc")
        ndowndocc = dmrg_obj.getMeasurement("ndowndocc")
        doccnup = dmrg_obj.getMeasurement("doccnup")
        doccndown = dmrg_obj.getMeasurement("doccndown")

        # TODO remove this and directly put it into mat merge
        dmup_labels = dm_up[0]
        dmup = dm_up[1]
        dmdown_labels = dm_down[0]
        dmdown = dm_down[1]
        nupnup_labels = nupnup[0]
        nupnup = nupnup[1]
        nupndown_labels = nupndown[0]
        nupndown = nupndown[1]
        ndownnup_labels = ndownnup[0]
        ndownnup = ndownnup[1]
        ndownndown_labels = ndownndown[0]
        ndownndown = ndownndown[1]
        doccdocc_labels = doccdocc[0]
        doccdocc = doccdocc[1]
        tuwd_labels = transfer_up_while_down[0]
        tuwd = transfer_up_while_down[1]
        tdwu_labels = transfer_down_while_up[0]
        tdwu = transfer_down_while_up[1]
        transferpair_labels = transfer_pair[0]
        transferpair = transfer_pair[1]
        spinflip_labels = spinflip[0]
        spinflip = spinflip[1]
        tuwd_at2_labels = transfer_up_while_down_at_2[0]
        tuwd_at2 = transfer_up_while_down_at_2[1]
        tuwd_at1_labels = transfer_up_while_down_at_1[0]
        tuwd_at1 = transfer_up_while_down_at_1[1]
        tdwu_at2_labels = transfer_down_while_up_at_2[0]
        tdwu_at2 = transfer_down_while_up_at_2[1]
        tdwu_at1_labels = transfer_down_while_up_at_1[0]
        tdwu_at1 = transfer_down_while_up_at_1[1]
        nupdocc_labels = nupdocc[0]
        nupdocc = nupdocc[1]
        doccnup_labels = doccnup[0]
        doccnup = doccnup[1]
        ndowndocc_labels = ndowndocc[0]
        ndowndocc = ndowndocc[1]
        doccndown_labels = doccndown[0]
        doccndown = doccndown[1]

        # build matrices from labels and corresponding vectors
        data.dmUpSym = self.__mat_measurement(dmup_labels, dmup)
        data.dmDownSym = self.__mat_measurement(dmdown_labels, dmdown)
        data.nupNupSym = self.__mat_measurement(nupnup_labels, nupnup)
        data.nupNdownSym = self.__mat_measurement(nupndown_labels, nupndown)
        data.ndownNupSym = self.__mat_measurement(ndownnup_labels, ndownnup)
        data.ndownNdownSym = self.__mat_measurement(ndownndown_labels, ndownndown)
        data.doccDoccSym = self.__mat_measurement(doccdocc_labels, doccdocc)
        data.transferUpWhileDownSym = self.__mat_measurement(tuwd_labels, tuwd)
        data.transferDownWhileUpSym = self.__mat_measurement(tdwu_labels, tdwu)
        data.transferPairSym = self.__mat_measurement(transferpair_labels, transferpair)
        data.spinflipSym = self.__mat_measurement(spinflip_labels, spinflip)
        data.transUpDown2Asym = self.__mat_merge_transpose(tuwd_at2_labels, tuwd_at2, tuwd_at1_labels, tuwd_at1)
        data.transUpDown1Asym = self.__mat_merge_transpose(tuwd_at1_labels, tuwd_at1, tuwd_at2_labels, tuwd_at2)
        data.transDownUp2Asym = self.__mat_merge_transpose(tdwu_at2_labels, tdwu_at2, tdwu_at1_labels, tdwu_at1)
        data.transDownUp1Asym = self.__mat_merge_transpose(tdwu_at1_labels, tdwu_at1, tdwu_at2_labels, tdwu_at2)
        data.nupDoccAsym = self.__mat_merge_transpose(nupdocc_labels, nupdocc, doccnup_labels, doccnup)
        data.doccNupAsym = self.__mat_merge_transpose(doccnup_labels, doccnup, nupdocc_labels, nupdocc)
        data.ndownDoccAsym = self.__mat_merge_transpose(ndowndocc_labels, ndowndocc, doccndown_labels, doccndown)
        data.doccNdownAsym = self.__mat_merge_transpose(doccndown_labels, doccndown, ndowndocc_labels, ndowndocc)
        return data

    def make_one_ordm(self, data: Datasets) -> np.ndarray:
        """
        Evaluate one orbital RDM from qcmaquis measurements.

        Parameters
        ----------
        data : Datasets
            the dataset with all measurements

        Returns
        -------
        one_ordm : np.ndarray
            The one ORBITAL reduced density matrix
        """
        one_ordm = np.zeros((self.L, 4))
        for i in range(self.L):
            one_ordm[i, 0] = data.nup[i] - data.nupdown[i]
            one_ordm[i, 1] = data.ndown[i] - data.nupdown[i]
            one_ordm[i, 2] = 1 - data.nup[i] - data.ndown[i] + data.nupdown[i]
            one_ordm[i, 3] = data.nupdown[i]
        return one_ordm

    def make_two_ordm(self, data: Datasets) -> np.ndarray:
        """
        Evaluate two orbital RDM from qcmaquis measurements.

        Parameters
        ----------
        data : Datasets
            the dataset with all measurements

        Returns
        -------
        two_ordm : np.ndarray
            The two ORBITAL reduced density matrix
        """
        two_ordm = np.zeros((self.L, self.L, 16, 16))

        # p, q is the orbital index from whole space
        # pylint: disable=C0103
        # pylint: disable=W503
        for p in range(self.L):
            for q in range(p + 1, self.L):
                two_ordm[p, q, 0, 0] = (
                    1.0
                    + data.nupdown[p]
                    + data.nupdown[q]
                    + data.doccDoccSym[p, q]
                    - data.ndown[p]
                    - data.ndownDoccAsym[p, q]
                    - data.ndown[q]
                    - data.doccNdownAsym[p, q]
                    + data.ndownNdownSym[p, q]
                    - data.nup[p]
                    - data.nupDoccAsym[p, q]
                    + data.nupNdownSym[p, q]
                    - data.nup[q]
                    - data.doccNupAsym[p, q]
                    + data.ndownNupSym[p, q]
                    + data.nupNupSym[p, q]
                )
                two_ordm[p, q, 1, 1] = (
                    -data.nupdown[p]
                    - data.doccDoccSym[p, q]
                    + data.ndown[p]
                    + data.ndownDoccAsym[p, q]
                    + data.doccNdownAsym[p, q]
                    - data.ndownNdownSym[p, q]
                    + data.doccNupAsym[p, q]
                    - data.ndownNupSym[p, q]
                )
                two_ordm[p, q, 2, 2] = (
                    -data.nupdown[p]
                    - data.doccDoccSym[p, q]
                    + data.doccNdownAsym[p, q]
                    + data.nup[p]
                    + data.nupDoccAsym[p, q]
                    - data.nupNdownSym[p, q]
                    + data.doccNupAsym[p, q]
                    - data.nupNupSym[p, q]
                )
                two_ordm[p, q, 3, 3] = (
                    data.nupdown[p]
                    - data.doccNdownAsym[p, q]
                    - data.doccNupAsym[p, q]
                    + data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 4, 4] = (
                    -data.nupdown[q]
                    - data.doccDoccSym[p, q]
                    + data.ndownDoccAsym[p, q]
                    + data.ndown[q]
                    + data.doccNdownAsym[p, q]
                    - data.ndownNdownSym[p, q]
                    + data.nupDoccAsym[p, q]
                    - data.nupNdownSym[p, q]
                )
                two_ordm[p, q, 5, 5] = (
                    data.ndownNdownSym[p, q]
                    - data.ndownDoccAsym[p, q]
                    - data.doccNdownAsym[p, q]
                    + data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 6, 6] = (
                    data.nupNdownSym[p, q]
                    - data.doccNdownAsym[p, q]
                    - data.nupDoccAsym[p, q]
                    + data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 7, 7] = (
                    data.doccNdownAsym[p, q]
                    - data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 8, 8] = (
                    -data.nupdown[q]
                    - data.doccDoccSym[p, q]
                    + data.ndownDoccAsym[p, q]
                    + data.nupDoccAsym[p, q]
                    + data.nup[q]
                    + data.doccNupAsym[p, q]
                    - data.ndownNupSym[p, q]
                    - data.nupNupSym[p, q]
                )
                two_ordm[p, q, 9, 9] = (
                    data.ndownNupSym[p, q]
                    - data.ndownDoccAsym[p, q]
                    - data.doccNupAsym[p, q]
                    + data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 10, 10] = (
                    data.nupNupSym[p, q]
                    - data.nupDoccAsym[p, q]
                    - data.doccNupAsym[p, q]
                    + data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 11, 11] = (
                    data.doccNupAsym[p, q]
                    - data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 12, 12] = (
                    data.nupdown[q]
                    - data.nupDoccAsym[p, q]
                    - data.ndownDoccAsym[p, q]
                    + data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 13, 13] = (
                    data.ndownDoccAsym[p, q]
                    - data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 14, 14] = (
                    data.nupDoccAsym[p, q]
                    - data.doccDoccSym[p, q]
                )
                two_ordm[p, q, 15, 15] = data.doccDoccSym[p, q]
                two_ordm[p, q, 1, 4] = (
                    data.dmDownSym[p, q]
                    - data.transDownUp1Asym[p, q]
                    - data.transDownUp2Asym[p, q]
                    + data.transferDownWhileUpSym[p, q]
                )
                two_ordm[p, q, 4, 1] = two_ordm[p, q, 1, 4]
                two_ordm[p, q, 2, 8] = (
                    data.dmUpSym[p, q]
                    - data.transUpDown1Asym[p, q]
                    - data.transUpDown2Asym[p, q]
                    + data.transferUpWhileDownSym[p, q]
                )
                two_ordm[p, q, 8, 2] = two_ordm[p, q, 2, 8]
                two_ordm[p, q, 3, 6] = (
                    data.transDownUp1Asym[p, q]
                    - data.transferDownWhileUpSym[p, q]
                )
                two_ordm[p, q, 6, 3] = two_ordm[p, q, 3, 6]
                two_ordm[p, q, 3, 9] = (
                    -data.transUpDown1Asym[p, q]
                    + data.transferUpWhileDownSym[p, q]
                )
                two_ordm[p, q, 9, 3] = two_ordm[p, q, 3, 9]
                two_ordm[p, q, 6, 9] = data.spinflipSym[p, q]
                two_ordm[p, q, 9, 6] = two_ordm[p, q, 6, 9]
                two_ordm[p, q, 3, 12] = data.transferPairSym[p, q]
                two_ordm[p, q, 12, 3] = two_ordm[p, q, 3, 12]
                two_ordm[p, q, 6, 12] = -(
                    -data.transUpDown2Asym[p, q]
                    + data.transferUpWhileDownSym[p, q]
                )
                two_ordm[p, q, 12, 6] = two_ordm[p, q, 6, 12]
                two_ordm[p, q, 9, 12] = -(
                    data.transDownUp2Asym[p, q]
                    - data.transferDownWhileUpSym[p, q]
                )
                two_ordm[p, q, 12, 9] = two_ordm[p, q, 9, 12]
                two_ordm[
                    p, q, 7, 13
                ] = -data.transferUpWhileDownSym[p, q]
                two_ordm[p, q, 13, 7] = two_ordm[p, q, 7, 13]
                two_ordm[
                    p, q, 11, 14
                ] = -data.transferDownWhileUpSym[p, q]
                two_ordm[p, q, 14, 11] = two_ordm[p, q, 11, 14]
        return two_ordm

    def make_s1(self, one_ordm: np.ndarray):
        """Evaluate single orbital entropy from one orbital RDM.

        Parameters
        ----------
        one_ordm : np.ndarray
            one orbital RDM
        """
        self.s1_entropy = np.zeros((self.L))
        for site in range(self.L):
            s1_entropy = 0
            for alpha in range(len(one_ordm[site])):
                eigenvalue = one_ordm[site, alpha]
                if eigenvalue > 0:
                    s1_entropy = s1_entropy - eigenvalue * np.log(eigenvalue)
            self.s1_entropy[site] = s1_entropy

    def make_s2(self, two_ordm: np.ndarray):
        """Evaluate two orbital entropy from two orbital RDM.

        Parameters
        ----------
        two_ordm : np.ndarray
            two orbital RDM
        """
        self.s2_entropy = np.zeros((self.L, self.L))
        for site_1 in range(self.L):
            for site_2 in range(site_1 + 1, self.L):
                sub_matrix = two_ordm[site_1][site_2][:][:]
                eigenvalue, _ = np.linalg.eig(sub_matrix)  # type: ignore[attr-defined]
                s2_entropy = 0
                for alpha in range(16):
                    if eigenvalue[alpha] > 0:
                        s2_entropy = s2_entropy - (eigenvalue[alpha] * np.log(eigenvalue[alpha]))
                self.s2_entropy[site_1, site_2] = s2_entropy.real
                self.s2_entropy[site_2, site_1] = s2_entropy.real

    def make_mutual_information(self):
        """Evaluate the mutual information.

        Raises
        ------
        AttributeError
            if s1 and s2 entropies are not evaluated
        """
        if not (self.s1_entropy.size and self.s2_entropy.size):
            raise AttributeError("Evaluate s2 and s2 entropy before evaluating mutual information")
        self.mutual_information = np.zeros((self.L, self.L))
        for site_1 in range(self.L):
            for site_2 in range(site_1 + 1, self.L):
                self.mutual_information[site_1, site_2] = 0.5 * (
                    self.s1_entropy[site_1] + self.s1_entropy[site_2] - self.s2_entropy[site_1, site_2]
                )
                self.mutual_information[site_2, site_1] = self.mutual_information[site_1, site_2]

    def make_diagnostics(self, dmrg_obj):
        """Generate orbital RDMs from QCMaquis output file and evaluate entropies.

        Raises
        ------
        AttributeError
            if hdf5 converter has not read the output yet
        """
        data = self.get_measurements(dmrg_obj)

        self.make_s1(self.make_one_ordm(data))
        self.make_s2(self.make_two_ordm(data))
        self.make_mutual_information()

        # handle fiedler ordering
        sort_key = np.argsort(self.orbital_order)
        self.s1_entropy = self.s1_entropy[sort_key]
        self.s2_entropy = self.s2_entropy[sort_key][:, sort_key]
        self.mutual_information = self.mutual_information[sort_key][:, sort_key]
