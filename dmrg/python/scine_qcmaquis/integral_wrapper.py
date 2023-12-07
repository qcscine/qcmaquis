import re
from dataclasses import dataclass
from enum import Enum
from io import TextIOWrapper
from typing import Dict, List, Tuple, Union

import numpy as np
from _dmrg import ComplexTCIntegralMap, IntegralMap, TCIntegralMap


class IntegralType(Enum):
    """Set type for reading integrals.

    Attributes
    ----------
    CONVENTIONAL
    TRANSCORRLEATED
    """
    CONVENTIONAL = "conventional"
    """Enable 8-fold symmetry in the two-body integrals."""
    TRANSCORRLEATED = "transcorrelated"
    """Enable 2-fold symmetry in the two-body integrals."""


class IntegralMapWrapper:
    """Wrapper for IntegralMap Python Interface.

    Attributes
    ----------
    _integral_map : IntegralMap
        the QcMaquis integral map
    _type : IntegralType, default = IntegralType.CONVENTIONAL
        the type of integrals, responsible for symmetries
    _parser : IntegralsParser
        handler to parse integrals
    """
    # TODO: add __slots__

    def __init__(self, ):
        """Constructor."""
        self._integral_map = IntegralMap()
        """The binded integral map."""
        self._type = IntegralType.CONVENTIONAL
        """The integral type."""
        self._parser = IntegralsParser()
        """Handle integral parsing."""

    def set_type(self, integral_type: IntegralType):
        """Set the integral type.

        Parameters
        ----------
        integral_type : IntegralType
            the type of integrals, responsible for symmetries
        """
        if self._type != integral_type:
            self._type = integral_type
            if integral_type == IntegralType.CONVENTIONAL:
                self._integral_map = IntegralMap()
            elif integral_type == IntegralType.TRANSCORRLEATED:
                # TODO: Check if this has to be complex
                self._integral_map = ComplexTCIntegralMap()
            else:
                raise NotImplementedError(f"IntegralType: <{integral_type}> is unavailable atm.")

    # TODO
    def fill_from_fcidump(self, fcidump: str):
        """Fill IntegralMap from an FCIDUMP."""
        raise NotImplementedError("integrals from fcidump are not yet supported")

    def fill_from_pyscf(self, core_value: float, one_body: np.ndarray, two_body: np.ndarray, norb: int):
        """Fill IntegralMap from a PySCF wavefunction.

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
        self._parser.set_core(core_value)
        self._parser.parse_one_body(one_body, norb)
        self._parser.parse_two_body(two_body, norb)
        self.update_from_parsing()

    def get(self) -> IntegralMap:
        """Get IntegralMap Python Interface.

        Return
        ------
        self._integral_map : IntegralMap
            the filled integral map
        """
        return self._integral_map

    def set(self, integral_map: Union[IntegralMap, TCIntegralMap]):
        """Replace integral_map with new map.

        Parameters
        ----------
        integral_map : IntegralMap
            the new, filled integral_map
        """
        self._integral_map = integral_map

    def update_from_parsing(self, transcorrelated=False):
        """Smth."""
        if transcorrelated is False:
            self._integral_map = IntegralMap()
            for key in self._parser.get_unique_indices():
                new_key = (key[0], key[1], key[2], key[3])
                self._integral_map.set(new_key, float(self._parser.get_integral_value(key)))
        else:
            self._integral_map = TCIntegralMap()
            for key in self._parser.get_unique_indices():
                new_key = (key[0], key[1], key[2], key[3], key[4], key[5])
                self._integral_map.set(new_key, float(self._parser.get_integral_value(key)))


class IntegralsParser:
    """Parse Integrals.

    Attributes
    ----------
    _unique_term : Dict(Tuple[int, int, int, int], float)
        Store integrals.
    _integral_utils : IntegralUtils
        Handler for integral symmtries
    fcidump_values : FcidumpValues
        values from the FCIDUMP header

    Note
    ----
    When storing the integrals, the key is a tuple of indices with a value
    corresponding to the integral.
    """

    # @dataclass
    class FcidumpValues:
        """Store data values from FCIDUMP header.

        Attributes
        ----------
        norb: int
            number of orbitals
        nelec: int
            number of electrons
        ms2: int
            spin polarization
        orbsym: List[int]
            symmetry of each orbital
        isym: int
            point group
        transcorrelated : bool, default = False
            integrals are transcorrelated
        unrestricted : bool, default = False
            integrals are in spin orbitals
        """

        __slots__ = ("norb", "nelec", "ms2", "orbsym", "isym", "transcorrelated", "unrestricted")

        # TODO: use regex for fcidump parsing

        def __init__(self):
            self.norb: int = 0
            """Number of orbitals."""
            self.nelec: int
            """Number of electrons."""
            self.ms2: int
            """Spin polarization."""
            self.orbsym: List[int]
            """Symmetry of each orbital."""
            self.isym: int
            """Point group."""
            self.transcorrelated = False
            """Integrals are transcorrelated."""
            self.unrestricted = False
            """Integrals are in spin orbital basis."""

    # TODO add __slots__
    # TODO interface to Hamiltonian in CC code
    def __init__(self):
        """Constructor."""
        self._unique_term: Dict(Tuple[int, int, int, int], float) = {}
        """Store integrals."""
        self._integral_utils = IntegralUtils()
        """Handler for symmtries."""
        self.fcidump_values = self.FcidumpValues()
        """Store FCIDUMP header data."""

    def get_unique_indices(self) -> Tuple[int, int, int, int]:
        """Get unique integrals."""
        return self._unique_term.keys()

    def get_integral_value(self, key: Tuple[int, int, int, int]) -> float:
        """Get integral_value."""
        return self._unique_term[key]

    def parse_fcidump(self, fcidump: str):
        """Parse fcidump."""
        fcidump = open(fcidump, "r")
        self._parse_fcidump_header(fcidump)
        self._parse_fcidump_body(fcidump)

    def _parse_fcidump_header(self, file: TextIOWrapper):
        parse_header = False
        line = file.readline().lower()
        while line:
            if "&fci" in line:
                parse_header = True

            if parse_header is True:
                self._find_keywords(line)

            if "&end" in line:
                parse_header = False
                break

            line = file.readline().lower()

    def _find_keywords(self, line: str):
        if "norb" in line:
            self.fcidump_values.norb = int(re.search(r"norb\s*=\s*(\d+)", line).group(1))
        if "nelec" in line:
            self.fcidump_values.nelec = int(re.search(r"nelec\s*=\s*(\d+)", line).group(1))
        if "ms2" in line:
            self.fcidump_values.ms2 = int(re.search(r"ms2\s*=\s*(\d+)", line).group(1))
        if "isym" in line:
            self.fcidump_values.isym = int(re.search(r"isym\s*=\s*(\d+)", line).group(1))
        if "transcorrelated" in line:
            self.fcidump_values.transcorrelated = True
        if "unrestricted" in line:
            self.fcidump_values.unrestricted = True
        if "orbsym" in line:
            self.fcidump_values.orbsym = [
                int(x) for x in re.search(r"orbsym\s*=\s*([\d+,]+)\s*,", line)
                .group(1).split(",") if x.strip().isdigit()
            ]

    def _parse_fcidump_body(self, file: TextIOWrapper):
        line = file.readline()
        while line:
            self._add_term(line)
            line = file.readline()

    def _add_term(self, line: str):
        line = line.split()
        value = float(line[0])
        # pylint: disable=invalid-name
        p = int(line[1])
        q = int(line[2])
        r = int(line[3])
        s = int(line[4])
        # pylint: enable=invalid-name
        # chemist -> physics
        #       1  2  1  2
        if p != 0 and q != 0 and r != 0 and s != 0:
            term = (p, r, q, s)
            if self._is_unique((p, q, r, s)):
                self._unique_term[term] = value
        else:
            term = (p, q, 0, 0)
            self._unique_term[term] = value

    def _is_unique(self, indices: List[int]) -> bool:
        tmp = self._integral_utils.get_symmetric_indices(indices, "eight")
        for term in tmp:
            if tuple(term) in self._unique_term:
                return False
        return True

    def _is_unique_one_body(self, indices: Tuple[int, int, int, int]) -> bool:
        if indices[2] != 0 or indices[3] != 0:
            raise ValueError("wrong indices for one body")
        symm_indices = (indices[1], indices[0], 0, 0)
        tmp = [indices, symm_indices]
        # if any([term in self._unique_term for term in tmp]):
        if any(term in self._unique_term for term in tmp):
            return False
        return True

    def set_core(self, core_value):
        """From Pyscf."""
        self._unique_term[(0, 0, 0, 0)] = core_value

    def parse_one_body(self, one_body_ints, norb):
        """From Pyscf."""
        for i in range(norb):
            for j in range(norb):
                if self._is_unique_one_body((i + 1, j + 1, 0, 0)):
                    self._unique_term[(i + 1, j + 1, 0, 0)] = one_body_ints[i, j]

    def parse_two_body(self, two_body_ints, norb):
        """From Pyscf."""
        for i in range(1, norb + 1):
            for j in range(1, norb + 1):
                for k in range(1, norb + 1):
                    for l in range(1, norb + 1):
                        if self._is_unique([i, j, k, l]):
                            self._unique_term[(i, j, k, l)] = two_body_ints[i - 1, j - 1, k - 1, l - 1]


class IntegralNotation(Enum):
    """Indicate integral notation."""
    PHYSICS = 1
    CHEMISTRY = 2


class IntegralUtils:
    """Converte Integrals from extern to qcmaquis notation."""
    # TODO add __slots__

    def __init__(self):
        self._notation = IntegralNotation.CHEMISTRY

    def set_notation(self, notation: IntegralNotation):
        """Set notation."""
        self._notation = notation

    def _permute_particle_block(
        self,
        result_list: List[List[int]],
        block_index_1: int,
        block_index_2: int,
    ) -> List[List[int]]:
        """Generate symmetry by permuting two orbital blocks,
        e.g. in chemistry notation (pq | rs | tu)=(rs | pq | tu).

        Parameters
        ----------
        result_list : List[List[int]]
            List containing List with indices to permute
        block_index_1 : int
            first index of the first block to permute
        block_index_2 : int
            first index of the second block to permute with

        Returns
        -------
        List[tuple(int)]
            appended List containing new list of indices
        """
        for index_list in result_list:
            new_index_list = index_list.copy()

            blub = new_index_list[block_index_1]
            new_index_list[block_index_1] = new_index_list[block_index_2]
            new_index_list[block_index_2] = blub

            blub = new_index_list[block_index_1 + 1]
            new_index_list[block_index_1 + 1] = new_index_list[block_index_2 + 1]
            new_index_list[block_index_2 + 1] = blub

            if new_index_list not in result_list:
                result_list.append(new_index_list)
        return result_list

    # Copied from full cc
    def _permute_same_particle(
        self,
        result_list: List[List[int]],
        particle_index_1: int,
        particle_index_2: int,
    ) -> List[List[int]]:
        """Generate symmetry by permuting two orbitals with same particle,
        e.g. in chemistry notation (pq | rs | tu)=(pq | sr| tu).

        Parameters
        ----------
        result_list : List[List[int]]
            List containing List with indices to permute
        block_index_1 : int
            index of the first particle to permute
        block_index_2 : int
            index of the second particle to permute with

        Returns
        -------
        List[tuple(int)]
            appended List containing new list of indices
        """
        for index_list in result_list:
            new_index_list = index_list.copy()

            blub = new_index_list[particle_index_1]
            new_index_list[particle_index_1] = new_index_list[particle_index_2]
            new_index_list[particle_index_2] = blub

            if new_index_list not in result_list:
                result_list.append(new_index_list)
        return result_list

    def get_symmetric_indices(
        self,
        index_list: List[int],
        symmetry: str,
    ) -> List[List[int]]:
        """Get all list of list of indices corresponding to same integral
        based on provided symmetry.

        Parameters
        ----------
        index_list : List[int]
            List with indices to permute
        symmetry : str
            Viable symmetries are: 'one', 'two', 'four', 'eight', 'fourtyeight'

        Returns
        -------
        List[List[int]]
            a List with Lists of indices
        """
        if self._notation is not IntegralNotation.CHEMISTRY:
            raise ValueError("Only Chemistry notation supported yet")

        result_list = []
        result_list.append(index_list)
        if symmetry == "one":
            pass
        # transcorr
        elif symmetry == "two":
            # (pq | rs) = (rs| pq)
            result_list = self._permute_particle_block(result_list, 0, 2)
        elif symmetry == "four":
            # (pq | rs) = (qp | rs) = (qp | sr) = (pq | sr)
            result_list = self._permute_same_particle(result_list, 0, 1)
            result_list = self._permute_same_particle(result_list, 2, 3)
        # Conventional
        elif symmetry == "eight":
            result_list = self._permute_particle_block(result_list, 0, 2)
            result_list = self._permute_same_particle(result_list, 0, 1)
            result_list = self._permute_same_particle(result_list, 2, 3)
        # transcorr
        elif symmetry == "fourtyeight":
            result_list = self._permute_particle_block(result_list, 0, 2)
            result_list = self._permute_particle_block(result_list, 2, 4)
            result_list = self._permute_particle_block(result_list, 4, 0)
            result_list = self._permute_same_particle(result_list, 0, 1)
            result_list = self._permute_same_particle(result_list, 2, 3)
            result_list = self._permute_same_particle(result_list, 4, 5)
        else:
            raise ValueError(f"{symmetry} is not supported!")
        return result_list
