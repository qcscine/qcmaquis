import os
import shutil

from typing import Any, List, Optional, Tuple, Union
from itertools import permutations

import numpy as np
import pyscf

from scine_qcmaquis import MaquisDmrg
from scine_qcmaquis.dmrg_wrapper import ParametersWrapper


class QcMaquis:
    """QCMaquis interface for pyscf.

    This interface allows qcmaquis to act as a fcisolver in pyscf.mcscf methods.

    TODO
    ----
    For MRPT higher order rdms are required, and the corresponding functions have
    to be implemented here, in order to work automatically by injection of this
    interface into pyscf.

    Attributes
    ----------
    mol : Any
    verbose : int
    orbsym : List[int]
    method : str
    n_states : Optional[int]
    fiedler : bool, default = False
    log : Any
    dmrg : Optional[MaquisDmrg]
    file_path : str
    measure_entropies : bool
    checkpoint_name : str
    results_name : str
    energy_threshold : float
    nsweeps : int
    bond_dim : int
    initial_truncation_thresh : float
    final_truncation_thresh : float
    """

    # pylint: disable=W0613
    def __init__(
        self,
        mol: Any,
        verbose: Optional[int] = None,
        fiedler: bool = False,
        **kwargs: Any,
    ) -> None:
        """Initialize interface.

        Parameters
        ----------
        mol : Any
            pyscf molecule
        **kwargs : Any
            other potential arguments, but unused
        """
        # Pyscf stuff
        self.mol: Any = mol
        """Pyscf molecule"""
        self.verbose = verbose if verbose is not None else mol.verbose
        """Output verbosity"""
        self.orbsym: List[int] = []
        """Orbital symmetries"""
        self.log: Any = None
        """Pyscf logger."""

        # DMRG stuff
        self.method: str = "conventional"
        """DMRG method, e.g. conventional, ..."""
        self.n_states: Optional[int] = None
        """Number of states."""
        self.fiedler: bool = fiedler
        """Enable fiedler ordering"""

        # QCMaquis related stuff
        self.dmrg: Optional[MaquisDmrg] = None
        """QCMaquis python interface, will be set automatically."""
        self.file_path: str = os.getcwd()
        """Path for qcmaquis dumps."""
        self.parameters = ParametersWrapper()
        """QCMaquis parameter handler"""

        # dont change this
        self._dmrgscf = False
        """Flag to indicate DMRGSCF, important for checkpoint files"""
        self._dmrgscf_checkpoint_name = "checkpoint_DMRGSCF.h5"
        """Name of DMRGSCF checkpoint file"""
        self._dmrgscf_results_name = "results_DMRGSCF.h5"
        """Name of DMRGSCF results file"""

        # self.measure_entropies: bool = False
        # """Enable measurements for s1, s2 and mut inf."""
        # self.orbital_optimization: bool = True
        # """Enable measurements for 1rdm and 2rdm"""
        # self.checkpoint_name: str = "qcmaquis_checkpoint.h5"
        # """Name of qcmaquis checkpoint."""
        # self.results_name: str = "qcmaquis_result_file.h5"
        # """Name of qcmaquis results file."""

        # DMRG parameters
        # self.energy_threshold: float = 1e-6
        # """Energy threshold."""
        # self.nsweeps: int = 100
        # """Max number of sweeps."""
        # self.bond_dim: int = 250
        # """Bond dimension."""
        # self.initial_truncation_thresh: float = 1e-16  # 1e-6
        # """Truncation threshold for first sweeps"""
        # self.final_truncation_thresh: float = 1e-16  # 1e-10
        # """Truncation threshold for final sweeps"""

        # Set other parameters directly through qcmaquis parameters
        # self.feast_window = None
        # self.feast_states = None
        # self.nroots = 1
        # self.stdout = mol.stdout

        # Do I need this?
        # self._keys = set(self.__dict__.keys())

    # pylint: enable =W0613

    def dump_flags(self, verbose: Optional[int] = None):
        """Log settings for pyscf logger.

        Parameters
        ----------
        verbose : Optional[int]
            pyscf verbosity
        """
        self.log = pyscf.lib.logger.new_logger(self, verbose)
        self.log.info("************** QcMaquis flags **************")
        self.log.info("method           = %s", str(self.method))
        self.log.info("fiedler ordering = %s", str(self.fiedler))
        if self.file_path:
            self.log.info("checkpoint path  = %s", str(self.file_path))
            self.log.info(
                "checkpoint name  = %s", str(self.parameters._checkpoint_path)
            )
            self.log.info("results name     = %s", str(self.parameters._results_path))
            self.log.info("storage dir      = %s", str(self.parameters._storage_dir))
        else:
            self.log.info("skipping checkpoints")

        if self.n_states:
            self.log.info("Number of state  = %s", str(self.n_states))

    def _gen_rdm_permutations(self, n_particles: int, is_hermitian=True):
        """Generate all equivalent elements for n_particle-RDMs"""
        indices = []
        creators_original = [i for i in range(n_particles)]
        for creators in permutations(creators_original, n_particles):
            annhilators = tuple(
                [(i + n_particles) % (2 * n_particles) for i in creators]
            )
            indices.append(creators + annhilators)
            if is_hermitian:
                indices.append(annhilators + creators)

        return indices

    def _get_rdm1(self, norb: int) -> np.ndarray:
        """Getter for 1 rdm.

        Transforms qcmaquis 1 particle dm into pyscf compatible format.

        QCMaquis format: dm[p,q] = < p^+ q >
        PySCF format:    dm[p,q] = < q^+ p>

        Parameters
        ----------
        norb : int
            number of orbitals

        Returns
        -------
        rdm1 : np.ndarray
            1-particle reduced density matrix
        """
        # in case of feast there is no rdm
        try:
            if self.dmrg is not None:
                if self.verbose > 4:
                    # 2rdm
                    maquis_rdm1 = self.dmrg.get_one_rdm()
                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        # 2rdm
                        maquis_rdm1 = self.dmrg.get_one_rdm()
            else:
                maquis_rdm1 = ([[0 for _ in range(2)]], [0])

        except RuntimeError:
            maquis_rdm1 = ([[0 for _ in range(2)]], [0])

        # convert 1 rdm from qcmaquis to pyscf format
        rdm2 = np.zeros((norb,) * 2)
        for i, vec in enumerate(maquis_rdm1[0]):
            for permutation in self._gen_rdm_permutations(1, True):
                rdm2[tuple(vec[i] for i in permutation)] = maquis_rdm1[1][i]
        rdm2 = rdm2.T
        return rdm2

    def _get_rdm2(self, norb: int) -> np.ndarray:
        """Getter for 2 rdm.

        Transforms qcmaquis 2 particle dm into pyscf compatible format.

        QCMaquis format: dm[p,r,...,s,q] = < p^+ r^+ ... s q >
        PySCF format:    dm[p,q,r,s,...] = < p^+ r^+ ... s q >

        Parameters
        ----------
        norb : int
            number of orbitals

        Returns
        -------
        rdm2 : np.ndarray
            2-particle reduced density matrix
        """
        # in case of feast there is no rdm
        try:
            if self.dmrg is not None:
                if self.verbose > 4:
                    # 2rdm
                    maquis_rdm2 = self.dmrg.get_two_rdm()
                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        # 2rdm
                        maquis_rdm2 = self.dmrg.get_two_rdm()
            else:
                maquis_rdm2 = ([[0 for _ in range(4)]], [0])

        except RuntimeError:
            maquis_rdm2 = ([[0 for _ in range(4)]], [0])

        # convert 2 rdm from qcmaquis to pyscf format
        rdm2 = np.zeros((norb,) * 4)
        for i, vec in enumerate(maquis_rdm2[0]):
            for permutation in self._gen_rdm_permutations(2, True):
                rdm2[tuple(vec[i] for i in permutation)] = maquis_rdm2[1][i]
        rdm2 = rdm2.transpose(0, 3, 1, 2)
        return rdm2

    def _get_rdm3(self, norb: int) -> np.ndarray:
        """Getter for 3 rdm.

        Transforms qcmaquis 3 particle dm into pyscf compatible format.

        QCMaquis format: dm[p,r,...,s,q] = < p^+ r^+ ... s q >
        PySCF format:    dm[p,q,r,s,...] = < p^+ r^+ ... s q >

        Parameters
        ----------
        norb : int
            number of orbitals

        Returns
        -------
        rdm3 : np.ndarray
            3-particle reduced density matrix
        """
        # in case of feast there is no rdm
        try:
            if self.dmrg is not None:
                if self.verbose > 4:
                    # 3rdm
                    maquis_rdm3 = self.dmrg.get_three_rdm()
                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        # 3rdm
                        maquis_rdm3 = self.dmrg.get_three_rdm()
            else:
                maquis_rdm3 = ([[0 for _ in range(6)]], [0])

        except RuntimeError:
            maquis_rdm3 = ([[0 for _ in range(6)]], [0])

        # convert 3 rdm from qcmaquis to pyscf format
        rdm3 = np.zeros((norb,) * 6)
        for i, vec in enumerate(maquis_rdm3[0]):
            for permutation in self._gen_rdm_permutations(3, True):
                rdm3[tuple(vec[i] for i in permutation)] = maquis_rdm3[1][i]

        rdm3 = rdm3.transpose(0, 5, 1, 4, 2, 3)
        return rdm3

    def _get_rdm4(self, norb: int) -> np.ndarray:
        """Getter for 4 rdm.

        Transforms qcmaquis 4 particle dm into pyscf compatible format.

        QCMaquis format: dm[p,r,...,s,q] = < p^+ r^+ ... s q >
        PySCF format:    dm[p,q,r,s,...] = < p^+ r^+ ... s q >

        Parameters
        ----------
        norb : int
            number of orbitals

        Returns
        -------
        rdm4 : np.ndarray
            4-particle reduced density matrix
        """
        # in case of feast there is no rdm
        try:
            if self.dmrg is not None:
                if self.verbose > 4:
                    # 4rdm
                    maquis_rdm4 = self.dmrg.get_four_rdm()
                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        # 4rdm
                        maquis_rdm4 = self.dmrg.get_four_rdm()
            else:
                maquis_rdm4 = ([[0 for _ in range(8)]], [0])

        except RuntimeError:
            maquis_rdm4 = ([[0 for _ in range(8)]], [0])

        # convert 4 rdm from qcmaquis to pyscf format
        rdm4 = np.zeros((norb,) * 8)
        for i, vec in enumerate(maquis_rdm4[0]):
            for permutation in self._gen_rdm_permutations(4, True):
                rdm4[tuple(vec[i] for i in permutation)] = maquis_rdm4[1][i]

        rdm4 = rdm4.transpose(0, 7, 1, 6, 2, 5, 3, 4)
        return rdm4

    def _make_dm1234(
        self, norb: int
    ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        dm1 = self._get_rdm1(norb)
        dm2 = self._get_rdm2(norb)
        dm3 = self._get_rdm3(norb)
        dm4 = self._get_rdm4(norb)
        return dm1, dm2, dm3, dm4

    # TODO: enable excited states
    def _set_excited_state_options(self):
        """Set excited state settings."""
        pass
        # Consistency with excited states
        # if self.n_states is None:
        #     pass
        # elif self.n_states >= 1:
        #     self.n_states = None

        # FEAST
        # if self.method.lower() == "feast":
        #     if self.feast_states is not None and self.feast_window is not None:
        #         self.dmrg.set_feast(self.feast_window, self.feast_states)
        #     else:
        #         raise ValueError("For a feast calculation a feast_window and feast_states")

        # ORTHO
        # if self.method.lower() == "ortho":
        #     if self.n_states is None:
        #         raise ValueError("For a feast calculation a feast_window and feast_states")

    def _check_spin(self, nelec: Union[int, Tuple[int, int]]) -> Tuple[int, int]:
        """Check spin based on numpy electron definition.

        Parameters
        ----------
        nelec : Union[int, Tuple[int]]
            either number of electrons or number of alpha and beta electrons

        Returns
        -------
        nelec : int
            total number of electrons
        spin2 : int
            number of unpaired electrons (2S)
        """
        # check spin
        if isinstance(nelec, (int, np.integer)):
            spin2 = 0
        else:
            spin2 = nelec[0] - nelec[1]
            nelec = sum(nelec)

        return nelec, spin2

    def _get_energy(self) -> Union[float, List[float]]:
        """Get energy from interface.

        Returns
        -------
        energy: Union[float, List[float]]
            dmrg energy for each state
        """
        if self.dmrg is not None:
            energy = self.dmrg.get_energy()
            # convert complex to float
            # if self.method.lower() == "feast":
            #     for i, ener in enumerate(energy):
            #         energy[i] = ener.real
            return energy
        raise ValueError("Run DMRG before requesting energies")

    def kernel(self, h1e, eri, norb, nelec, ci0=None, ecore=0, **kwargs):
        """Kernel function for pyscf.

        Create the interface(at the moment integrals are required in order to
        initialize the interface) and run qcmaquis calculation.
        """

        # orbsym is unused
        self.orbsym = np.zeros(norb, np.int32)
        eri = pyscf.ao2mo.restore(1, eri, norb)

        # initialize qcmaquis
        self.dmrg = MaquisDmrg()
        # onerdm is required for pyscf
        self.parameters.set_orbital_optimization()
        self.dmrg.replace_parameters(self.parameters)
        self.dmrg.set_integrals(ecore, h1e, eri, norb)
        # TODO: excited states are not supported yet
        self._set_excited_state_options()
        nelec, spin2 = self._check_spin(nelec)

        # Always run DMRG if no file path is set
        if not self.file_path:
            # make sure no chkpfile and result file is set
            self.parameters.erase("chkpfile", verbose=False)
            self.parameters.erase("resultfile", verbose=False)
            if self.verbose > 4:
                self.dmrg.run(
                    norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler
                )
            else:
                with pyscf.lib.capture_stdout() as stdout:
                    self.dmrg.run(
                        norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler
                    )

        else:
            self._check_file_path()
            cases = self.check_checkpoint_and_results_file()

            #
            if self._dmrgscf:
                self.dmrg.replace_parameters(self.parameters)
                if self.verbose > 4:
                    self.dmrg.run(
                        norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler
                    )
                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        self.dmrg.run(
                            norb,
                            nelec,
                            spin2,
                            n_states=self.n_states,
                            fiedler=self.fiedler,
                        )

            # load from checkpoint (but only DMRGCI)
            elif cases == "DMRGCI":
                self.parameters.erase("integrals")
                self.dmrg.replace_parameters(self.parameters)
                print(
                    f"""No DMRG calculation required. MPS is loaded from: {self.parameters.get_checkpoint_path()}"""
                )
                if self.verbose > 4:
                    self.dmrg.init_dmrg(
                        self.parameters.get_checkpoint_path(), norb, nelec, spin2
                    )

                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        self.dmrg.init_dmrg(
                            self.parameters.get_checkpoint_path(), norb, nelec, spin2
                        )

            # cannot load DMRGSCF from checkpoint yet, so we just run dmrgscf from scratch
            elif cases == "DMRGSCF":
                self._dmrgscf = True
                self._check_file_path()
                cases = self.check_checkpoint_and_results_file()
                self.dmrg.replace_parameters(self.parameters)
                if self.verbose > 4:
                    self.dmrg.run(
                        norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler
                    )
                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        self.dmrg.run(
                            norb,
                            nelec,
                            spin2,
                            n_states=self.n_states,
                            fiedler=self.fiedler,
                        )
            else:
                raise RuntimeError("How did we get here")

        energy = self._get_energy()
        fakewfn_by_rdm2 = self._get_rdm2(norb)

        return energy, fakewfn_by_rdm2

    def make_rdm12(
        self, fakewfn_by_rdm2: np.ndarray, ncas: int, nelec: int, **kwargs: Any
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Make 1rdm and 2rdm.

        TODO
        ----
        Check if this function is required for pyscf, but I think so.
        """
        if not isinstance(nelec, (int, np.integer)):
            nelec = sum(nelec)
        rdm2 = fakewfn_by_rdm2
        rdm1 = np.einsum("ijkk->ij", rdm2) / (nelec - 1)
        return rdm1, rdm2

    def make_rdm1(
        self, fcivec: Any, norb: int, nelec: int, link_index: Any = None, **kwargs: Any
    ) -> np.ndarray:
        """Make 1rdm.

        TODO
        ----
        Check if this function is required for pyscf, but I think so.
        """
        return QcMaquis.make_rdm12(self, fcivec, norb, nelec, **kwargs)[0]

    def get_entropies(self):
        """Getter for chementropies from qcmaquis.

        NOTE
        ----
        Not working yet, since the function in cpp interface is not working yet.
        """
        return self.dmrg.get_entropies()

    def _check_file_path(self):
        """Check if file path is valid and check if checkpoint, result files are in that path."""
        # check filepath
        if self.file_path:
            if not self.file_path.endswith("/"):
                self.file_path += "/"
            if not self.parameters.get_checkpoint_path().startswith(self.file_path):
                self.parameters.set_checkpoint_path(
                    self.file_path + self.parameters.get_checkpoint_path()
                )

            if not self.parameters.get_result_path().startswith(self.file_path):
                self.parameters.set_result_path(
                    self.file_path + self.parameters.get_result_path()
                )

    def check_checkpoint_and_results_file(self):
        """Check if checkpoint and results file exist.

        TODO: rewrite this!
        In case they exist, create new checkpoint/resultsfile.
        If these files exists, they get deleted and newly created.
        Due to this, DMRGSCF works without reading the checkpoint/results from the previous iteration.
        And we can analyze DMRG for the initial orbitals, as well as the optimized orbitals.
        """
        dmrgscf_checkpoint = ""
        for i in self.parameters.get_checkpoint_path().split("/")[:-1]:
            dmrgscf_checkpoint += f"{i}/"
        dmrgscf_checkpoint += self._dmrgscf_checkpoint_name

        dmrgscf_results = ""
        for i in self.parameters.get_result_path().split("/")[:-1]:
            dmrgscf_results += f"{i}/"
        dmrgscf_results += self._dmrgscf_results_name

        # Already did dmrgscf. Read from checkpoint
        if not self._dmrgscf and os.path.exists(dmrgscf_checkpoint):
            self.parameters.set_checkpoint_path(dmrgscf_checkpoint)
            self.parameters.set_result_path(dmrgscf_results)
            return "DMRGSCF"

        if not self._dmrgscf and os.path.exists(self.parameters.get_checkpoint_path()):
            self.parameters.set_checkpoint_path(self.parameters.get_checkpoint_path())
            self.parameters.set_result_path(self.parameters.get_result_path())
            return "DMRGCI"

        self._dmrgscf = True
        # check checkpoint path first dmrgscf iteration this is false
        if os.path.exists(self.parameters.get_checkpoint_path()):
            # erase dmrgscf checkpoint files for new dmrgscf iteration
            if self.parameters.get_checkpoint_path().endswith(
                self._dmrgscf_checkpoint_name
            ):
                shutil.rmtree(self.parameters.get_checkpoint_path())
            # replace dmrgci checkpioint name with dmrgscf checkpoint
            else:
                self.parameters.set_checkpoint_path(dmrgscf_checkpoint)

        # check result path (same logic as for checkpoint)
        if os.path.exists(self.parameters.get_result_path()):
            if self.parameters.get_result_path().endswith(self._dmrgscf_results_name):
                os.remove(self.parameters.get_result_path())
            else:
                self.parameters.set_result_path(dmrgscf_results)
