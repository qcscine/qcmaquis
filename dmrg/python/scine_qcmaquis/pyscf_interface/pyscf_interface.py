import os
import shutil
from typing import Any, List, Optional, Tuple, Union

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
    def __init__(self, mol: Any, **kwargs: Any) -> None:
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
        self.verbose = mol.verbose
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
        self.fiedler: bool = False
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
    # pylint: enable =W0613

    def dump_flags(self, verbose: Optional[int] = None):
        """Log settings for pyscf logger.

        Parameters
        ----------
        verbose : Optional[int]
            pyscf verbosity
        """
        self.log = pyscf.lib.logger.new_logger(self, verbose)
        self.log.info('************** QcMaquis flags **************')
        self.log.info('method           = %s', str(self.method))
        self.log.info('fiedler ordering = %s', str(self.fiedler))
        if self.file_path:
            self.log.info('checkpoint path  = %s', str(self.file_path))
            self.log.info('checkpoint name  = %s', str(self.parameters._checkpoint_path))
            self.log.info('results name     = %s', str(self.parameters._results_path))
            self.log.info('storage dir      = %s', str(self.parameters._storage_dir))
        else:
            self.log.info('skipping checkpoints')

        if self.n_states:
            self.log.info('Number of state  = %s', str(self.n_states))

    def _get_rdm2(self, norb: int) -> np.ndarray:
        """Getter for 2 rdm.

        Transforms qcmaquis 2 particle dm into pyscf compatible format.

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
                    # 1rdm, 2rdm
                    _, maquis_rdm2 = self.dmrg.get_reduced_density_matrices()
                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        # 1rdm, 2rdm
                        _, maquis_rdm2 = self.dmrg.get_reduced_density_matrices()
            else:
                maquis_rdm2 = ([[0, 0, 0, 0]], [0])

        except RuntimeError:
            maquis_rdm2 = ([[0, 0, 0, 0]], [0])

        # convert 2 rdm from qcmaquis to pyscf format
        rdm2 = np.zeros((norb,) * 4)
        for i, vec in enumerate(maquis_rdm2[0]):
            rdm2[vec[0], vec[1], vec[2], vec[3]] = maquis_rdm2[1][i]
            rdm2[vec[2], vec[3], vec[0], vec[1]] = maquis_rdm2[1][i]
            rdm2[vec[1], vec[0], vec[3], vec[2]] = maquis_rdm2[1][i]
            rdm2[vec[3], vec[2], vec[1], vec[0]] = maquis_rdm2[1][i]
        rdm2 = rdm2.transpose(0, 3, 1, 2)
        return rdm2

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
        nelec, spin2 = self._check_spin(nelec)

        # Always run DMRG if no file path is set
        if not self.file_path:
            # make sure no chkpfile and result file is set
            self.parameters.erase("chkpfile", verbose=False)
            self.parameters.erase("resultfile", verbose=False)
            if self.verbose > 4:
                self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)
            else:
                with pyscf.lib.capture_stdout() as stdout:
                    self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)

        else:
            self._check_file_path()
            cases = self.check_checkpoint_and_results_file()

            #
            if self._dmrgscf:
                self.dmrg.replace_parameters(self.parameters)
                if self.verbose > 4:
                    self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)
                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)

            # load from checkpoint (but only DMRGCI)
            elif cases == "DMRGCI":
                self.parameters.erase("integrals")
                self.dmrg.replace_parameters(self.parameters)
                print(f"""No DMRG calculation required. MPS is loaded from: {self.parameters.get_checkpoint_path()}""")
                if self.verbose > 4:
                    self.dmrg.init_dmrg(self.parameters.get_checkpoint_path(), norb, nelec, spin2)

                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        self.dmrg.init_dmrg(self.parameters.get_checkpoint_path(), norb, nelec, spin2)

            # cannot load DMRGSCF from checkpoint yet, so we just run dmrgscf from scratch
            elif cases == "DMRGSCF":
                self._dmrgscf = True
                self._check_file_path()
                cases = self.check_checkpoint_and_results_file()
                self.dmrg.replace_parameters(self.parameters)
                if self.verbose > 4:
                    self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)
                else:
                    with pyscf.lib.capture_stdout() as stdout:
                        self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)
            else:
                raise RuntimeError("How did we get here")

        energy = self._get_energy()
        fakewfn_by_rdm2 = self._get_rdm2(norb)

        return energy, fakewfn_by_rdm2

    def make_rdm12(self, fakewfn_by_rdm2: np.ndarray, ncas: int, nelec: int, **kwargs: Any) -> Tuple[np.ndarray, np.ndarray]:
        """Make 1rdm and 2rdm.

        TODO
        ----
        Check if this function is required for pyscf, but I think so.
        """
        if not isinstance(nelec, (int, np.integer)):
            nelec = sum(nelec)
        rdm2 = fakewfn_by_rdm2
        rdm1 = np.einsum('ijkk->ij', rdm2) / (nelec - 1)
        return rdm1, rdm2

    def make_rdm1(self, fcivec: Any, norb: int, nelec: int, link_index: Any = None, **kwargs: Any) -> np.ndarray:
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
                self.parameters.set_checkpoint_path(self.file_path + self.parameters.get_checkpoint_path())

            if not self.parameters.get_result_path().startswith(self.file_path):
                self.parameters.set_result_path(self.file_path + self.parameters.get_result_path())

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
            if self.parameters.get_checkpoint_path().endswith(self._dmrgscf_checkpoint_name):
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
