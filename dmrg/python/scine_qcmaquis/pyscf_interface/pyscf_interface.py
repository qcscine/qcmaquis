import os
import sys
from functools import reduce

import numpy
import numpy as np
import pyscf
from pyscf import ao2mo

from .. import MaquisDmrg


# TODO this is not needed here, but maybe important somewhere else -> move it to better place
def make_cas_hamil(mf, norb, nelecs):
    mo_coeff = mf.mo_coeff
    nocc_tot = int(sum(mf.mo_occ) / 2)
    ncore = int(nocc_tot - nelecs / 2)
    # nocc = nocc_tot
    eris = mf._eri
    ncas = norb

    # 2e
    # eri_cas = casci.get_h2eff(mocas)
    if mo_coeff.shape[1] != ncas:
        mo_coeff_tmp = mo_coeff[:, ncore:ncore + ncas]
    eri = ao2mo.full(eris, mo_coeff_tmp)
    eri = ao2mo.restore(1, eri, norb)

    # 1e
    # h1eff, energy_core = casci.get_h1eff(mocas)
    mo_core = mo_coeff[:, :ncore]
    mo_cas = mo_coeff[:, ncore:ncore + ncas]
    hcore = mf.get_hcore()
    energy_core = mf.energy_nuc()

    core_dm = np.dot(mo_core, mo_core.conj().T) * 2
    corevhf = mf.get_veff(None, core_dm)
    energy_core += np.einsum('ij,ji', core_dm, hcore).real
    energy_core += np.einsum('ij,ji', core_dm, corevhf).real * .5
    h1eff = reduce(np.dot, (mo_cas.conj().T, hcore + corevhf, mo_cas))

    return energy_core, h1eff, eri


class QcMaquis:
    """QCMaquis interface for pyscf.

    This interface allows qcmaquis to act as a fcisolver in pyscf.mcscf methods.

    TODO
    ----
    For MRPT higher order rdms are required, and the corresponding functions have
    to be implemented here, in order to work automatically by injection of this
    interface into pyscf.
    """

    def __init__(self, mol, **kwargs):
        """Initialize interface."""
        self.mol = mol
        """Pyscf molecule"""
        self.verbose = mol.verbose
        """Output verbosity"""
        self.orbsym = []
        """Orbital symmetries"""
        # DMRG related stuff
        self.method = "conventional"
        """DMRG method, e.g. conventional, ..."""
        self.n_states = None
        """Number of states."""
        self.fiedler = False
        """Enable fiedler ordering"""
        # self.entropies = None

        self.log = None
        """Pyscf logger."""

        # QCMaquis related stuff
        self.dmrg = None
        """QCMaquis python interface, will be set automatically."""
        self.file_path = os.getcwd()
        """Path for qcmaquis dumps."""
        self.measure_entropies = False
        """Enable measurements for s1, s2 and mut inf."""
        self.checkpoint_name = "qcmaquis_checkpoint.h5"
        """Name of qcmaquis checkpoint."""
        self.results_name = "qcmaquis_result_file.h5"
        """Name of qcmaquis results file."""

        # DMRG parameters
        self.energy_threshold = 1e-6
        """Energy threshold."""
        self.nsweeps = 100
        """Max number of sweeps."""
        self.bond_dim = 250
        """Bond dimension."""
        self.initial_truncation_thresh = 1e-16  # 1e-6

        self.final_truncation_thresh = 1e-16  # 1e-10

        # Set other parameters directly through qcmaquis parameters

        # self.feast_window = None
        # self.feast_states = None
        # self.nroots = 1
        # self.stdout = mol.stdout

        self._keys = set(self.__dict__.keys())

    def dump_flags(self, verbose=None):
        """Log settings for pyscf logger."""
        self.log = pyscf.lib.logger.new_logger(self, verbose)
        self.log.info('************** QcMaquis flags **************')
        self.log.info('method           = %s', str(self.method))
        self.log.info('energy_threshold = %s', str(self.energy_threshold))
        self.log.info('n sweeps         = %s', str(self.nsweeps))
        self.log.info('bond dimension   = %s', str(self.bond_dim))
        self.log.info('fiedler ordering = %s', str(self.fiedler))
        self.log.info('entropies        = %s', str(self.measure_entropies))
        self.log.info('checkpoint name  = %s', str(self.checkpoint_name))
        self.log.info('checkpoint path  = %s', str(self.file_path))
        self.log.info('results name     = %s', str(self.results_name))

    def _get_rdm2(self, norb):
        """Getter for 2 rdm.

        Transforms qcmaquis 2 particle dm into pyscf compatible format.
        """
        # in case of feast there is no rdm
        try:
            if self.verbose > 4:
                maquis_rdm1, maquis_rdm2 = self.dmrg.get_reduced_density_matrices()
            else:
                with pyscf.lib.capture_stdout() as stdout:
                    maquis_rdm1, maquis_rdm2 = self.dmrg.get_reduced_density_matrices()
        # in case measurement failed
        # TODO: check if this is correct
        except RuntimeError:
            maquis_rdm2 = ([[0, 0, 0, 0]], [0])
        # sometimes there is another issue than runtime
        except:
            maquis_rdm2 = ([[0, 0, 0, 0]], [0])

        # convert 2 rdm from qcmaquis to pyscf format
        rdm2 = numpy.zeros((norb,) * 4)
        for i, vec in enumerate(maquis_rdm2[0]):
            rdm2[vec[0], vec[1], vec[2], vec[3]] = maquis_rdm2[1][i]
            rdm2[vec[2], vec[3], vec[0], vec[1]] = maquis_rdm2[1][i]
            rdm2[vec[1], vec[0], vec[3], vec[2]] = maquis_rdm2[1][i]
            rdm2[vec[3], vec[2], vec[1], vec[0]] = maquis_rdm2[1][i]
        rdm2 = rdm2.transpose(0, 3, 1, 2)
        return rdm2

    # TODO enable excited states
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

    def _check_spin(self, nelec):
        """Check spin based on numpy electron definition."""
        # check spin
        if isinstance(nelec, (int, numpy.integer)):
            spin2 = 0
        else:
            spin2 = (nelec[0] - nelec[1])
            nelec = sum(nelec)

        return nelec, spin2

    def _get_energy(self):
        """Get energy from interface."""
        energy = self.dmrg.get_energy()
        # convert complex to float
        # if self.method.lower() == "feast":
        #     for i, ener in enumerate(energy):
        #         energy[i] = ener.real
        return energy

    def _set_parameters(self):
        self.dmrg.set_bond_dimension(self.bond_dim)
        self.dmrg.set_parameter("conv_thresh", self.energy_threshold)
        self.dmrg.set_parameter("nsweeps", self.nsweeps)
        self.dmrg.set_parameter("truncation_initial", self.initial_truncation_thresh)
        self.dmrg.set_parameter("truncation_final", self.final_truncation_thresh)

    def kernel(self, h1e, eri, norb, nelec, ci0=None, ecore=0, **kwargs):
        """Kernel function for pyscf.

        Create the interface (at the moment integrals are required in order to
        initialize the interface) and run qcmaquis calculation.
        """

        self.orbsym = numpy.zeros(norb, numpy.int32)
        eri = pyscf.ao2mo.restore(1, eri, norb)

        # initialize qcmaquis
        self.dmrg = MaquisDmrg()
        if self.file_path:
            self.dmrg._parameters.set_checkpoint_path(self.file_path + "/" + self.checkpoint_name)
            self.dmrg._parameters.set_result_path(self.file_path + "/" + self.results_name)

        # enable 1 and 2 rdm
        self.dmrg.set_orbital_optimization()
        self._set_parameters()

        # enable chementropy measurement
        if self.measure_entropies is True:
            self.dmrg.set_entropies()

        # set eris
        self.dmrg.set_integrals(ecore, h1e, eri, norb)

        self._set_excited_state_options()

        nelec, spin2 = self._check_spin(nelec)

        if not os.path.exists(self.results_name):
            if self.verbose > 4:
                self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)
            else:
                with pyscf.lib.capture_stdout() as stdout:
                    self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)
        else:
            print(f"""No DMRG calculation required.
MPS is loaded from: {self.checkpoint_name} in {os.getcwd()}""")
            if self.verbose > 4:
                self.dmrg.init_dmrg(self.checkpoint_name, norb, nelec, spin2)
            else:
                with pyscf.lib.capture_stdout() as stdout:
                    self.dmrg.init_dmrg(self.checkpoint_name, norb, nelec, spin2)

        energy = self._get_energy()
        fakewfn_by_rdm2 = self._get_rdm2(norb)

        return energy, fakewfn_by_rdm2

    def make_rdm12(self, fakewfn_by_rdm2, ncas, nelec, **kwargs):
        """Make 1rdm and 2rdm.

        TODO
        ----
        Check if this function is required for pyscf, but I think so.
        """
        if not isinstance(nelec, (int, numpy.integer)):
            nelec = sum(nelec)
        rdm2 = fakewfn_by_rdm2
        rdm1 = numpy.einsum('ijkk->ij', rdm2) / (nelec - 1)
        return rdm1, rdm2

    def make_rdm1(self, fcivec, norb, nelec, link_index=None, **kwargs):
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


# def au_to_ev(energy):
#     return energy * 27.211324570273
#
#
# def ev_to_au(energy):
#     return energy / 27.211324570273
#
#
# if __name__ == '__main__':
#     from pyscf import gto, lib, mcscf, scf
#
#     # allow only 4 threads
#     lib.num_threads(2)
#
#     # b = 1.4
#     mol = gto.Mole()
#
#     from pyscf import mrpt
#     b = 1.4
#     mol = gto.Mole()
#     mol.build(
#         verbose=1,
#         # output='out-casscf',
#         atom=[['H', (0., 0., i)] for i in range(8)],
#         basis={'H': '6-31g'},
#         symmetry=False,
#         spin=2,
#     )
#     m = scf.RHF(mol)
#     m.scf()
#
#     mc = mcscf.CASCI(m, 4, 4)
#     mc.fcisolver = QcMaquis(mol)
#     emc_0 = mc.casci()[0]
#
#     mc = mcscf.CASSCF(m, 4, 4)
#     mc.max_cycle_macro = 20
#     mc.fcisolver = QcMaquis(mol)
#     emc_1 = mc.mc2step()[0]
#
#     # nevpt2_ener = mrpt.NEVPT2(mc).kernel()
#
#     print("**********************************")
#     b = 1.4
#     mol = gto.Mole()
#     mol.build(
#         verbose=1,
#         # output='out-casscf',
#         atom=[['H', (0., 0., i)] for i in range(8)],
#         basis={'H': '6-31g'},
#         symmetry=False,
#         spin=2,
#     )
#     m = scf.RHF(mol)
#     m.scf()
#
#     mc = mcscf.CASCI(m, 4, 4)
#     emc_0ref = mc.casci()[0]
#
#     mc = mcscf.CASSCF(m, 4, 4)
#     mc.max_cycle_macro = 20
#     emc_1ref = mc.mc2step()[0]
#
#     # nevpt2_ener_2 = mrpt.NEVPT2(mc).kernel()
#
#     print('Maquis-CI  = %.15g CASCI  = %.15g' % (emc_0, emc_0ref))
#     print('Diff = %.15g' % (emc_0 - emc_0ref))
#     print('Maquis-SCF = %.15g CASSCF = %.15g' % (emc_1, emc_1ref))
#     print('Diff = %.15g' % (emc_1 - emc_1ref))
#     # print('NEVPT2 maquis = %.15g NEVPT2 = %.15g' % (nevpt2_ener, nevpt2_ener_2))
#     # print('Diff = %.15g' % (nevpt2_ener - nevpt2_ener_2))
