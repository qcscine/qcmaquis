"""Hihi."""
import sys
from functools import reduce

import numpy
import numpy as np
import pyscf
from pyscf import ao2mo

from dmrg import MaquisDmrg


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
    print(h1eff.shape)

    return energy_core, h1eff, eri


class QcMaquis:
    """Huhu."""

    def __init__(self, mol, **kwargs):
        """Huhu."""
        self.mol = mol
        self.verbose = mol.verbose
        self.orbsym = []
        # Maquis related stuff
        self.method = "conventional"
        self.feast_window = None
        self.feast_states = None
        self.n_states = None
        self.nroots = 1
        self.fiedler = False

        self.log = None

        # self.stdout = mol.stdout
        # Do I need the following options?
        # self.wfn_irrep = 0
        # self.spin_2s = 0  # spin = 2*s, 0 means singlet
        # self.dmrg_states = [200, 500, 1000, 1000]
        # self.dmrg_noise = [1, 1, 1, 0]
        # self.dmrg_e_convergence = 1e-16
        # self.dmrg_noise_factor = 0.03
        # self.dmrg_maxiter_noise = 5
        # self.dmrg_maxiter_silent = 100

        self._keys = set(self.__dict__.keys())

    def dump_flags(self, verbose=None):
        """Haha."""
        self.log = pyscf.lib.logger.new_logger(self, verbose)
        self.log.info('******** QcMaquis flags ********')
        self.log.info('method = %s', str(self.method))
        if self.feast_window is not None and self.feast_states is not None:
            self.log.info('feast window = %s', str(self.feast_window))
            self.log.info('feast states = %s', str(self.feast_states))
        elif self.n_states is not None:
            self.log.info('n orthogonal states = %s', str(self.n_states))

    def _get_rdm2(self, norb):
        # in case of feast there is no rdm
        try:
            with pyscf.lib.capture_stdout() as stdout:
                maquis_rdm1, maquis_rdm2 = self.dmrg.get_reduced_density_matrices()
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

    def _set_excited_state_options(self):
        # Consistency with excited states
        if self.n_states is None:
            pass
        elif self.n_states >= 1:
            self.n_states = None

        # FEAST
        if self.method.lower() == "feast":
            if self.feast_states is not None and self.feast_window is not None:
                self.dmrg.set_feast(self.feast_window, self.feast_states)
            else:
                raise ValueError("For a feast calculation a feast_window and feast_states")

        # ORTHO
        if self.method.lower() == "ortho":
            if self.n_states is None:
                raise ValueError("For a feast calculation a feast_window and feast_states")

    def _check_spin(self, nelec):
        # check spin
        if isinstance(nelec, (int, numpy.integer)):
            spin2 = 0
        else:
            spin2 = (nelec[0] - nelec[1])
            nelec = sum(nelec)

        return nelec, spin2

    def _get_energy(self):
        energy = self.dmrg.get_energy()
        # convert complex to float
        if self.method.lower() == "feast":
            for i, ener in enumerate(energy):
                energy[i] = ener.real
        return energy

    def kernel(self, h1e, eri, norb, nelec, ci0=None, ecore=0, **kwargs):
        """Haha."""

        self.orbsym = numpy.zeros(norb, numpy.int32)
        eri = pyscf.ao2mo.restore(1, eri, norb)

        # initialize qcmaquis
        self.dmrg = MaquisDmrg()
        # self.dmrg._parameters.set("symmetry", "su2u1pg")
        # enable 1 and 2 rdm
        self.dmrg.set_orbital_optimization()

        # set eris
        self.dmrg.set_integrals(ecore, h1e, eri, norb)

        self._set_excited_state_options()

        nelec, spin2 = self._check_spin(nelec)

        if self.verbose > 4:
            self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)
        else:
            with pyscf.lib.capture_stdout() as stdout:
                self.dmrg.run(norb, nelec, spin2, n_states=self.n_states, fiedler=self.fiedler)

        energy = self._get_energy()
        fakewfn_by_rdm2 = self._get_rdm2(norb)

        return energy, fakewfn_by_rdm2

    def make_rdm12(self, fakewfn_by_rdm2, ncas, nelec, **kwargs):
        """Haha."""
        if not isinstance(nelec, (int, numpy.integer)):
            nelec = sum(nelec)
        rdm2 = fakewfn_by_rdm2
        rdm1 = numpy.einsum('ijkk->ij', rdm2) / (nelec - 1)
        return rdm1, rdm2

    def make_rdm1(self, fcivec, norb, nelec, link_index=None, **kwargs):
        """Haha."""
        return QcMaquis.make_rdm12(self, fcivec, norb, nelec, **kwargs)[0]


def au_to_ev(energy):
    return energy * 27.211324570273


def ev_to_au(energy):
    return energy / 27.211324570273


if __name__ == '__main__':
    from pyscf import gto, lib, mcscf, scf

    # allow only 4 threads
    lib.num_threads(2)

    # b = 1.4
    mol = gto.Mole()

    from pyscf import mrpt
    b = 1.4
    mol = gto.Mole()
    mol.build(
        verbose=1,
        # output='out-casscf',
        atom=[['H', (0., 0., i)] for i in range(8)],
        basis={'H': '6-31g'},
        symmetry=False,
        spin=2,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASCI(m, 4, 4)
    mc.fcisolver = QcMaquis(mol)
    emc_0 = mc.casci()[0]

    mc = mcscf.CASSCF(m, 4, 4)
    mc.max_cycle_macro = 20
    mc.fcisolver = QcMaquis(mol)
    emc_1 = mc.mc2step()[0]

    # nevpt2_ener = mrpt.NEVPT2(mc).kernel()

    print("**********************************")
    b = 1.4
    mol = gto.Mole()
    mol.build(
        verbose=1,
        # output='out-casscf',
        atom=[['H', (0., 0., i)] for i in range(8)],
        basis={'H': '6-31g'},
        symmetry=False,
        spin=2,
    )
    m = scf.RHF(mol)
    m.scf()

    mc = mcscf.CASCI(m, 4, 4)
    emc_0ref = mc.casci()[0]

    mc = mcscf.CASSCF(m, 4, 4)
    mc.max_cycle_macro = 20
    emc_1ref = mc.mc2step()[0]

    # nevpt2_ener_2 = mrpt.NEVPT2(mc).kernel()

    print('Maquis-CI  = %.15g CASCI  = %.15g' % (emc_0, emc_0ref))
    print('Diff = %.15g' % (emc_0 - emc_0ref))
    print('Maquis-SCF = %.15g CASSCF = %.15g' % (emc_1, emc_1ref))
    print('Diff = %.15g' % (emc_1 - emc_1ref))
    # print('NEVPT2 maquis = %.15g NEVPT2 = %.15g' % (nevpt2_ener, nevpt2_ener_2))
    # print('Diff = %.15g' % (nevpt2_ener - nevpt2_ener_2))
