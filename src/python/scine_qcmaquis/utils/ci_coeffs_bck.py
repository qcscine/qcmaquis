from typing import List, Tuple

# # TODO: this is not needed here, but maybe important somewhere else -> move it to better place
# def make_cas_hamil(mf: Any, norb: int, nelecs: int) -> Tuple[float, np.ndarray, np.ndarray]:
#     """Create the CAS Hamiltonian for n orbitals and n electrons from a HF wave function.
#
#     Note
#     ----
#     This is not tested or used.
#
#     Parameters
#     ----------
#     mf : Any
#         Pyscf wave function
#     norb : int
#         Number of orbitals
#     nelecs : int
#         Number of electrons
#
#     Returns
#     -------
#     Tuple(float, np.ndarray, np.ndarray)
#         core energy, effective one body operator, eris
#     """
#     mo_coeff = mf.mo_coeff
#     nocc_tot = int(sum(mf.mo_occ) / 2)
#     ncore = int(nocc_tot - nelecs / 2)
#     # pylint: disable = W0212
#     eris = mf._eri
#     # pylint: enable = W0212
#     ncas = norb
#
#     # 2e
#     if mo_coeff.shape[1] != ncas:
#         mo_coeff_tmp = mo_coeff[:, ncore:ncore + ncas]
#     eri = ao2mo.full(eris, mo_coeff_tmp)
#     eri = ao2mo.restore(1, eri, norb)
#
#     # 1e
#     mo_core = mo_coeff[:, :ncore]
#     mo_cas = mo_coeff[:, ncore:ncore + ncas]
#     hcore = mf.get_hcore()
#     energy_core = mf.energy_nuc()
#
#     core_dm = np.dot(mo_core, mo_core.conj().T) * 2
#     corevhf = mf.get_veff(None, core_dm)
#     energy_core += np.einsum('ij,ji', core_dm, hcore).real
#     energy_core += np.einsum('ij,ji', core_dm, corevhf).real * .5
#     h1eff = reduce(np.dot, (mo_cas.conj().T, hcore + corevhf, mo_cas))
#
#     return energy_core, h1eff, eri

# def make_hf_occ(norb: int, nocc: int) -> List[int]:
#     """
#     Create restricted hartree fock determinant.
#
#     Create determinant from number of occupied orbitals and total number of
#     orbitals.
#
#     norb: int
#         total number of orbitals
#     nocc: int
#         number of occupied orbitals
#     """
#     hf_list = []
#     for i in range(norb):
#         if i < nocc:
#             hf_list.append(1)
#         else:
#             hf_list.append(0)
#     return hf_list
#
#
# def exciste_elec(from_orb: int, to_orb: int, hf_occ: List[int]) -> List[int]:
#     """
#     Promote one electron from the occupied to the virtual space.
#
#     e.g. create a new occupied orbital int the virtual space and remove it from
#     the occupied.
#
#     from_orb: int
#         index of the occupied orbital to remove
#     to_orb: int
#         index of the virtual orbital to occupy
#     hf_occ: list[int]
#         occupation number vector as base for the excitation
#     """
#     if hf_occ[from_orb] == 0:
#         raise ValueError(f"orb {from_orb} is not occupied")
#     if hf_occ[to_orb] == 1:
#         raise ValueError(f"orb {to_orb} is occupied")
#     new_list = hf_occ.copy()
#
#     new_list[from_orb] = 0
#     new_list[to_orb] = 1
#     return new_list
#
#
# def permute_alpha_beta_for_sign(qcmaquis_string: str) -> str:
#     """
#     Permute alpha and beta orbitals for correct sign.
#
#     The orbital order of PySCF is different from the one in qcmaquis.
#     In order to get to the correct order, orbitals have to be
#     permuted and hence the sign may change.
#
#     qcmaquis_string: str
#         4 double, 3 alpha, 2 beta and 1 unoccupied
#         e.g. 4,4,4,3,2,1,1
#     """
#     sign = 1
#
#     alpha_string = ""
#     beta_string = ""
#     for orbital in qcmaquis_string:
#         if orbital == "4":
#             alpha_string += "1"
#             beta_string += "1"
#         elif orbital == "3":
#             alpha_string += "1"
#             beta_string += "0"
#         elif orbital == "2":
#             alpha_string += "0"
#             beta_string += "1"
#         elif orbital == "1":
#             alpha_string += "0"
#             beta_string += "0"
#
#     for index, orbital in enumerate(beta_string):
#         if orbital == "1":
#             for blub in range(index + 1, len(alpha_string)):
#                 if alpha_string[blub] == "1":
#                     sign *= -1
#     return sign
#
#
# def make_qcmaquis_string(alpha, beta):
#     """
#     Make the qcmaquis string.
#
#     alpha and beta are occupation number vectors for alpha and beta orbitals.
#     """
#     alpha = [i * 2 for i in alpha]
#     qcmaquis_list = [sum(x) + 1 for x in zip(alpha, beta)]
#     qcmaquis_str = ",".join(str(i) for i in qcmaquis_list)
#     sign = permute_alpha_beta_for_sign(qcmaquis_str)
#     return qcmaquis_str, sign
#
#
# def make_ref(nocc: int, norb: int,) -> Tuple[str, int]:
#     """
#     Make the reference string for qcmaquis.
#
#     nocc: int
#         number of spatial occupied orbitals
#     norb: int
#         number of spatial orbitals
#     c0: float
#         ci coefficient
#     """
#     hf_list = make_hf_occ(norb, nocc)
#     qcmaquis_str, sign = make_qcmaquis_string(hf_list, hf_list)
#     return (qcmaquis_str, sign)
#
#
# def make_singles_bb(nocc, norb, i, a):
#     """
#     Make the beta single excitation string for qcmaquis.
#
#     nocc:
#     norb:
#     c0:
#     """
#     hf_list = make_hf_occ(norb, nocc)
#     single_aa = exciste_elec(i, nocc + a, hf_list)
#     qcmaquis_str, sign = make_qcmaquis_string(single_aa, hf_list)
#     return (qcmaquis_str, sign)
#
#
# def make_singles_aa(nocc, norb, i, a):
#     """
#     Make the alpha single excitation string for qcmaquis.
#
#     nocc:
#     norb:
#     c0:
#     """
#     hf_list = make_hf_occ(norb, nocc)
#     single_aa = exciste_elec(i, nocc + a, hf_list)
#     qcmaquis_str, sign = make_qcmaquis_string(hf_list, single_aa)
#     return (qcmaquis_str, sign)
#
#
# def make_doubles_bb(nocc, norb, i, j, a, b):
#     """
#     Make the beta double excitation string for qcmaquis.
#
#     nocc:
#     norb:
#     c0:
#     """
#     hf_list = make_hf_occ(norb, nocc)
#     single_aa = exciste_elec(i, nocc + a, hf_list)
#     doubles_aa = exciste_elec(j, nocc + b, single_aa)
#     qcmaquis_str, sign = make_qcmaquis_string(hf_list, doubles_aa)
#     return (qcmaquis_str, sign)
#
#
# def make_doubles_aa(nocc, norb, i, j, a, b):
#     """
#     Make the alpha double excitation string for qcmaquis.
#
#     nocc:
#     norb:
#     c0:
#     """
#     hf_list = make_hf_occ(norb, nocc)
#     single_aa = exciste_elec(i, nocc + a, hf_list)
#     doubles_aa = exciste_elec(j, nocc + b, single_aa)
#     qcmaquis_str, sign = make_qcmaquis_string(doubles_aa, hf_list)
#     return (qcmaquis_str, sign)
#
#
# def make_doubles_ab(nocc, norb, i, j, a, b):
#     """
#     Make the alpha beta double excitation string for qcmaquis.
#
#     nocc:
#     norb:
#     c0:
#     """
#     hf_list = make_hf_occ(norb, nocc)
#     single_aa = exciste_elec(i, nocc + a, hf_list)
#     single_bb = exciste_elec(j, nocc + b, hf_list)
#     qcmaquis_str, sign = make_qcmaquis_string(single_aa, single_bb)
#     return (qcmaquis_str, sign)
#
#
# def make_doubles_ba(nocc, norb, i, j, a, b):
#     """
#     Make the beta alpha double excitation string for qcmaquis.
#
#     nocc:
#     norb:
#     c0:
#     """
#     hf_list = make_hf_occ(norb, nocc)
#     single_aa = exciste_elec(i, nocc + a, hf_list)
#     single_bb = exciste_elec(j, nocc + b, hf_list)
#     qcmaquis_str, sign = make_qcmaquis_string(single_bb, single_aa)
#     return (qcmaquis_str, sign)
#
#
# if __name__ == "__main__":
#     print(make_ref(4, 10))
#     print(make_singles_aa(4, 10, 1, 1))
#     print(make_singles_bb(4, 10, 1, 1))
#     print(make_doubles_aa(4, 10, 1, 2, 1, 2))
#     print(make_doubles_bb(4, 10, 1, 2, 1, 2))
#     print(make_doubles_ab(4, 10, 1, 2, 1, 2))
#     print(make_doubles_ba(4, 10, 1, 2, 1, 2))
#     pass
