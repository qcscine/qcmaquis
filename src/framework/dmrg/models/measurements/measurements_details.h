/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */
#ifndef MEASUREMENTS_DETAILS_H
#define MEASUREMENTS_DETAILS_H

#include "dmrg/models/model.h"
namespace measurements_details {

template <class symm, class = void>
class checkpg {
 public:
  template <class matrix>
  bool operator()(
      term_descriptor<typename matrix::value_type> const& term,
      std::shared_ptr<TagHandler<matrix, symm>> tag_handler, Lattice const& lat
  ) {
    return true;
  }
};

template <class symm>
class checkpg<symm, symm_traits::enable_if_pg_t<symm>> {
 public:
  using charge = typename symm::charge;
  using subcharge = typename symm::subcharge;

  template <class matrix>
  bool operator()(
      term_descriptor<typename matrix::value_type> const& term,
      std::shared_ptr<TagHandler<matrix, symm>> tag_handler, Lattice const& lat
  ) {
    using op_t = typename TagHandler<matrix, symm>::op_t;

    charge acc = symm::IdentityCharge;
    for (std::size_t p = 0; p < term.size(); ++p) {
      charge local = symm::IdentityCharge;
      if (tag_handler->is_fermionic(term.operator_tag(p))) {
        // stknecht: this check does not work properly for U1DG. FIXME!
        if (symm_traits::HasU1DG<symm>::value) {
          if (p % 2 == 0) {
            symm::irrep(local) =
                lat.get_prop<subcharge>("type", term.position(p));
          } else {
            symm::irrep(local) =
                symm::adjoin(lat.get_prop<subcharge>("type", term.position(p)));
          }
        } else {
          symm::irrep(local) =
              lat.get_prop<subcharge>("type", term.position(p));
        }
      }

      // maquis::cout << " index " << p << " --> accumulated charge (before) "
      // << acc << " local charge " << local << std::endl;
      acc = symm::fuse(acc, local);
      // maquis::cout << " index " << p << " --> accumulated charge (after ) "
      // << acc << " local charge " << local << std::endl;
    }

    return (acc == symm::IdentityCharge);
  }
};

template <class T>
inline T get_indx_contr(std::vector<T> const& positions) {
  T contr_indx;
  return (
      ((positions[0] + 1 - 1) * (positions[0] + 1) * (positions[0] + 1 + 1) *
       (positions[0] + 1 + 2)) /
          24 +
      ((positions[1] + 1 - 1) * (positions[1] + 1) * (positions[1] + 1 + 1)) /
          6 +
      ((positions[2] + 1 - 1) * (positions[2] + 1)) / 2 + positions[3] + 1
  );
};

template <class T>
inline bool compare_norm(std::vector<T> const& pos) {
  std::vector<T> positions_lhs{pos[0], pos[1], pos[2], pos[3]};
  std::vector<T> positions_rhs{pos[4], pos[5], pos[6], pos[7]};

  // reverse sorting to ensure maximum norm
  std::sort(positions_lhs.begin(), positions_lhs.end(), std::greater<T>());
  std::sort(positions_rhs.begin(), positions_rhs.end(), std::greater<T>());

  T norm_lhs = get_indx_contr(positions_lhs);
  T norm_rhs = get_indx_contr(positions_rhs);

  // maquis::cout << "lhs norm "  << norm_lhs << " <--> rhs norm " << norm_rhs
  // << std::endl;

  return (norm_rhs > norm_lhs);
}

// Function to handle n-RDM index permutations
// it is written in a generic way to avoid copy-paste
// for now it is used to
// a) get the total number of permutations for a given n-RDM or its slice.
// b) construct the index tuples to iterate over
// Params:
// L: number of orbitals
// bra_neq_ket -- true for transition RDMs, false otherwise Unsupported
// (ignored) for 4-RDM. positions_first: optional -- fixed first 4 indices to
// obtain the slice of 4-RDM (or 2-3 for 3-RDM, unused for lower RDM) fun: a
// function (wrapped in a class that has return_type and operator() defined)
// that returns F::return_type that gets executed in the middle of a loop
// (e.g. increments a counter or adds an index vector)
// examples for fun see below
template <class F, int N>
struct iterate_rdm_indices {
  using pos_t = Lattice::pos_t;
  typename F::return_type operator()(
      F fun, pos_t L, bool bra_neq_ket = false,
      const std::vector<pos_t>& positions_first = std::vector<pos_t>()
  ) {
    throw std::runtime_error(
        "iterate_rdm_indices not implemented for this number of indices"
    );
  }
};

// specialization for 4-RDMs
template <class F>
struct iterate_rdm_indices<F, 4> {
  using pos_t = Lattice::pos_t;
  typename F::return_type operator()(
      F fun, pos_t L, bool bra_neq_ket = false,
      const std::vector<pos_t>& positions_first = std::vector<pos_t>()
  ) {
    pos_t p4_start = 0;
    pos_t p3_start = 0;
    pos_t p1_start = L - 1;
    pos_t p2_start = p1_start;
    pos_t p4_end = L - 1;
    pos_t p3_end = L - 1;
    pos_t p1_end = 0;
    pos_t p2_end = 0;
    pos_t p_max = L;
    if (positions_first.size() == 4) {
      p4_start = positions_first[0];
      p3_start = positions_first[1];
      p1_start = positions_first[2];
      p2_start = positions_first[3];
      p4_end = positions_first[0] + 1;
      p3_end = positions_first[1] + 1;
      p1_end = positions_first[2];
      p2_end = positions_first[3];
    }
    for (pos_t p4 = p4_start; p4 < p4_end; ++p4) {
      for (pos_t p3 = p3_start; p3 < p3_end; ++p3) {
        for (pos_t p1 = p1_start; p1 >= p1_end; --p1) {
          if (p4 > p3 || p4 > p1 || p3 > p1) {
            continue;
          }

          if (positions_first.empty()) {
            p2_start = p1;
            p2_end = 0;
          }

          for (pos_t p2 = p2_start; p2 >= p2_end; --p2) {
            if (p3 > p2) {
              continue;
            }

            // N-rdm are zero for > 2 identical creation ops
            if ((p1 == p2 && p1 == p3) || (p1 == p2 && p1 == p4) ||
                (p1 == p3 && p1 == p4) || (p2 == p3 && p2 == p4)) {
              continue;
            }

            bool double_equal = (p1 == p2 && p3 == p4);            // case 1
            bool ij_equal = (p1 == p2 && p2 != p3 && p3 != p4);    // case 2
            bool jk_equal = (p1 != p2 && p2 == p3 && p3 != p4);    // case 3
            bool kl_equal = (p1 != p2 && p2 != p3 && p3 == p4);    // case 4
            bool none_equal = (p1 != p2 && p2 != p3 && p3 != p4);  // case 5

            for (pos_t p5 = p1; p5 >= 0; --p5) {
              for (pos_t p6 = p1; p6 >= 0; --p6) {
                for (pos_t p7 = 0; p7 < p_max; ++p7) {
                  // set restrictions on index p6
                  if ((double_equal || ij_equal) && p6 > p5) {
                    continue;
                  }

                  // set restrictions on index p7
                  if (double_equal) {
                    if (p7 > p5) {
                      continue;
                    }
                  } else {
                    if (p7 > p1) {
                      continue;
                    }
                  }

                  if (p5 == p6 && p5 == p7) {
                    continue;
                  }

                  pos_t p8_end = 0;
                  if (double_equal) {
                    // set restrictions on index p8
                    p8_end = p5 + 1;
                  } else if (kl_equal) {
                    p8_end = p7 + 1;
                  } else {
                    p8_end = p1 + 1;
                  }

                  for (pos_t p8 = 0; p8 < p8_end; ++p8) {
                    // eighth index must be different if p5 == p6 or p5 == p7 or
                    // p6 == p7
                    if ((p5 == p6 && p8 == p5) || (p5 == p7 && p8 == p5) ||
                        (p6 == p7 && p8 == p6)) {
                      continue;
                    }

                    // case 1
                    if (double_equal && p8 > p7 &&
                        (p8 < p6 || p8 == p5 || p8 == p6 || p5 == p6 ||
                         p5 == p7 || p6 == p7)) {
                      continue;
                    }
                    if (double_equal && p8 > p7 && p7 > p6 &&
                        (p8 < p6 || p8 == p5 || p8 == p6 || p5 == p6 ||
                         p5 == p7 || p6 == p7)) {
                      continue;
                    }
                    if (double_equal && p8 > p7 && p7 > p6 &&
                        (p8 > p6 || p8 == p5 || p8 == p6 || p5 == p6 ||
                         p5 == p7 || p6 == p7)) {
                      continue;
                    }
                    if (double_equal && p8 < p7 && p7 > p6 &&
                        (p8 == p5 || p8 == p6 || p7 == p5 || p7 == p6)) {
                      continue;
                    }
                    if (double_equal && p8 < p7 && p7 > p6 && p8 < p6) {
                      continue;
                    }
                    // case 2/3/4: 2x2 equal indices
                    if ((ij_equal || jk_equal || kl_equal) && p5 == p6 &&
                        p7 == p8 && p7 > p6) {
                      continue;
                    }

                    // case 2/4: 2 equal indices
                    if ((ij_equal || kl_equal) &&
                        (p5 == p7 || p5 == p8 || p6 == p7 || p6 == p8)) {
                      continue;
                    }

                    // case 3
                    if (jk_equal) {
                      // 2 equal indices
                      if (p5 == p7 || p6 == p7 || p6 == p8) {
                        continue;
                      }
                      if (p5 == p6 && p7 != p8 && p6 > p7) {
                        continue;
                      }
                      if (p5 == p8 && p7 > p6) {
                        continue;
                      }
                      // none equal
                      if (std::min(p5, p6) != std::min(p7, p8) && p7 != p8 &&
                          p7 > p6) {
                        continue;
                      }
                    }

                    // case 5
                    if (none_equal) {
                      if ((p5 == p6 && p7 == p8 && p5 < p7) ||
                          (p5 == p7 && p6 == p8 && p5 < p6) ||
                          (p5 == p8 && p6 == p7 && p5 < p6)) {
                        continue;
                      }
                    }

                    // defines position vector for spin-free 4-RDM element
                    std::vector<pos_t> positions{p1, p2, p3, p4,
                                                 p5, p6, p7, p8};

                    // check norm of lhs and rhs - skip if norm of rhs > lhs
                    if (compare_norm(positions)) {
                      continue;
                    }

                    // execute functor
                    fun(positions);
                  }
                }
              }
            }
          }
        }
      }
    }
    return fun.get();
  }
};

// Same for 3-RDM indices
// positions[p, q, r, s, t, u] = p_1+ q_2+ r_3+ s_1 t_2 u_3
template <class F>
struct iterate_rdm_indices<F, 3> {
  using pos_t = Lattice::pos_t;
  typename F::return_type operator()(
      F fun, pos_t L, bool bra_neq_ket = false,
      const std::vector<pos_t>& positions_first = std::vector<pos_t>()
  ) {
    // simple version for full transition RDM
    // BUG: Does not work for some reason
    if (positions_first.empty() && bra_neq_ket) {
      for (pos_t p1 = 0; p1 < L; ++p1) {
        for (pos_t p2 = p1; p2 < L; ++p2) {
          for (pos_t p3 = p2; p3 < L; ++p3) {
            bool three_identical_creation_ops = (p1 == p2 && p1 == p3);
            if (three_identical_creation_ops) {
              continue;
            }
            for (pos_t p4 = 0; p4 < L; ++p4) {
              for (pos_t p5 = p4; p5 < L; ++p5) {
                for (pos_t p6 = p5; p6 < L; ++p6) {
                  bool three_identical_annhilation_ops = (p4 == p5 && p4 == p6);
                  if (three_identical_annhilation_ops) {
                    continue;
                  }
                  std::vector<pos_t> positions{p1, p2, p3, p4, p5, p6};
                  fun(positions);
                }
              }
            }
          }
        }
      }
      return fun.get();
    }

    // generic version with slicing
    // NOTE: not perfect, contains redundant elements
    pos_t p1_start = 0;
    pos_t p2_start = 0;
    pos_t p3_start = 0;
    pos_t p1_end = L;
    pos_t p2_end = L;
    pos_t p3_end = L;
    pos_t p4_end = L;
    pos_t p5_end = L;
    pos_t p6_end = L;

    if (positions_first.size() == 2) {
      p1_start = positions_first[0];
      p2_start = positions_first[1];
      p1_end = positions_first[0] + 1;
      p2_end = positions_first[1] + 1;
    }

    // Leon: allow both 2-index and 3-index 3-RDM measurement splitting
    // warning: if two indices are used, the order in the input file is p1,p2;
    // but for three it's p2,p1,p3 !

    if (positions_first.size() == 3) {
      p1_start = positions_first[0];
      p2_start = positions_first[1];
      p3_start = positions_first[2];
      p1_end = positions_first[0] + 1;
      p2_end = positions_first[1] + 1;
      p3_end = positions_first[2] + 1;
    }

    for (pos_t p1 = p1_start; p1 < p1_end; ++p1) {
      for (pos_t p2 = p2_start; p2 < p2_end; ++p2) {
        for (pos_t p3 = p3_start; p3 < p3_end; ++p3) {
          for (pos_t p4 = 0; p4 < p4_end; ++p4) {
            for (pos_t p5 = 0; p5 < p5_end; ++p5) {
              // index restrictions
              if ((p1 < p2) || (p3 < std::min(p1, p2))) {
                continue;
              }
              // N-rdm are zero for > 2 identical creation ops
              bool three_identical_creation_ops = (p1 == p2 && p1 == p3);
              if (three_identical_creation_ops) {
                continue;
              }
              if ((!bra_neq_ket && p4 < std::min(p1, p2)) ||
                  (!bra_neq_ket && p5 < std::min(p1, p2))) {
                continue;
              }

              for (pos_t p6 = std::min(p4, p5); p6 < p6_end; ++p6) {
                // N-rdms are zero for > 2 identical annhilation ops
                bool three_identical_annhilation_ops = p4 == p5 && p4 == p6;
                if (three_identical_annhilation_ops) {
                  continue;
                }

                // do with the indices what's required to do -- define a
                // vector with all positions and pass it on to the functor
                std::vector<pos_t> positions{p1, p2, p3, p4, p5, p6};
                fun(positions);
              }
            }
          }
        }
      }
    }
    return fun.get();
  }
};

// 2-RDMs
// positions[p,q,r,s] = p_1+ q_2+ r_2 s_1
template <class F>
struct iterate_rdm_indices<F, 2> {
  using pos_t = Lattice::pos_t;
  typename F::return_type operator()(
      F fun, pos_t L, bool bra_neq_ket = false,
      const std::vector<pos_t>& positions_first = std::vector<pos_t>()
  ) {
    // Permutation symmetry for bra == ket: pqrs == qpsr == rspq == srqp
    // if bra != ket, pertmutation symmetry is only pqrs == qpsr
    if (!bra_neq_ket) {
      for (pos_t p1 = 0; p1 < L; ++p1) {
        for (pos_t p2 = p1; p2 < L; ++p2) {
          for (pos_t p3 = p1; p3 < L; ++p3) {
            for (pos_t p4 = p1; p4 < L; ++p4) {
              // if two indices are the same the other two must be sorted
              if (((p1 == p2) && (p3 > p4)) || ((p1 == p3) && (p2 > p4)) ||
                  ((p1 == p4) && (p2 > p3))) {
                continue;
              }
              std::vector<pos_t> positions{p1, p2, p3, p4};
              fun(positions);
            }
          }
        }
      }
    } else {
      for (pos_t p1 = 0; p1 < L; ++p1) {
        for (pos_t p2 = 0; p2 < L; ++p2) {
          // if bra != ket, pertmutation symmetry is only pqrs == qpsr
          for (pos_t p3 = 0; p3 < L; ++p3) {
            for (pos_t p4 = p3; p4 < L; ++p4) {
              if ((p3 == p4) && (p1 < p2)) {
                continue;
              }
              std::vector<pos_t> positions{p1, p2, p3, p4};
              fun(positions);
            }
          }
        }
      }
    }
    // original should verify that transition rdm matches reference
    // for (pos_t p1 = 0; p1 < L; ++p1) {
    //   for (pos_t p2 = 0; p2 < L; ++p2) {
    //     // Permutation symmetry for bra == ket: pqrs == qpsr == rspq ==
    //     srqp
    //     // if bra != ket, pertmutation symmetry is only pqrs == qpsr
    //     for (pos_t p3 = (bra_neq_ket) ? 0 : std::min(p1, p2); p3 < L; ++p3)
    //     {
    //       for (pos_t p4 = p3; p4 < L; ++p4) {
    //         std::vector<pos_t> positions{p1, p2, p3, p4};
    //         fun(positions);
    //       }
    //     }
    //   }
    // }
    return fun.get();
  }
};

// Helper class for counting all n-RDM permutations using the
// handle_<n>rdm_indices function
template <class I = Lattice::pos_t, class Dummy = std::vector<I>>
class nrdm_counter {
 public:
  using return_type = I;
  nrdm_counter() : counter_(0) {}
  return_type get() { return counter_; }
  void operator()(const Dummy& d) { counter_++; }

 private:
  return_type counter_;
};

// Helper class for accumulating iterators of indices
template <class I = Lattice::pos_t>
class nrdm_iterator {
 public:
  using vec_type = std::vector<I>;
  using return_type = std::vector<std::vector<I>>;
  nrdm_iterator() : indexes_() {};
  return_type get() { return indexes_; }
  void operator()(const vec_type& vec) { indexes_.push_back(vec); }

 private:
  return_type indexes_;
};

// nice wrapper functions for the above classes
template <int N, class I = Lattice::pos_t>
I get_nrdm_permutations(
    I L, bool bra_neq_ket = false,
    const std::vector<I>& positions_first = std::vector<I>()
) {
  return iterate_rdm_indices<nrdm_counter<I>, N>()(
      nrdm_counter<I>(), L, bra_neq_ket, positions_first
  );
}

template <int N, class I = Lattice::pos_t>
typename nrdm_iterator<I>::return_type iterate_nrdm(
    I L, bool bra_neq_ket = false,
    const std::vector<I>& positions_first = std::vector<I>()
) {
  return iterate_rdm_indices<nrdm_iterator<I>, N>()(
      nrdm_iterator<I>(), L, bra_neq_ket, positions_first
  );
}

}  // namespace measurements_details
#endif
