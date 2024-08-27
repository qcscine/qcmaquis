/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef MPOTENSOR_DETAIL_H
#define MPOTENSOR_DETAIL_H

#include <numeric>
#include <utility>
#include <vector>

#include <boost/utility.hpp>
#include <boost/type_traits.hpp>

#include "dmrg/models/OperatorHandlers/OpTable.h"

template <class Matrix, class SymmGroup>
class MPOTensor;

namespace MPOTensor_detail {
template <class T, bool C>
struct const_type {
  using type = T;
};

template <class T>
struct const_type<T, true> {
  using type = const T;
};

template <class Matrix, class SymmGroup, bool Const>
class term_descriptor {
  using value_type = typename Matrix::value_type;
  using op_t = typename OPTable<Matrix, SymmGroup>::op_t;
  using tag_type = typename OPTable<Matrix, SymmGroup>::tag_type;
  using internal_value_type =
      typename MPOTensor<Matrix, SymmGroup>::internal_value_type;
  using op_table_ptr = typename MPOTensor<Matrix, SymmGroup>::op_table_ptr;

 public:
  term_descriptor() = default;
  term_descriptor(
      typename const_type<internal_value_type, Const>::type &term_descs,
      op_table_ptr op_tbl_
  )
      : term_descriptors(term_descs), operator_table(op_tbl_) {}

  std::size_t size() const { return term_descriptors.size(); }
  typename const_type<op_t, Const>::type &op(std::size_t i = 0) {
    return (*operator_table)[term_descriptors[i].first];
  }
  typename const_type<value_type, Const>::type &scale(std::size_t i = 0) {
    return term_descriptors[i].second;
  }

 private:
  typename const_type<internal_value_type, Const>::type &term_descriptors;
  op_table_ptr operator_table;
};

template <class ConstIterator>
class IteratorWrapper {
  using internal_iterator = ConstIterator;

 public:
  using iterator_category = std::forward_iterator_tag;
  using self_type = IteratorWrapper<ConstIterator>;
  using value_type =
      typename std::iterator_traits<internal_iterator>::value_type;

  IteratorWrapper(internal_iterator i) : it_(i) {}

  void operator++() { ++it_; }
  void operator++(int) { it_++; }
  bool operator!=(self_type const &rhs) { return it_ != rhs.it_; }

  value_type index() const { return *it_; }
  value_type operator*() const {
    throw std::runtime_error(
        "direct MPOTensor access via row iterators currently not implemented\n"
    );
    return *it_;
  }

 private:
  internal_iterator it_;
};

template <class ConstIterator>
class row_proxy : public std::pair<ConstIterator, ConstIterator> {
  using internal_iterator = ConstIterator;
  using base = std::pair<internal_iterator, internal_iterator>;

 public:
  using const_iterator = IteratorWrapper<ConstIterator>;
  row_proxy(internal_iterator b, internal_iterator e) : base(b, e) {}

  const_iterator begin() const { return const_iterator(base::first); }
  const_iterator end() const { return const_iterator(base::second); }
};

template <class Tuple>
struct row_cmp {
  bool operator()(Tuple const &i, Tuple const &j) const {
    if (std::get<0>(i) == std::get<0>(j)) {
      return std::get<1>(i) < std::get<1>(j);
    } else {
      return std::get<0>(i) < std::get<0>(j);
    }
  }
};

template <class Tuple>
struct col_cmp {
  bool operator()(Tuple const &i, Tuple const &j) const {
    if (std::get<1>(i) == std::get<1>(j)) {
      return std::get<0>(i) < std::get<0>(j);
    } else {
      return std::get<1>(i) < std::get<1>(j);
    }
  }
};

class Hermitian {
  using index_type = std::size_t;

  friend Hermitian operator*(Hermitian const &, Hermitian const &);

 public:
  Hermitian(index_type ld, index_type rd) {
    LeftHerm.resize(ld);
    RightHerm.resize(rd);
    LeftPhase = std::vector<int>(ld, 1);
    RightPhase = std::vector<int>(rd, 1);

    std::iota(LeftHerm.begin(), LeftHerm.end(), static_cast<index_type>(0));
    std::iota(RightHerm.begin(), RightHerm.end(), static_cast<index_type>(0));
  }

  Hermitian(
      std::vector<index_type> lh, std::vector<index_type> rh,
      std::vector<int> lp, std::vector<int> rp
  )
      : LeftHerm(std::move(lh)),
        RightHerm(std::move(rh)),
        LeftPhase(std::move(lp)),
        RightPhase(std::move(rp)) {}

  bool left_skip(index_type b1) const { return LeftHerm[b1] < b1; }
  bool right_skip(index_type b2) const { return RightHerm[b2] < b2; }

  index_type left_conj(index_type b1) const { return LeftHerm[b1]; }
  index_type right_conj(index_type b2) const { return RightHerm[b2]; }

  std::size_t left_size() const { return LeftHerm.size(); }
  std::size_t right_size() const { return RightHerm.size(); }

  int left_phase(std::size_t i) const { return LeftPhase[i]; }
  int right_phase(std::size_t i) const { return RightPhase[i]; }

 private:
  std::vector<index_type> LeftHerm;
  std::vector<index_type> RightHerm;

  std::vector<int> LeftPhase;
  std::vector<int> RightPhase;
};

inline Hermitian operator*(Hermitian const &a, Hermitian const &b) {
  return Hermitian(a.LeftHerm, b.RightHerm, a.LeftPhase, b.RightPhase);
}

template <class Matrix, class SymmGroup>
symm_traits::disable_if_su2_t<SymmGroup, int> get_spin(
    MPOTensor<Matrix, SymmGroup> const &mpo,
    typename MPOTensor<Matrix, SymmGroup>::index_type k, bool left
) {
  return 0;
}

template <class Matrix, class SymmGroup>
symm_traits::enable_if_su2_t<SymmGroup, int> get_spin(
    MPOTensor<Matrix, SymmGroup> const &mpo,
    typename MPOTensor<Matrix, SymmGroup>::index_type k, bool left
) {
  if (left) {
    return mpo.left_spin(k).get();
  } else {
    return mpo.right_spin(k).get();
  }
}
}  // namespace MPOTensor_detail

#endif
