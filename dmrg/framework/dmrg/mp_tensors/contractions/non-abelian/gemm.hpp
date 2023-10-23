/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef CONTRACTIONS_SU2_GEMM_HPP
#define CONTRACTIONS_SU2_GEMM_HPP

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/indexing.h"

namespace SU2 {

/** @brief Returns the factor correcting the MPO-MPS network for the Hermitian conjugate case */
template <class T, class SymmGroup>
T conjugate_correction(typename SymmGroup::charge lc, typename SymmGroup::charge rc, typename SymmGroup::subcharge tensor_spin)
{
    assert( SymmGroup::spin(lc) >= 0);
    assert( SymmGroup::spin(rc) >= 0);
    typename SymmGroup::subcharge S = std::min(SymmGroup::spin(rc), SymmGroup::spin(lc));
    typename SymmGroup::subcharge spin_diff = SymmGroup::spin(rc) - SymmGroup::spin(lc);
    T ret = 0.;
    if (tensor_spin == 0) {
        ret = T(1.);
    }
    else if (tensor_spin == 1) {
        if (spin_diff > 0)
            ret = -T( sqrt((S + 1.)/(S + 2.)) );
        else if (spin_diff < 0)
            ret = T( sqrt((S + 2.)/(S + 1.)) );
    }
    else if (tensor_spin == 2) {
        if (spin_diff > 0)
            ret = -T( sqrt( (S + 1.) / (S + 3.)) );
        else if (spin_diff < 0)
            ret = -T( sqrt((S + 3.) / (S + 1.)) );
        else
            ret = T(1.);
    }
    else {
        throw std::runtime_error("hermitian conjugate for reduced tensor operators only implemented up to rank 1");
    }
    return ret;
}
/** @brief block_matrix-block_matrix multiplication */
template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
void gemm(block_matrix<Matrix1, SymmGroup> const & A,
          block_matrix<Matrix2, SymmGroup> const & B,
          block_matrix<Matrix3, SymmGroup> & C,
          int spin)
{
    typedef typename SymmGroup::charge charge;
    typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
    typedef typename Matrix3::value_type value_type;
    C.clear();
    assert(B.basis().is_sorted());
    const_iterator B_begin = B.basis().begin();
    const_iterator B_end = B.basis().end();
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        charge ar = A.basis().right_charge(k);
        const_iterator it = B.basis().left_lower_bound(ar);
        for ( ; it != B_end && it->lc == ar; ++it)
        {
            std::size_t matched_block = std::distance(B_begin, it);
            if (!(spin == -1) && !::SU2::triangle(SymmGroup::spin(A.basis().left_charge(k)), spin, SymmGroup::spin(it->rc)))
                continue;
            std::size_t c_block = C.find_block(A.basis().left_charge(k), it->rc);
            if (c_block == C.n_blocks())
                c_block = C.insert_block(Matrix3(num_rows(A[k]), it->rs), A.basis().left_charge(k), it->rc);
            boost::numeric::bindings::blas::gemm(value_type(1), A[k], B[matched_block], value_type(1), C[c_block]);
        }
    }
}

/** @brief block_matrix multipication with left trimming (see abelian version of [gemm_trim_left]) */
template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
void gemm_trim_left(block_matrix<Matrix1, SymmGroup> const & A,
                    block_matrix<Matrix2, SymmGroup> const & B,
                    block_matrix<Matrix3, SymmGroup> & C,
                    Index<SymmGroup> const & refIndex,
                    std::vector<typename Matrix1::value_type> conj_scales = std::vector<typename Matrix1::value_type>())
{
    using charge = typename SymmGroup::charge;
    using const_iterator = typename DualIndex<SymmGroup>::const_iterator;
    using value_type = typename Matrix3::value_type;

    if (conj_scales.size() != A.n_blocks())
        conj_scales = std::vector<value_type>(A.n_blocks(), 1.);

    C.clear();
    assert(B.basis().is_sorted());

    const_iterator B_begin = B.basis().begin();
    const_iterator B_end = B.basis().end();
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {

        //assert(B.basis().left_has(A.basis().left_charge(k)));
        if (!refIndex.has(A.basis().left_charge(k)))
            continue;

        charge ar = A.basis().right_charge(k);
        const_iterator it = B.basis().left_lower_bound(ar);

        for ( ; it != B_end && it->lc == ar; ++it)
        {
            std::size_t matched_block = std::distance(B_begin, it);
            std::size_t c_block = C.find_block(A.basis().left_charge(k), it->rc);
            if (c_block == C.n_blocks())
                c_block = C.insert_block(Matrix3(num_rows(A[k]), it->rs), A.basis().left_charge(k), it->rc);
            boost::numeric::bindings::blas::gemm(conj_scales[k], A[k], B[matched_block], value_type(1), C[c_block]);
        }
    }
}

/** @brief Same as [gemm_trim_right], but does not do any multiplication and returns only the dimensions */
template<class Matrix2, class SymmGroup>
DualIndex<SymmGroup> gemm_trim_right_pretend(DualIndex<SymmGroup> const& A,
                                             block_matrix<Matrix2, SymmGroup> const & B,
                                             Index<SymmGroup> const & refIndex)
{
    typedef typename SymmGroup::charge charge;
    typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
    assert(B.basis().is_sorted());
    DualIndex<SymmGroup> ret;
    const_iterator B_end = B.basis().end();
    for (std::size_t k = 0; k < A.size(); ++k) {
        charge ar = A.right_charge(k);
        const_iterator it = B.basis().left_lower_bound(ar);
        for ( ; it != B_end && it->lc == ar; ++it)
        {
            if (!refIndex.has(it->rc))
                continue;
            charge lc = A.left_charge(k);
            charge rc = it->rc;
            if (!ret.has(lc, rc))
                ret.insert(typename DualIndex<SymmGroup>::value_type(lc, rc, A.left_size(k), it->rs));
        }
    }
    return ret;
}

/** @brief block_matrix multipication with right trimming (see abelian version of [gemm_trim_right]) */
template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
void gemm_trim_right(block_matrix<Matrix1, SymmGroup> const & A,
                     block_matrix<Matrix2, SymmGroup> const & B,
                     block_matrix<Matrix3, SymmGroup> & C,
                     Index<SymmGroup> const & refIndex,
                     std::vector<typename Matrix1::value_type> conj_scales = std::vector<typename Matrix1::value_type>())
{
    typedef typename SymmGroup::charge charge;
    typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
    typedef typename Matrix3::value_type value_type;

    if (conj_scales.size() != B.n_blocks())
        conj_scales = std::vector<value_type>(B.n_blocks(), 1.);
    C.clear();
    assert(B.basis().is_sorted());
    const_iterator B_begin = B.basis().begin();
    const_iterator B_end = B.basis().end();
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        charge ar = A.basis().right_charge(k);
        const_iterator it = B.basis().left_lower_bound(ar);
        for ( ; it != B_end && it->lc == ar; ++it)
        {
            std::size_t matched_block = std::distance(B_begin, it);
            //assert(A.basis().left_has(it->rc));
            if (!refIndex.has(it->rc))
                continue;
            std::size_t c_block = C.find_block(A.basis().left_charge(k), it->rc);
            if (c_block == C.n_blocks())
                c_block = C.insert_block(Matrix3(num_rows(A[k]), it->rs), A.basis().left_charge(k), it->rc);
            boost::numeric::bindings::blas::gemm(conj_scales[matched_block], A[k], B[matched_block], value_type(1), C[c_block]);
        }
    }
}

template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
void gemm_trim(block_matrix<Matrix1, SymmGroup> const & A,
               block_matrix<Matrix2, SymmGroup> const & B,
               block_matrix<Matrix3, SymmGroup> & C,
               std::vector<typename Matrix1::value_type> conj_scales,
               bool conjugate_a)
{
    typedef typename SymmGroup::charge charge;
    typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
    typedef typename Matrix3::value_type value_type;
    assert(B.basis().is_sorted());
    assert( (conjugate_a && A.n_blocks() == conj_scales.size()) || (!conjugate_a && B.n_blocks() == conj_scales.size()));
    C.clear();
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        charge al = A.basis().left_charge(k);
        charge ar = A.basis().right_charge(k);
        std::size_t matched_block = B.basis().position(ar, al);
        if (matched_block == B.n_blocks()) continue;
        std::size_t c_block = C.find_block(al, al);
        if (c_block == C.n_blocks())
            c_block = C.insert_block(Matrix3(num_rows(A[k]), num_cols(B[matched_block])), al, al);
        //boost::numeric::bindings::blas::gemm(value_type(1), A[k], B[matched_block], value_type(1), C[c_block]);
        boost::numeric::bindings::blas::gemm(conj_scales[ (conjugate_a) ? k : matched_block], A[k], B[matched_block],
                                             value_type(1), C[c_block]);
    }
}

} // gemm

#endif
