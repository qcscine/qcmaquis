/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MPOTENSOR_H
#define MPOTENSOR_H

#include <iostream>
#include <set>
#include <iterator>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/indexing.h"
#include "utils/function_objects.h"
#include "dmrg/models/OperatorHandlers/OpTable.h"
#include "dmrg/mp_tensors/mpotensor_detail.h"

template<class Matrix, class SymmGroup>
class MPOTensor
{
public:
    using index_type = std::size_t;
    using value_type = typename Matrix::value_type;
    using scalar_type = typename maquis::traits::scalar_type<Matrix>::type;

    using tag_type = typename OPTable<Matrix, SymmGroup>::tag_type;
    using op_t = typename OPTable<Matrix, SymmGroup>::op_t;
    using op_table_ptr = std::shared_ptr<OPTable<Matrix, SymmGroup>>;
    using pv_type = std::pair<tag_type, value_type>;

    using internal_value_type = std::vector<pv_type>;

private:
    using CSCMatrix = boost::numeric::ublas::compressed_matrix<internal_value_type, boost::numeric::ublas::column_major, 0>;

    using RowIndex = std::vector<std::set<index_type>>;

public:
    using row_proxy = MPOTensor_detail::row_proxy<typename RowIndex::value_type::const_iterator>;
    using col_proxy = boost::numeric::ublas::matrix_column<const CSCMatrix>;

    using prempo_t = std::vector<boost::tuple<std::size_t, std::size_t, tag_type, value_type>>;
    using spin_desc_t = SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type>;
    using spin_index = std::vector<spin_desc_t>;

public:
    MPOTensor(index_type = 1, index_type = 1, prempo_t = prempo_t(), op_table_ptr = op_table_ptr(),
              MPOTensor_detail::Hermitian = MPOTensor_detail::Hermitian(1,1),
              spin_index const & lspins = spin_index(), spin_index const & rspins = spin_index());

    index_type row_dim() const;
    index_type col_dim() const;

    // tagged operator ()
    // warning: this method allows to (indirectly) change the op in the table, all tags pointing to it will
    //          get a modified matrix!
    //          better design needed
    void set(index_type li, index_type ri, op_t const & op, value_type scale_ = 1.0);

    // tagged operator() const
    MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true>
    at(index_type left_index, index_type right_index) const;

    // warning: this method allows to (indirectly) change the op in the table, all tags pointing to it will
    //          get a modified matrix!
    //          better design needed
    MPOTensor_detail::term_descriptor<Matrix, SymmGroup, false>
    at(index_type left_index, index_type right_index);

    row_proxy row(index_type row_i) const;
    col_proxy column(index_type col_i) const;

    tag_type tag_number(index_type left_index, index_type right_index, size_t index = 0) const;
    op_table_ptr get_operator_table() const;

    void multiply_by_scalar(value_type);
    void divide_by_scalar(value_type);

    bool has(index_type left_index, index_type right_index) const;

    spin_desc_t left_spin(index_type left_index) const;
    spin_desc_t right_spin(index_type right_index) const;
    spin_index const & row_spin_dim() const;
    spin_index const & col_spin_dim() const;
    index_type num_row_non_zeros(index_type row_i) const;
    index_type num_col_non_zeros(index_type col_i) const;
    index_type num_one_rows() const;
    index_type num_one_cols() const;
    MPOTensor_detail::Hermitian herm_info;

private:
    index_type left_i, right_i;
    spin_index left_spins, right_spins;
    std::vector<index_type> row_non_zeros, col_non_zeros;
    index_type num_one_rows_, num_one_cols_;
    CSCMatrix col_tags;
    RowIndex row_index;
    op_table_ptr operator_table;
};

// TODO: add swap

#include "dmrg/mp_tensors/mpotensor.hpp"


#endif
