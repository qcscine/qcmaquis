/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

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
#include "dmrg/models/tag_table.h"
#include "dmrg/mp_tensors/mpotensor_detail.h"


template <class Matrix, class SymmGroup> class column_iterator;
template <class Matrix, class SymmGroup> class compressor;
template <class Matrix, class SymmGroup> class MPOIndexer;

template<class Matrix, class SymmGroup>
class MPOTensor
{
public:
    typedef std::size_t index_type;
    typedef typename Matrix::value_type value_type;
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;

    typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
    typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
    typedef boost::shared_ptr<OPTable<Matrix, SymmGroup> > op_table_ptr;

private:
    typedef std::pair<tag_type, value_type> internal_value_type;

    typedef boost::numeric::ublas::compressed_matrix< internal_value_type,
                                                      boost::numeric::ublas::column_major
                                                      , 0, boost::numeric::ublas::unbounded_array<index_type> 
                                                    > CSCMatrix;

    typedef std::vector<std::set<index_type> > RowIndex;
    
public:
    typedef MPOTensor_detail::row_proxy<typename RowIndex::value_type::const_iterator> row_proxy;
    typedef boost::numeric::ublas::matrix_column<const CSCMatrix> col_proxy;

private:
    typedef std::vector<boost::tuple<std::size_t, std::size_t, tag_type, value_type> > prempo_t;
    
public:
    MPOTensor(index_type = 1, index_type = 1, prempo_t const & = prempo_t(), op_table_ptr = op_table_ptr());
    
    index_type row_dim() const;
    index_type col_dim() const;
    
    //block_matrix<Matrix, SymmGroup> const & operator()(index_type left_index,
    //                                                   index_type right_index) const;
    //block_matrix<Matrix, SymmGroup> & operator()(index_type left_index,
    //                                             index_type right_index);

    // tagged operator ()
    void set(index_type li, index_type ri, op_t const & op, value_type scale_ = 1.0);

    // tagged operator() const
    MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup>
    at(index_type left_index, index_type right_index) const;

    row_proxy row(index_type row_i) const;
    col_proxy column(index_type col_i) const;

    tag_type tag_number(index_type left_index, index_type right_index) const;
    op_table_ptr get_operator_table() const;

    void multiply_by_scalar(const scalar_type&);
    void divide_by_scalar(const scalar_type&);
    
    bool has(index_type left_index, index_type right_index) const;

    // These will be removed soon
    friend class column_iterator<Matrix, SymmGroup>;
    friend class compressor<Matrix, SymmGroup>;
    friend class MPOIndexer<Matrix, SymmGroup>;
    
private:
    CSCMatrix col_tags;
    RowIndex row_index;
    op_table_ptr operator_table;
    
    index_type left_i, right_i;
};

#include "dmrg/mp_tensors/mpotensor.hpp"


#endif
