/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef BLOCK_MATRIX_H
#define BLOCK_MATRIX_H

#include <sstream>
#include <algorithm>
#include <numeric>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"

#include "utils/timings.h"
#include "types/utils/traits.hpp"

#include <boost/ptr_container/ptr_vector.hpp>

template<class Matrix, class SymmGroup>
class block_matrix
{
private:
    typedef typename SymmGroup::charge charge;
public:
    typedef typename Matrix::size_type size_type;
    typedef typename Matrix::value_type value_type;
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    
    block_matrix(Index<SymmGroup> rows = Index<SymmGroup>(),
                 Index<SymmGroup> cols = Index<SymmGroup>());

    Index<SymmGroup> const & left_basis() const;
    Index<SymmGroup> const & right_basis() const;
    
    block_matrix(charge rc, charge cc, Matrix const & m)
    {
        rows_.push_back(std::make_pair(rc, num_rows(m)));
        cols_.push_back(std::make_pair(cc, num_columns(m)));
        data_.push_back(m);
    }
    std::string description() const;
    
    block_matrix &       operator=(block_matrix rhs);
    Matrix &             operator[](size_type c);
    Matrix const &       operator[](size_type c) const;
    value_type &         operator()(std::pair<charge, size_type> const & r,
                                    std::pair<charge, size_type> const & c);
    value_type const &   operator()(std::pair<charge, size_type> const & r,
                                    std::pair<charge, size_type> const & c) const;
    block_matrix &       operator+=(block_matrix const & rhs);
    block_matrix &       operator-=(block_matrix const & rhs);
    block_matrix const & operator*=(scalar_type v);
    block_matrix const & operator/=(scalar_type v);

    size_type n_blocks() const;
    bool has_block(charge r, charge c) const;
    bool has_block(std::pair<charge, size_type> const & r,
                   std::pair<charge, size_type> const & c) const;
    
    void insert_block(Matrix const &, charge, charge);
    void remove_block(charge r, charge c);
    void remove_block(std::size_t which);

    void remove_rows_from_block(size_type block, size_type r, size_type k = 1);
    void remove_cols_from_block(size_type block, size_type r, size_type k = 1);
    
    scalar_type trace() const;
    void inplace_transpose();
    void inplace_conjugate();
    void clear();
    template<class Generator>
    void generate(Generator g);

    void match_and_add_block(Matrix const &, charge, charge);
    
    void reserve(charge, charge, std::size_t, std::size_t);
    inline void reserve_pos(charge, charge, std::size_t, std::size_t);
    void allocate_blocks();
    
    void resize_block(charge r, charge c,
                      size_type new_r, size_type new_c,
                      bool pretend = false);
    
    friend void swap(block_matrix & x, block_matrix & y)
    {
        swap(x.data_, y.data_);
        swap(x.rows_, y.rows_);
        swap(x.cols_, y.cols_);
    }
    
    Matrix const & operator()(charge r, charge c) const
    {
        assert( has_block(r, c) );
        assert( rows_.position(r) == cols_.position(c) );
        return data_[rows_.position(r)];
    }
    
    Matrix & operator()(charge r, charge c)
    {
        assert( has_block(r, c) );
        assert( rows_.position(r) == cols_.position(c) );
        return data_[rows_.position(r)];
    }
    
#ifdef HAVE_ALPS_HDF5
    void load(alps::hdf5::archive & ar);
    void save(alps::hdf5::archive & ar) const;
#endif
    
private:
    boost::ptr_vector<Matrix> data_;
    Index<SymmGroup> rows_, cols_;
};    

#include "dmrg/block_matrix/block_matrix.hpp"
template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> operator*(typename block_matrix<Matrix,SymmGroup>::scalar_type v,
                                          block_matrix<Matrix, SymmGroup> bm)
{
    bm *= v;
    return bm;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> operator*(block_matrix<Matrix, SymmGroup> bm,
                                          typename block_matrix<Matrix,SymmGroup>::scalar_type v)
{
    bm *= v;
    return bm;
}

template<class Matrix, class SymmGroup>
std::size_t size_of(block_matrix<Matrix, SymmGroup> const & m)
{
    size_t r = 0;
    for (size_t i = 0; i < m.n_blocks(); ++i)
        r += size_of(m[i]);
    return r;
}


#endif
