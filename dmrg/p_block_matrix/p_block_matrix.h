#ifndef BLOCK_MATRIX_H
#define BLOCK_MATRIX_H

#include <sstream>
#include <algorithm>
#include <numeric>

#include "p_block_matrix/indexing.h"
#include "p_block_matrix/symmetry.h"
#include "ambient/interfaces/i_block_matrix.h"


template<class Matrix, class SymmGroup>
class p_block_matrix : public i_block_matrix
{
private:
    typedef typename SymmGroup::charge charge;
public:
    typedef typename Matrix::size_type size_type;
    typedef typename Matrix::value_type value_type;
    
    p_block_matrix(Index<SymmGroup> rows, Index<SymmGroup> cols);

    /**
      * The i_block_matrix implementation:
      * @returns memory pointer to the actual matrix memory
      *
      */
    void* memory_pointer();

    Index<SymmGroup> const & left_basis() const;
    Index<SymmGroup> const & right_basis() const;
    
    p_block_matrix(charge rc, charge cc, Matrix const & m)
    {
        rows_.push_back(std::make_pair(rc, num_rows(m)));
        cols_.push_back(std::make_pair(cc, num_columns(m)));
        data_.push_back(m);
    }
    std::string description() const;
    
    p_block_matrix &     operator=(p_block_matrix rhs);
    Matrix &             operator[](size_type c);
    Matrix const &       operator[](size_type c) const;
    value_type &         operator()(std::pair<charge, size_type> const & r,
                                    std::pair<charge, size_type> const & c);
    value_type const &   operator()(std::pair<charge, size_type> const & r,
                                    std::pair<charge, size_type> const & c) const;
    p_block_matrix &       operator+=(p_block_matrix const & rhs);
    p_block_matrix &       operator-=(p_block_matrix const & rhs);
    p_block_matrix const & operator*=(value_type v);
    p_block_matrix const & operator/=(value_type v);

    size_type n_blocks() const;
    bool has_block(charge r, charge c) const;
    bool has_block(std::pair<charge, size_type> const & r,
                   std::pair<charge, size_type> const & c) const;
    
    void insert_block(Matrix const &, charge, charge);
    void remove_block(charge r, charge c);
    
    void remove_rows_from_block(size_type block, size_type r, size_type k);
    void remove_cols_from_block(size_type block, size_type r, size_type k);
    
    value_type trace() const;
    void inplace_transpose();
    void inplace_conjugate();
    void clear();
    template<class Generator>
    void generate(Generator g);

    void match_and_add_block(Matrix const &, charge, charge);
    
    void reserve(charge, charge, std::size_t, std::size_t);
    void allocate_blocks();
    
    void resize_block(charge r, charge c,
                      size_type new_r, size_type new_c,
                      bool pretend = false);
    
    friend void swap(p_block_matrix & x, p_block_matrix & y)
    {
        swap(x.data_, y.data_);
        swap(x.rows_, y.rows_);
        swap(x.cols_, y.cols_);
    }
    
    Matrix const & operator()(charge r, charge c) const
    {
        assert( rows_.position(r) == cols_.position(c) );
        return data_[rows_.position(r)];
    }
    
    Matrix & operator()(charge r, charge c)
    {
        assert( rows_.position(r) == cols_.position(c) );
        return data_[rows_.position(r)];
    }
    
#ifdef HAVE_ALPS_HDF5
    void serialize(alps::hdf5::iarchive & ar);
    void serialize(alps::hdf5::oarchive & ar) const;
#endif
    
private:
    std::vector<Matrix> data_;
    Index<SymmGroup> rows_, cols_;
};    

#include "p_block_matrix/p_block_matrix.hpp"
template<class Matrix, class SymmGroup>
p_block_matrix<Matrix, SymmGroup> operator*(typename Matrix::value_type v,
                                            p_block_matrix<Matrix, SymmGroup> bm)
{
    bm *= v;
    return bm;
}

template<class Matrix, class SymmGroup>
p_block_matrix<Matrix, SymmGroup> operator*(p_block_matrix<Matrix, SymmGroup> bm,
                                            typename Matrix::value_type v)
{
    bm *= v;
    return bm;
}

#endif
