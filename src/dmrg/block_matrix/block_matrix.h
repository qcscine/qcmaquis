#ifndef BLOCK_MATRIX_H
#define BLOCK_MATRIX_H

#include <sstream>
#include <algorithm>
#include <numeric>

#include <boost/tuple/tuple.hpp>

#include "block_matrix/indexing.h"
#include "block_matrix/symmetry.h"
#include "block_matrix/i_block_matrix.h"


template<class Matrix, class SymmGroup>
class block_matrix : implements i_block_matrix
{
private:
    typedef typename SymmGroup::charge charge;
public:
    typedef typename Matrix::size_type size_type;
    typedef typename Matrix::value_type value_type;
    
    block_matrix(Index<SymmGroup> rows, Index<SymmGroup> cols);
    block_matrix(charge c, Matrix const & m);

    /**
      * The i_block_matrix implementation:
      * @returns memory pointer to the actual matrix memory
      *
      */
    void* memory_pointer();

    Index<SymmGroup> const & left_basis() const;
    Index<SymmGroup> const & right_basis() const;
    
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
    block_matrix const & operator*=(value_type v);
    block_matrix const & operator/=(value_type v);

    size_type n_blocks() const;
    bool has_block(charge r, charge c);
    bool has_block(std::pair<charge, size_type> const & r,
                   std::pair<charge, size_type> const & c);
    void insert_block(boost::tuple<Matrix const &, charge, charge> const & block);
    void remove_rows_from_block(size_type block, size_type r, size_type k);
    void remove_cols_from_block(size_type block, size_type r, size_type k);
    value_type trace() const;
    void inplace_transpose();
    void inplace_conjugate();
    void clear();
    template<class Generator>
    void generate(Generator g);

    friend void swap(block_matrix & x, block_matrix & y)
    {
        swap(x.data_, y.data_);
        swap(x.rows_, y.rows_);
        swap(x.cols_, y.cols_);
    }
    
private:
    std::vector<Matrix> data_;
    Index<SymmGroup> rows_, cols_;
};    

#include "block_matrix/block_matrix.hpp"
#endif
