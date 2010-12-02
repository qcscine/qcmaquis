#ifndef BLOCK_MATRIX_H
#define BLOCK_MATRIX_H

#include <sstream>
#include <algorithm>
#include <numeric>

#include <boost/tuple/tuple.hpp>
#include <boost/lambda/bind.hpp>

#include "indexing.h"
#include "symmetry.h"

#include "function_objects.h"

template<class Matrix, class SymmGroup>
class block_matrix
{
private:
    typedef typename SymmGroup::charge charge;
    
public:
    typedef typename Matrix::size_type size_type;
    typedef typename Matrix::value_type value_type;
    
    block_matrix(Index<SymmGroup> rows = Index<SymmGroup>(), Index<SymmGroup> cols = Index<SymmGroup>())
    : rows_(rows)
    , cols_(cols)
    {
        assert(rows_.size() == cols_.size());
        for (size_type k = 0; k < rows_.size(); ++k)
            data_.push_back(Matrix(rows_[k].second, cols_[k].second));
    }
    
    block_matrix(charge c, Matrix const & m)
    {
        rows_.push_back(std::make_pair(c, num_rows(m)));
        cols_.push_back(std::make_pair(c, num_columns(m)));
        data_.push_back(m);
    }
    
    block_matrix & operator+=(block_matrix const & rhs)
    {
        for (size_type k = 0; k < rhs.n_blocks(); ++k)
        {
            charge rhs_rc = rhs.rows_[k].first;
            size_type goesto = rows_.position(rhs_rc);
            if (goesto == rows_.size()) { // it's a new block
                size_type i1 = rows_.insert(rhs.rows_[k]);
                size_type i2 = cols_.insert(rhs.cols_[k]);
                assert(i1 == i2);
                data_.insert(data_.begin() + i1, rhs.data_[k]);
            } else { // this block exists already -> pass to Matrix
                data_[goesto] += rhs.data_[k];
            }
        }
        return *this;
    }
    
    void insert_block(boost::tuple<Matrix, charge, charge> const & block)
    {
        Matrix const & mtx = boost::tuples::get<0>(block);
        
        std::pair<charge, size_type>
        p1 = std::make_pair(boost::tuples::get<1>(block), mtx.num_rows()),
        p2 = std::make_pair(boost::tuples::get<2>(block), mtx.num_columns());
        
        size_type i1 = rows_.insert(p1);
        size_type i2 = cols_.insert(p2);
        assert(i1 == i2);
        data_.insert(data_.begin() + i1, mtx);
    }
    
    friend void swap(block_matrix & x, block_matrix & y)
    {
        std::swap(x.data_, y.data_);
        std::swap(x.rows_, y.rows_);
        std::swap(x.cols_, y.cols_);
    }
    
    block_matrix & operator=(block_matrix rhs)
    {
        swap(*this, rhs);
        return *this;
    }
    
    Index<SymmGroup> left_basis() const { return rows_; }
    Index<SymmGroup> right_basis() const { return cols_; }
    
    size_type n_blocks() const { return data_.size(); }
    
    std::string description() const
    {
        std::ostringstream oss;
        oss << rows_ << cols_;
        return oss.str();
    }
    
    Matrix & operator[](size_type c) { return data_[c]; }
    Matrix const & operator[](size_type c) const { return data_[c]; }
    
    value_type & operator()(std::pair<charge, size_type> const & r,
                            std::pair<charge, size_type> const & c)
    {
        assert( rows_.position(r.first) == cols_.position(c.first) );
        return data_[rows_.position(r.first)](r.second, c.second);
    }
    
    value_type const & operator()(std::pair<charge, size_type> const & r,
                                  std::pair<charge, size_type> const & c) const
    {
        assert( rows_.position(r.first) == cols_.position(c.first) );
        return data_[rows_.position(r.first)](r.second, c.second);
    }
    
    void remove_rows_from_block(size_type block, size_type r, size_type k = 1)
    {
        remove_rows(data_[block], r, k);
        rows_[block].second -= k;
    }
    
    void remove_cols_from_block(size_type block, size_type r, size_type k = 1)
    {
        remove_columns(data_[block], r, k);
        cols_[block].second -= k;
    }
    
    block_matrix operator*=(value_type v)
    {
        std::for_each(data_.begin(), data_.end(), boost::lambda::_1 *= v);
        return *this;
    }
    
    value_type trace() const
    {
        std::vector<value_type> vt(n_blocks());
        std::transform(data_.begin(), data_.end(), vt.begin(), functors::ftrace());
        return std::accumulate(vt.begin(), vt.end(), value_type());
    }
    
    void inplace_transpose()
    {
        std::transform(data_.begin(), data_.end(), data_.begin(), functors::ftranspose());
        std::swap(rows_, cols_);
    }
    
    void inplace_conjugate()
    {
        std::transform(data_.begin(), data_.end(), data_.begin(), functors::fconjugate());
    }
    
    template<class Generator>
    void fill_with_random(Generator & g)
    {
        for (std::size_t k = 0; k < n_blocks(); ++k)
            std::generate(elements(data_[k]).first, elements(data_[k]).second, g);
    }
    
private:
    std::vector<Matrix> data_;
    Index<SymmGroup> rows_, cols_;
};    

template<class Matrix, class SymmGroup>    
std::ostream& operator<<(std::ostream& os, block_matrix<Matrix, SymmGroup> const & m)
{
    os << "Left HS: " << std::endl << m.left_basis();
    os << "Right HS: " << std::endl << m.right_basis();
    for (std::size_t k = 0; k < m.n_blocks(); ++k)
        os << "Block (" << m.left_basis()[k].first << "," << m.right_basis()[k].first << "):" << std::endl << m[k];
    os << std::endl;
    return os;
}

#endif
