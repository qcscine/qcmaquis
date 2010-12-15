#include "block_matrix/block_matrix.h"

#include "utils/function_objects.h"
#include <boost/lambda/bind.hpp>

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>::block_matrix(Index<SymmGroup> rows = Index<SymmGroup>(), Index<SymmGroup> cols = Index<SymmGroup>())
: rows_(rows), cols_(cols)
{
    assert(rows_.size() == cols_.size());
    for (size_type k = 0; k < rows_.size(); ++k)
        data_.push_back(Matrix(rows_[k].second, cols_[k].second));
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>::block_matrix(typename SymmGroup::charge c, Matrix const & m)
{
    rows_.push_back(std::make_pair(c, num_rows(m)));
    cols_.push_back(std::make_pair(c, num_columns(m)));
    data_.push_back(m);
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & block_matrix<Matrix, SymmGroup>::operator+=(block_matrix const & rhs)
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

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & block_matrix<Matrix, SymmGroup>::operator-=(block_matrix const & rhs)
{
    for (size_type k = 0; k < rhs.n_blocks(); ++k)
    {
        charge rhs_rc = rhs.rows_[k].first;
        size_type goesto = rows_.position(rhs_rc);
        if (goesto == rows_.size()) { // it's a new block
            size_type i1 = rows_.insert(rhs.rows_[k]);
            size_type i2 = cols_.insert(rhs.cols_[k]);
            assert(i1 == i2);
            data_.insert(data_.begin() + i1, -1*rhs.data_[k]);
        } else { // this block exists already -> pass to Matrix
            data_[goesto] += -1*rhs.data_[k];
        }
    }
    return *this;
}

template<class Matrix, class SymmGroup>
void* block_matrix<Matrix, SymmGroup>::memory_pointer(){

    return NULL;
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::insert_block(boost::tuple<Matrix const &, charge, charge> const & block)
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

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & block_matrix<Matrix, SymmGroup>::operator=(block_matrix rhs)
{
    swap(*this, rhs);
    return *this;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> const & block_matrix<Matrix, SymmGroup>::left_basis() const { return rows_; }

template<class Matrix, class SymmGroup>
Index<SymmGroup> const & block_matrix<Matrix, SymmGroup>::right_basis() const { return cols_; }

template<class Matrix, class SymmGroup>
typename Matrix::size_type block_matrix<Matrix, SymmGroup>::n_blocks() const { return data_.size(); }

template<class Matrix, class SymmGroup>
std::string block_matrix<Matrix, SymmGroup>::description() const
{
    std::ostringstream oss;
    oss << rows_ << cols_;
    return oss.str();
}

template<class Matrix, class SymmGroup>
Matrix & block_matrix<Matrix, SymmGroup>::operator[](size_type c) { return data_[c]; }

template<class Matrix, class SymmGroup>
Matrix const & block_matrix<Matrix, SymmGroup>::operator[](size_type c) const { return data_[c]; }

template<class Matrix, class SymmGroup>
bool block_matrix<Matrix, SymmGroup>::has_block(charge r, charge c)
{
    return rows_.has(r) && cols_.has(c) && rows_.position(r) == cols_.position(c);
}

template<class Matrix, class SymmGroup>
bool block_matrix<Matrix, SymmGroup>::has_block(std::pair<charge, size_type> const & r,
               std::pair<charge, size_type> const & c)
{
    return has_block(r.first, c.first);
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type & block_matrix<Matrix, SymmGroup>::operator()(std::pair<charge, size_type> const & r,
                                                                          std::pair<charge, size_type> const & c)
{
    assert( rows_.position(r.first) == cols_.position(c.first) );
    return data_[rows_.position(r.first)](r.second, c.second);
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type const & block_matrix<Matrix, SymmGroup>::operator()(std::pair<charge, size_type> const & r,
                                                                                std::pair<charge, size_type> const & c) const
{
    assert( rows_.position(r.first) == cols_.position(c.first) );
    return data_[rows_.position(r.first)](r.second, c.second);
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::remove_rows_from_block(size_type block, size_type r, size_type k = 1)
{
    remove_rows(data_[block], r, k);
    rows_[block].second -= k;
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::remove_cols_from_block(size_type block, size_type r, size_type k = 1)
{
    remove_columns(data_[block], r, k);
    cols_[block].second -= k;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & block_matrix<Matrix, SymmGroup>::operator*=(value_type v)
{
    std::for_each(data_.begin(), data_.end(), boost::lambda::_1 *= v);
    return *this;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & block_matrix<Matrix, SymmGroup>::operator/=(value_type v)
{
    std::for_each(data_.begin(), data_.end(), boost::lambda::_1 /= v);
    return *this;
}


template<class Matrix, class SymmGroup>
typename Matrix::value_type block_matrix<Matrix, SymmGroup>::trace() const
{
    std::vector<value_type> vt(n_blocks());
    std::transform(data_.begin(), data_.end(), vt.begin(), utils::functor_trace());
    return std::accumulate(vt.begin(), vt.end(), value_type());
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::inplace_transpose()
{
    std::transform(data_.begin(), data_.end(), data_.begin(), utils::functor_transpose());
    std::swap(rows_, cols_);
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::inplace_conjugate()
{
    std::transform(data_.begin(), data_.end(), data_.begin(), utils::functor_conjugate());
}

template<class Matrix, class SymmGroup>
template<class Generator>
void block_matrix<Matrix, SymmGroup>::generate(Generator g)
{
    for (std::size_t k = 0; k < n_blocks(); ++k)
        std::generate(elements(data_[k]).first, elements(data_[k]).second, g);
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::clear()
{
    data_.clear();
    rows_ = Index<SymmGroup>();
    cols_ = Index<SymmGroup>();
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, block_matrix<Matrix, SymmGroup> const & m)
{
    os << "Left HS: " << m.left_basis() << std::endl;
    os << "Right HS: " << m.right_basis() << std::endl;
    for (std::size_t k = 0; k < m.n_blocks(); ++k)
        os << "Block (" << m.left_basis()[k].first << "," << m.right_basis()[k].first << "):" << std::endl << m[k];
    os << std::endl;
    return os;
}
