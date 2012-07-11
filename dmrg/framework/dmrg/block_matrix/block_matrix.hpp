/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/block_matrix/block_matrix.h"

#include "utils/function_objects.h"

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>::block_matrix() 
{
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>::block_matrix(Index<SymmGroup> rows,
                                              Index<SymmGroup> cols)
: rows_(rows)
, cols_(cols)
{
    assert(rows_.size() == cols_.size());
    for (size_type k = 0; k < rows_.size(); ++k)
        data_.push_back(new Matrix(rows_[k].second, cols_[k].second));
}

// Remove by Tim 06/08/2012, presently not used in any DMRG/TE code
//template<class Matrix, class SymmGroup>
//block_matrix<Matrix, SymmGroup>::block_matrix(charge rc, charge cc, Matrix& m)
//{
//    rows_.push_back(std::make_pair(rc, num_rows(m)));
//    cols_.push_back(std::make_pair(cc, num_cols(m)));
//    data_.push_back(&m);
//}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & block_matrix<Matrix, SymmGroup>::operator+=(block_matrix const & rhs)
{
    for (size_type k = 0; k < rhs.n_blocks(); ++k)
    {
        charge rhs_rc = rhs.rows_[k].first;
        charge rhs_cc = rhs.cols_[k].first;
        if (this->has_block(rhs_rc, rhs_cc))
            (*this)(rhs_rc, rhs_cc) += rhs.data_[k];
        else
            this->insert_block(rhs.data_[k], rhs_rc, rhs_cc);
    }
    return *this;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & block_matrix<Matrix, SymmGroup>::operator-=(block_matrix const & rhs)
{
    for (size_type k = 0; k < rhs.n_blocks(); ++k)
    {
        charge rhs_rc = rhs.rows_[k].first;
        charge rhs_cc = rhs.cols_[k].first;
        if (this->has_block(rhs_rc, rhs_cc))
            (*this)(rhs_rc, rhs_cc) -= rhs.data_[k];
        else
            this->insert_block(-1*rhs.data_[k], rhs_rc, rhs_cc);
    }
    return *this;
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::insert_block(Matrix const & mtx, charge c1, charge c2)
{
    assert( !has_block(c1, c2) );
    
    std::pair<charge, size_type>
    p1 = std::make_pair(c1, num_rows(mtx)),
    p2 = std::make_pair(c2, num_cols(mtx));
    
    size_type i1 = rows_.insert(p1);
    cols_.insert(i1, p2);
    data_.insert(data_.begin() + i1, new Matrix(mtx));
    
    //rows_.push_back(p1);
    //cols_.push_back(p2);
    //data_.push_back(mtx);
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::insert_block(Matrix * mtx, charge c1, charge c2)
{
    assert( !has_block(c1, c2) );
    
    std::pair<charge, size_type>
    p1 = std::make_pair(c1, num_rows(*mtx)),
    p2 = std::make_pair(c2, num_cols(*mtx));
    
    size_type i1 = rows_.insert(p1);
    cols_.insert(i1, p2);
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
bool block_matrix<Matrix, SymmGroup>::has_block(charge r, charge c) const
{
    bool ret;
    std::size_t p1 = rows_.position(r);
    if (p1 == rows_.size())
        ret = false;
    else {
        std::size_t p2 = cols_.position(c);
        if (p2 == cols_.size())
            ret = false;
        else
            ret = (p1 == p2);
    }
    return ret;
}

template<class Matrix, class SymmGroup>
bool block_matrix<Matrix, SymmGroup>::has_block(std::pair<charge, size_type> const & r,
                                                std::pair<charge, size_type> const & c) const
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
// Remove by Tim 06/08/2012, presently not used in any DMRG/TE code
//template<class Matrix, class SymmGroup>
//void block_matrix<Matrix, SymmGroup>::remove_rows_from_block(size_type block, size_type r, size_type k)
//{ // we should add an assert block < data_.size()
//    remove_rows(data_[block], r, k);
//    rows_[block].second -= k;
//}

// Remove by Tim 06/08/2012, presently not used in any DMRG/TE code
//template<class Matrix, class SymmGroup>
//void block_matrix<Matrix, SymmGroup>::remove_cols_from_block(size_type block, size_type r, size_type k)
//{ // we should add an assert block < data_.size()
//    remove_cols(data_[block], r, k);
//    cols_[block].second -= k;
//}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & block_matrix<Matrix, SymmGroup>::operator*=(const scalar_type& v)
{
    for(size_type k = 0; k < n_blocks(); ++k) data_[k] *= v;
    return *this;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & block_matrix<Matrix, SymmGroup>::operator/=(const scalar_type& v)
{
    for(size_type k = 0; k < n_blocks(); ++k) data_[k] /= v;
    return *this;
}


template<class Matrix, class SymmGroup>
typename block_matrix<Matrix, SymmGroup>::scalar_type block_matrix<Matrix, SymmGroup>::trace() const
{
    std::vector<scalar_type> vt;
    std::transform(data_.begin(), data_.end(), back_inserter(vt), utils::functor_trace());
    return std::accumulate(vt.begin(), vt.end(), scalar_type(0)); // maquis::accumulate
}

template<class Matrix, class SymmGroup>
typename block_matrix<Matrix, SymmGroup>::real_type block_matrix<Matrix, SymmGroup>::norm() const
{
    std::vector<real_type> vt;
    std::transform(data_.begin(), data_.end(), back_inserter(vt), utils::functor_norm_square());
    return sqrt(std::accumulate(vt.begin(), vt.end(), real_type(0)));
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::transpose_inplace()
{
    std::for_each(data_.begin(), data_.end(), utils::functor_transpose_inplace());
    std::swap(rows_, cols_);
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::conjugate_inplace()
{
    std::for_each(data_.begin(), data_.end(), utils::functor_conj_inplace());
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::adjoint_inplace()
{
    std::for_each(data_.begin(), data_.end(), utils::functor_adjoint_inplace());
    std::swap(rows_, cols_);
}

namespace detail {
    // md: this is needed because block_matrix::generate is hiding the correct one
    template<class Matrix, class Generator>
    void generate_impl(Matrix & m, Generator g)
    { generate(m, g); }
}

template<class Matrix, class SymmGroup>
template<class Generator>
void block_matrix<Matrix, SymmGroup>::generate(Generator g)
{
    for(std::size_t k = 0; k < n_blocks(); ++k) detail::generate_impl(data_[k], g);
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

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::match_and_add_block(Matrix const & mtx, charge c1, charge c2)
{
    
    if (this->has_block(c1, c2))
    {
        if (num_rows(mtx) == num_rows((*this)(c1, c2)) &&
            num_cols(mtx) == num_cols((*this)(c1, c2)))
            (*this)(c1, c2) += mtx;
        else if (num_rows(mtx) > num_rows((*this)(c1, c2)) &&
                 num_cols(mtx) > num_cols((*this)(c1, c2)))
        {
            resize_block(c1, c2, num_rows(mtx), num_cols(mtx));
            (*this)(c1, c2) += mtx;
        } else {
            std::size_t maxrows = std::max(num_rows(mtx),
                                           num_rows((*this)(c1, c2)));
            std::size_t maxcols = std::max(num_cols(mtx),
                                           num_cols((*this)(c1, c2)));
            
            Matrix cpy(mtx); // only in this case do we need to copy the argument matrix
            
            resize_block(c1, c2, maxrows, maxcols); // I think useless 
            resize(cpy, maxrows, maxcols);
            
            (*this)(c1, c2) += cpy;
        }
    } else
        insert_block(mtx, c1, c2);
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::resize_block(charge r, charge c,
                                                   size_type new_r, size_type new_c,
                                                   bool pretend)
{
    if (!pretend)
        resize((*this)(r,c), new_r, new_c);
    rows_[rows_.position(r)].second = new_r;
    cols_[cols_.position(c)].second = new_c;
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::remove_block(charge r, charge c)
{
    assert( has_block(r, c) );
    assert( rows_.position(r) == cols_.position(c) );
    
    std::size_t which = rows_.position(r);
    
    rows_.erase(rows_.begin() + which);
    cols_.erase(cols_.begin() + which);
    data_.erase(data_.begin() + which);
}
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::remove_block(std::size_t which)
{
    assert( which < data_.size() );

    rows_.erase(rows_.begin() + which);
    cols_.erase(cols_.begin() + which);
    data_.erase(data_.begin() + which);
}

#ifdef HAVE_ALPS_HDF5

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::load(alps::hdf5::archive & ar)
{
    ar >> alps::make_pvp("rows_", rows_);
    ar >> alps::make_pvp("cols_", cols_);
    std::vector<Matrix> tmp;
    ar >> alps::make_pvp("data_", tmp);
    for (typename std::vector<Matrix>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
        data_.push_back(new Matrix(*it));
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::save(alps::hdf5::archive & ar) const
{
    ar << alps::make_pvp("rows_", rows_);
    ar << alps::make_pvp("cols_", cols_);
    std::vector<Matrix> tmp(data_.begin(), data_.end());
    ar << alps::make_pvp("data_", tmp);
}

#endif

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::reserve(charge c1, charge c2,
                                              std::size_t r, std::size_t c)
{
    if (this->has_block(c1, c2))
    {
        std::size_t maxrows = std::max(rows_[rows_.position(c1)].second, r);
        std::size_t maxcols = std::max(cols_[cols_.position(c2)].second, c);
    
        rows_[rows_.position(c1)].second = maxrows;
        cols_[cols_.position(c2)].second = maxcols;
    } else {
        std::pair<charge, size_type>
        p1 = std::make_pair(c1, r),
        p2 = std::make_pair(c2, c);
        
        assert(rows_.size() == cols_.size());
        assert(rows_.size() == data_.size());
        
        size_type i1 = rows_.insert(p1);
        cols_.insert(i1, p2);
        data_.insert(data_.begin() + i1, new Matrix(1,1)); 
        /*rows_.push_back(p1);
        cols_.push_back(p2);
        data_.push_back(Matrix());*/
    }
}

template<class Matrix, class SymmGroup>
inline void block_matrix<Matrix, SymmGroup>::reserve_pos(charge c1, charge c2,
                                                         std::size_t ri, std::size_t ci)
{ reserve(c1, c2, ri+1, ci+1); }


template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::allocate_blocks()
{
    for (std::size_t k = 0; k < n_blocks(); ++k)
        resize(data_[k], rows_[k].second, cols_[k].second);
}

template<class Matrix, class SymmGroup>
bool block_matrix<Matrix, SymmGroup>::reasonable() const
{
    for (size_t k=0; k<n_blocks(); ++k)
        if (num_rows((*this)[k]) != rows_[k].second || num_cols((*this)[k]) != cols_[k].second)
            return false;
    return true;
}


