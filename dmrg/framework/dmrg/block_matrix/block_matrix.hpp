/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "dmrg/block_matrix/block_matrix.h"

#include "utils/function_objects.h"
#include "utils/bindings.hpp"

#include <alps/type_traits/is_complex.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/ptr_container/serialize_ptr_vector.hpp>

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>::block_matrix() 
{
}

// final
template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>::block_matrix(Index<SymmGroup> rows,
                                              Index<SymmGroup> cols)
: rows_(rows)
, cols_(cols)
{
    assert(rows_.size() == cols_.size());

    basis_.resize(rows_.size());
    for (size_type k = 0; k < rows_.size(); ++k)
        basis_[k] = boost::make_tuple(rows_[k].first, cols_[k].first, rows_[k].second, cols_[k].second);

    for (size_type k = 0; k < rows_.size(); ++k)
        data_.push_back(new Matrix(rows_[k].second, cols_[k].second));
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>::block_matrix(block_matrix const& rhs)
: rows_(rhs.left_basis())
, cols_(rhs.right_basis())
, basis_(rhs.basis())
, data_(rhs.data_)
{
    #ifdef AMBIENT_TRACKING
    if(!rhs.label.empty()){
        this->label = rhs.label + "'";
        ambient_track_as(*this, this->label);
    }
    #endif
}

template<class Matrix, class SymmGroup>
template <class OtherMatrix>
block_matrix<Matrix, SymmGroup>::block_matrix(block_matrix<OtherMatrix,SymmGroup> const& rhs)
: rows_(rhs.left_basis())
, cols_(rhs.right_basis())
, basis_(rhs.basis())
{
    data_.reserve(rhs.n_blocks());
    for (size_type k = 0; k < rhs.n_blocks(); ++k)
        data_.push_back(new Matrix(rhs[k]));
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & block_matrix<Matrix, SymmGroup>::operator=(block_matrix rhs)
{
    swap(*this, rhs);
    return *this;
}

// final
template<class Matrix, class SymmGroup>
template<class OtherMatrix>
block_matrix<Matrix, SymmGroup> & block_matrix<Matrix, SymmGroup>::operator=(const block_matrix<OtherMatrix, SymmGroup> & rhs)
{
    rows_ = rhs.rows_;
    cols_ = rhs.cols_;
    basis_ = rhs.basis_;
    data_.resize(rhs.data_.size());
    for(int k = 0; k < data_.size(); k++){
        data_[k].resize(num_rows(rhs.data_[k]), num_cols(rhs.data_[k]));
        data_[k] = rhs.data_[k];
    }
    return *this;
}

// TODO: final
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

// TODO: final
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

// TODO: final
template<class Matrix, class SymmGroup>
typename block_matrix<Matrix, SymmGroup>::size_type block_matrix<Matrix, SymmGroup>::insert_block(Matrix const & mtx, charge c1, charge c2)
{
    assert( !has_block(c1, c2) );
    
    std::pair<charge, size_type>
    p1 = std::make_pair(c1, num_rows(mtx)),
    p2 = std::make_pair(c2, num_cols(mtx));
    
    size_type i1 = rows_.insert(p1);
    cols_.insert(i1, p2);
    basis_.insert(boost::make_tuple(c1, c2, num_rows(mtx), num_cols(mtx)));
    Matrix* block = new Matrix(mtx);
    data_.insert(data_.begin() + i1, block);
#ifdef AMBIENT_TRACKING
    ambient_track_as(*block, this->label);
#endif
    
    return i1;
}

// TODO: final
template<class Matrix, class SymmGroup>
typename block_matrix<Matrix, SymmGroup>::size_type block_matrix<Matrix, SymmGroup>::insert_block(Matrix * mtx, charge c1, charge c2)
{
    assert( !has_block(c1, c2) );
    
    std::pair<charge, size_type>
    p1 = std::make_pair(c1, num_rows(*mtx)),
    p2 = std::make_pair(c2, num_cols(*mtx));
    
    size_type i1 = rows_.insert(p1);
    cols_.insert(i1, p2);
    basis_.insert(boost::make_tuple(c1, c2, num_rows(*mtx), num_cols(*mtx)));
    data_.insert(data_.begin() + i1, mtx);
#ifdef AMBIENT_TRACKING
    ambient_track_as(*mtx, this->label);
#endif
    
    return i1;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> const & block_matrix<Matrix, SymmGroup>::left_basis() const { return rows_; }

template<class Matrix, class SymmGroup>
Index<SymmGroup> const & block_matrix<Matrix, SymmGroup>::right_basis() const { return cols_; }

template<class Matrix, class SymmGroup>
DualIndex<SymmGroup> const & block_matrix<Matrix, SymmGroup>::basis() const { return basis_; }

template<class Matrix, class SymmGroup>
typename Matrix::size_type block_matrix<Matrix, SymmGroup>::n_blocks() const { return data_.size(); }

template<class Matrix, class SymmGroup>
std::string block_matrix<Matrix, SymmGroup>::description() const
{
    std::ostringstream oss;
    oss << rows_ << cols_;
    return oss.str();
}

// final
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::shift_basis(block_matrix<Matrix, SymmGroup>::charge diff)
{
    rows_.shift(diff);
    cols_.shift(diff);
    basis_.shift(diff);
}

template<class Matrix, class SymmGroup>
Matrix & block_matrix<Matrix, SymmGroup>::operator[](size_type c) { return data_[c]; }

template<class Matrix, class SymmGroup>
Matrix const & block_matrix<Matrix, SymmGroup>::operator[](size_type c) const { return data_[c]; }

// TODO: internal
template<class Matrix, class SymmGroup>
typename block_matrix<Matrix, SymmGroup>::size_type block_matrix<Matrix, SymmGroup>::find_block(charge r, charge c) const
{
    std::size_t p1 = rows_.position(r);
    std::size_t p2 = cols_.position(c);
    
    if (p1 == p2 && p1 != rows_.size()) return p1;
    else                                return n_blocks();
}

// TODO: internal
template<class Matrix, class SymmGroup>
bool block_matrix<Matrix, SymmGroup>::has_block(charge r, charge c) const
{
    return basis_.position(r,c) != basis_.size();
}

// TODO: internal
template<class Matrix, class SymmGroup>
bool block_matrix<Matrix, SymmGroup>::has_block(std::pair<charge, size_type> const & r,
                                                std::pair<charge, size_type> const & c) const
{
    return has_block(r.first, c.first);
}

// TODO: internal
template<class Matrix, class SymmGroup>
typename Matrix::value_type & block_matrix<Matrix, SymmGroup>::operator()(std::pair<charge, size_type> const & r,
                                                                          std::pair<charge, size_type> const & c)
{
    assert( rows_.position(r.first) == cols_.position(c.first) );
    return data_[rows_.position(r.first)](r.second, c.second);
}

// TODO: internal
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
    // todo: check if "omp for" used in nested regions
    for(size_t k = 0; k < n_blocks(); ++k){
        select_proc(ambient::scope::balance(k,n_blocks()));
        data_[k] *= v;
    }
    return *this;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & block_matrix<Matrix, SymmGroup>::operator/=(const scalar_type& v)
{
    // todo: check if "omp for" used in nested regions
    for(size_t k = 0; k < n_blocks(); ++k){
        select_proc(ambient::scope::balance(k,n_blocks()));
        data_[k] /= v;
    }
    return *this;
}


template<class Matrix, class SymmGroup>
typename block_matrix<Matrix, SymmGroup>::scalar_type block_matrix<Matrix, SymmGroup>::trace() const
{
    std::vector<scalar_type> vt; vt.reserve(data_.size());
    std::transform(data_.begin(), data_.end(), back_inserter(vt), utils::functor_trace());
    return maquis::accumulate(vt.begin(), vt.end(), scalar_type(0.));
}

template<class Matrix, class SymmGroup>
typename block_matrix<Matrix, SymmGroup>::real_type block_matrix<Matrix, SymmGroup>::norm() const
{
    std::vector<real_type> vt; vt.reserve(data_.size());
    for(size_t k = 0; k < n_blocks(); ++k){
        select_proc(ambient::scope::balance(k,n_blocks()));
        vt.push_back(norm_square(data_[k]));
    }
    return maquis::sqrt(maquis::accumulate(vt.begin(), vt.end(), real_type(0.)));
}

// TODO: final
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::transpose_inplace()
{
    std::for_each(data_.begin(), data_.end(), utils::functor_transpose_inplace());
    std::swap(rows_, cols_);
    for (std::size_t i=0; i < basis_.size(); ++i) {
        std::swap(boost::tuples::get<0>(basis_[i]), boost::tuples::get<1>(basis_[i]));
        std::swap(boost::tuples::get<2>(basis_[i]), boost::tuples::get<3>(basis_[i]));
    }
}

template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::conjugate_inplace()
{
    std::for_each(data_.begin(), data_.end(), utils::functor_conj_inplace());
}

// TODO: final
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::adjoint_inplace()
{
    std::for_each(data_.begin(), data_.end(), utils::functor_adjoint_inplace());
    std::swap(rows_, cols_);
    for (std::size_t i=0; i < basis_.size(); ++i) {
        std::swap(boost::tuples::get<0>(basis_[i]), boost::tuples::get<1>(basis_[i]));
        std::swap(boost::tuples::get<2>(basis_[i]), boost::tuples::get<3>(basis_[i]));
    }
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

// TODO: final
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::clear()
{
    data_.clear();
    rows_ = Index<SymmGroup>();
    cols_ = Index<SymmGroup>();
    basis_ = DualIndex<SymmGroup>();
}

// final
template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, block_matrix<Matrix, SymmGroup> const & m)
{
    os << "Left HS: " << m.left_basis() << std::endl;
    os << "Right HS: " << m.right_basis() << std::endl;
    for (std::size_t k = 0; k < m.n_blocks(); ++k)
        os << "Block (" << m.left_basis()[k].first << "," << m.right_basis()[k].first << "):\n" << m[k] << std::endl;
    os << std::endl;
    return os;
}

// internal
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

// final
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::resize_block(charge r, charge c,
                                                   size_type new_r, size_type new_c,
                                                   bool pretend)
{
    if (!pretend)
        resize((*this)(r,c), new_r, new_c);
    rows_[rows_.position(r)].second = new_r;
    cols_[cols_.position(c)].second = new_c;

    std::size_t pos = basis_.position(r,c);
    boost::tuples::get<2>(basis_[pos]) = new_r;
    boost::tuples::get<3>(basis_[pos]) = new_c;
}

// final
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::remove_block(charge r, charge c)
{
    assert( has_block(r, c) );
    assert( rows_.position(r) == cols_.position(c) );
    
    std::size_t which = rows_.position(r);
    std::size_t which2 = basis_.position(r,c);
    
    rows_.erase(rows_.begin() + which);
    cols_.erase(cols_.begin() + which);
    basis_.erase(basis_.begin() + which2);
    data_.erase(data_.begin() + which);
}

// final
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::remove_block(std::size_t which)
{
    assert( which < data_.size() );

    rows_.erase(rows_.begin() + which);
    cols_.erase(cols_.begin() + which);
    basis_.erase(basis_.begin() + which);
    data_.erase(data_.begin() + which);
}

// final
template<class Matrix, class SymmGroup>
template<class Archive>
void block_matrix<Matrix, SymmGroup>::load(Archive & ar)
{
    ar["rows_"] >> rows_;
    ar["cols_"] >> cols_;
    basis_.resize(rows_.size());
    for (std::size_t s = 0; s < rows_.size(); ++s)
        basis_[s] = boost::make_tuple(rows_[s].first, cols_[s].first, rows_[s].second, cols_[s].second);

    data_.clear();
    if (alps::is_complex<typename Matrix::value_type>() && !ar.is_complex("data_"))
    {
        #ifdef USE_AMBIENT
        printf("ERROR: NOT TESTED!\n\n");
        #endif
        typedef typename alps::numeric::matrix<typename alps::numeric::real_type<typename Matrix::value_type>::type> LoadMatrix;
        std::vector<LoadMatrix> tmp;
        ar["data_"] >> tmp;
        for(typename std::vector<LoadMatrix>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
            data_.push_back(new Matrix(maquis::bindings::matrix_cast<Matrix>(*it)));
    } else {
        std::vector<Matrix> tmp;
        ar["data_"] >> tmp;
        // TODO: is swap here possible?
        for(typename std::vector<Matrix>::const_iterator it = tmp.begin(); it != tmp.end(); ++it)
            data_.push_back(new Matrix(*it));
    }
}

// final
template<class Matrix, class SymmGroup>
template<class Archive>
void block_matrix<Matrix, SymmGroup>::save(Archive & ar) const
{
    ar["rows_"] << rows_;
    ar["cols_"] << cols_;

    std::vector<Matrix> tmp(data_.begin(), data_.end());
    ar["data_"] << tmp;
}

// final
template<class Matrix, class SymmGroup>
template <class Archive>
void block_matrix<Matrix, SymmGroup>::serialize(Archive & ar, const unsigned int version)
{
    ar & rows_ & cols_ & data_;
}


// internal
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::reserve(charge c1, charge c2,
                                              std::size_t r, std::size_t c)
{
    if (this->has_block(c1, c2))
    {
        assert (rows_.position(c1) == cols_.position(c2));
        std::size_t pos = rows_.position(c1);
        std::size_t maxrows = std::max(rows_[pos].second, r);
        std::size_t maxcols = std::max(cols_[pos].second, c);
    
        rows_[pos].second = maxrows;
        cols_[pos].second = maxcols;
        boost::tuples::get<2>(basis_[pos]) = maxrows;
        boost::tuples::get<3>(basis_[pos]) = maxcols;
    } else {
        std::pair<charge, size_type>
        p1 = std::make_pair(c1, r),
        p2 = std::make_pair(c2, c);
        
        assert(rows_.size() == cols_.size());
        assert(rows_.size() == data_.size());
        
        size_type i1 = rows_.insert(p1);
        cols_.insert(i1, p2);
        basis_.insert(boost::make_tuple(c1, c2, r, c));
        Matrix* block = new Matrix(1,1);
        data_.insert(data_.begin() + i1, block); 
#ifdef AMBIENT_TRACKING
        ambient_track_as(*block, this->label);
#endif
    }
    assert( this->has_block(c1,c2) );
}

template<class Matrix, class SymmGroup>
inline void block_matrix<Matrix, SymmGroup>::reserve_pos(charge c1, charge c2,
                                                         std::size_t ri, std::size_t ci)
{ reserve(c1, c2, ri+1, ci+1); }


// internal
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::allocate_blocks()
{
    for (std::size_t k = 0; k < n_blocks(); ++k)
        resize(data_[k], rows_[k].second, cols_[k].second);
}

// internal
template<class Matrix, class SymmGroup>
bool block_matrix<Matrix, SymmGroup>::reasonable() const
{
    for (size_t k=0; k<n_blocks(); ++k)
        if (num_rows((*this)[k]) != rows_[k].second || num_cols((*this)[k]) != cols_[k].second)
            return false;
    return true;
}

template<class Matrix, class SymmGroup>
std::size_t block_matrix<Matrix, SymmGroup>::num_elements() const
{
    size_t ret = 0;
    for (size_t k = 0; k < n_blocks(); ++k)
        ret += num_rows(data_[k])*num_cols(data_[k]);
    return ret;
}

#ifdef USE_AMBIENT
template<class Matrix, class SymmGroup>
void block_matrix<Matrix, SymmGroup>::print_distribution() const
{
    if(!ambient::master()) return;
    double total = num_elements();
    printf("%.2f GB:", total*sizeof(typename Matrix::value_type)/1024/1024/1024);
    for(int p = 0; p < ambient::num_procs(); ++p){
        double part = 0;
        for(int i = 0; i < this->n_blocks(); ++i){
            if(!ambient::weak((*this)[i][0]) && ambient::get_owner((*this)[i][0]) == p)
                part += num_rows((*this)[i])*num_cols((*this)[i]);
        }
        printf(" %.1f%%", 100*part/total);
    }
    printf("\n");
}
#endif
