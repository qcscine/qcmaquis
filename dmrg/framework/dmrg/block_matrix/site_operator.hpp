/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#include "utils/function_objects.h"
#include "utils/bindings.hpp"

#include <boost/serialization/serialization.hpp>

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup>::SiteOperator() 
{
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup>::SiteOperator(Index<SymmGroup> const & rows,
                                              Index<SymmGroup> const & cols) : bm_(rows, cols)
{
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup>::SiteOperator(DualIndex<SymmGroup> const & basis)
: bm_(basis)
{
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup>::SiteOperator(block_matrix<Matrix,SymmGroup> const& rhs,
                                              typename SparseOperator<Matrix, SymmGroup, void>::spin_basis_type const& sb)
: spin_basis(sb), bm_(rhs), sparse_op(rhs, sb)
{
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> & SiteOperator<Matrix, SymmGroup>::operator=(SiteOperator rhs)
{
    swap(*this, rhs);
    return *this;
}

template<class Matrix, class SymmGroup>
template<class OtherMatrix>
SiteOperator<Matrix, SymmGroup> & SiteOperator<Matrix, SymmGroup>::operator=(const SiteOperator<OtherMatrix, SymmGroup> & rhs)
{
    block_matrix<Matrix, SymmGroup> cpy = rhs.bm_;
    sparse_op = rhs.sparse_op;
    swap(bm_, cpy);
    spin_ = rhs.spin();
    return *this;
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> & SiteOperator<Matrix, SymmGroup>::operator+=(SiteOperator const & rhs)
{
    assert (spin_.get() == rhs.spin().get() || n_blocks() == 0 || rhs.n_blocks() == 0);
    if (n_blocks() == 0) spin_ = rhs.spin();
    bm_ += rhs.bm_;
    return *this;
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> & SiteOperator<Matrix, SymmGroup>::operator-=(SiteOperator const & rhs)
{
    assert (spin_.get() == rhs.spin().get() || n_blocks() == 0 || rhs.n_blocks() == 0);
    if (n_blocks() == 0) spin_ = rhs.spin();
    bm_ -= rhs.bm_;
    return *this;
}

template<class Matrix, class SymmGroup>
typename SiteOperator<Matrix, SymmGroup>::size_type SiteOperator<Matrix, SymmGroup>::insert_block(Matrix const & mtx, charge c1, charge c2)
{
    return bm_.insert_block(mtx, c1, c2);  
}

template<class Matrix, class SymmGroup>
typename SiteOperator<Matrix, SymmGroup>::size_type SiteOperator<Matrix, SymmGroup>::insert_block(Matrix * mtx, charge c1, charge c2)
{
    return bm_.insert_block(mtx, c1, c2);  
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> SiteOperator<Matrix, SymmGroup>::left_basis() const
{ 
    return bm_.left_basis();
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> SiteOperator<Matrix, SymmGroup>::right_basis() const
{
    return bm_.right_basis();
}

template<class Matrix, class SymmGroup>
DualIndex<SymmGroup> const & SiteOperator<Matrix, SymmGroup>::basis() const { return bm_.basis(); }

template<class Matrix, class SymmGroup>
typename Matrix::size_type SiteOperator<Matrix, SymmGroup>::n_blocks() const { return bm_.n_blocks(); }

template<class Matrix, class SymmGroup>
std::string SiteOperator<Matrix, SymmGroup>::description() const
{
    return bm_.description();
}

template<class Matrix, class SymmGroup>
Matrix & SiteOperator<Matrix, SymmGroup>::operator[](size_type c) { return bm_[c]; }

template<class Matrix, class SymmGroup>
Matrix const & SiteOperator<Matrix, SymmGroup>::operator[](size_type c) const { return bm_[c]; }

template<class Matrix, class SymmGroup>
typename SiteOperator<Matrix, SymmGroup>::size_type SiteOperator<Matrix, SymmGroup>::find_block(charge r, charge c) const
{
    return bm_.find_block(r,c);
}

template<class Matrix, class SymmGroup>
bool SiteOperator<Matrix, SymmGroup>::has_block(charge r, charge c) const
{
    return bm_.has_block(r,c);
}

template<class Matrix, class SymmGroup>
bool SiteOperator<Matrix, SymmGroup>::has_block(std::pair<charge, size_type> const & r,
                                                std::pair<charge, size_type> const & c) const
{
    return has_block(r.first, c.first);
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type & SiteOperator<Matrix, SymmGroup>::operator()(std::pair<charge, size_type> const & r,
                                                                          std::pair<charge, size_type> const & c)
{
    return bm_(r, c);
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type const & SiteOperator<Matrix, SymmGroup>::operator()(std::pair<charge, size_type> const & r,
                                                                                std::pair<charge, size_type> const & c) const
{
    return bm_(r, c);
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> const & SiteOperator<Matrix, SymmGroup>::operator*=(const scalar_type& v)
{
    bm_ *= v;
    return *this;
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> const & SiteOperator<Matrix, SymmGroup>::operator/=(const scalar_type& v)
{
    bm_ /= v;
    return *this;
}

template<class Matrix, class SymmGroup>
typename SiteOperator<Matrix, SymmGroup>::real_type SiteOperator<Matrix, SymmGroup>::norm() const
{
    return bm_.norm();
}

template<class Matrix, class SymmGroup>
void SiteOperator<Matrix, SymmGroup>::transpose_inplace()
{
    bm_.transpose_inplace();
}

template<class Matrix, class SymmGroup>
void SiteOperator<Matrix, SymmGroup>::adjoint_inplace()
{
    bm_.adjoint_inplace();
}

template<class Matrix, class SymmGroup>
void SiteOperator<Matrix, SymmGroup>::clear()
{
    bm_.clear();
    //spin.clear();
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, SiteOperator<Matrix, SymmGroup> const & m)
{
    os << "Basis: " << m.basis() << std::endl;
    for (std::size_t k = 0; k < m.n_blocks(); ++k)
        os << "Block (" << m.basis()[k].lc << "," << m.basis()[k].rc
           << "):\n" << m[k] << std::endl;
    os << std::endl;
    return os;
}

template<class Matrix, class SymmGroup>
void SiteOperator<Matrix, SymmGroup>::match_and_add_block(Matrix const & mtx, charge c1, charge c2)
{
    bm_.match_and_add_block(mtx, c1, c2);
}

template<class Matrix, class SymmGroup>
void SiteOperator<Matrix, SymmGroup>::resize_block(charge r, charge c,
                                                   size_type new_r, size_type new_c,
                                                   bool pretend)
{
    bm_.resize_block(r, c, new_r, new_c, pretend);
}

template<class Matrix, class SymmGroup>
void SiteOperator<Matrix, SymmGroup>::resize_block(size_type pos,
                                                   size_type new_r, size_type new_c,
                                                   bool pretend)
{
    bm_.resize_block(pos, new_r, new_c, pretend);
}

template<class Matrix, class SymmGroup>
void SiteOperator<Matrix, SymmGroup>::remove_block(charge r, charge c)
{
    bm_.remove_block(r,c);
}

template<class Matrix, class SymmGroup>
void SiteOperator<Matrix, SymmGroup>::remove_block(std::size_t which)
{
    bm_.remove_block(which);
}


template<class Matrix, class SymmGroup>
template <class Archive>
void SiteOperator<Matrix, SymmGroup>::serialize(Archive & ar, const unsigned int version)
{
    ar & bm_;
}

namespace SiteOperator_detail {

    template <class Matrix, class SymmGroup>
    typename boost::disable_if<symm_traits::HasSU2<SymmGroup> >::type
    check_spin_basis(block_matrix<Matrix, SymmGroup> const & bm,
                     typename SparseOperator<Matrix, SymmGroup, void>::spin_basis_type &)
    {
    } 

    template <class Matrix, class SymmGroup>
    typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type
    check_spin_basis(block_matrix<Matrix, SymmGroup> const & bm,
                     typename SparseOperator<Matrix, SymmGroup, void>::spin_basis_type & spin_basis)
    {
        if (spin_basis.size() != bm.n_blocks())
        for(std::size_t b = 0; b < bm.n_blocks(); ++b)
            if (spin_basis.count(std::make_pair(bm.basis().left_charge(b), bm.basis().right_charge(b))) == 0)
                spin_basis[std::make_pair(bm.basis().left_charge(b), bm.basis().right_charge(b))]
                    = std::make_pair(std::vector<typename SymmGroup::subcharge>(num_rows(bm[b]), SymmGroup::spin(bm.basis().left_charge(b))),
                                     std::vector<typename SymmGroup::subcharge>(num_cols(bm[b]), SymmGroup::spin(bm.basis().right_charge(b)))
                                     );
    } 

}

template<class Matrix, class SymmGroup>
void SiteOperator<Matrix, SymmGroup>::update_sparse()
{
    SiteOperator_detail::check_spin_basis(bm_, spin_basis);
    sparse_op.update(bm_, spin_basis);
}
