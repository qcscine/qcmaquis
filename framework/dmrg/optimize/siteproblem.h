/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *
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

// +----------------------+
//  Site problem structure
// +----------------------+
//
// This structure contains left and right boundary, together with the
// MPO tensor of the site of interest

#ifndef SITE_PROBLEM_H
#define SITE_PROBLEM_H

template<class Matrix, class SymmGroup>
struct SiteProblem
{
    //TODO ALB the constructor with only one boundary is kept only for back-compatibility
    //TODO ALB should be generalized
    // Types definition
    typedef Boundary<typename storage::constrained<Matrix>::type, SymmGroup> boundary_type ;
    typedef std::vector<boundary_type> boundary_vector ;
    typedef std::vector<boundary_vector> boundaries ;
    typedef std::vector<boundary_type*> boundary_vector_ptr ;
    // Constructor with only one element
    SiteProblem(boundary_type & left_, boundary_type & right_, MPOTensor<Matrix, SymmGroup> const & mpo_)
            : mpo(mpo_)
    {
        left.push_back(&left_)   ;
        right.push_back(&right_) ;
        size = 1 ;
    }
    // Constructor with a vector of boundaries
    SiteProblem(boundaries & left_vec_ ,
                boundaries & right_vec_ ,
                MPOTensor<Matrix, SymmGroup> const & mpo_,
                std::size_t const & idx1,
                std::size_t const & idx2)
            : mpo(mpo_)
    {
        size = 0 ;
        assert (left_vec_.size() == right_vec_.size()) ;
        for (typename boundaries::iterator it1 = left_vec_.begin(), it2 = right_vec_.begin() ; it1 != left_vec_.end(), it2 != right_vec_.end() ; it1++, it2++) {
            left.push_back(&((*it1)[idx1]));
            right.push_back(&((*it2)[idx2]));
            size += 1 ;
        }
    }
    // Attributes (public)
    boundary_vector_ptr left;
    boundary_vector_ptr right;
    std::size_t size ;
    MPOTensor<Matrix, SymmGroup> const & mpo;
};

// Overloading of the ietl::mult function, to mimic a matrix-vector multiplication

namespace ietl {
    template<class Matrix, class SymmGroup>
    void mult(SiteProblem<Matrix, SymmGroup> const &H,
              MPSTensor <Matrix, SymmGroup> const &x,
              MPSTensor <Matrix, SymmGroup> &y,
              std::size_t const& idx = 0 ) {
        assert( idx < H.left.size() && idx < H.right.size() ) ;
        y = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left[idx]), *(H.right[idx]), H.mpo);
        x.make_left_paired();
    }
}

#endif
