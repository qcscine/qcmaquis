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

#include "dmrg/optimize/utils/bound_database.h"

#ifndef SITE_PROBLEM_H
#define SITE_PROBLEM_H

template<class Matrix, class SymmGroup>
struct SiteProblem
{
    //TODO ALB the constructor with only one boundary is kept only for back-compatibility
    // Types definition
    typedef Boundary<typename storage::constrained<Matrix>::type, SymmGroup>                      boundary_type ;
    typedef std::vector<boundary_type>                                                            boundary_vector ;
    typedef typename bound_database< MPS<Matrix, SymmGroup>, boundary_vector>::bound_database     bound_database ;
    typedef std::vector<boundary_vector>                                                          boundaries ;
    typedef std::vector<boundary_type*>                                                           boundary_vector_ptr ;
    typedef std::vector<boundary_vector_ptr>                                                      boundary_matrix_ptr ;
    // Constructor with only one element
    SiteProblem(boundary_type & left_, boundary_type & right_, MPOTensor<Matrix, SymmGroup> const & mpo_)
            : mpo(mpo_)
    {
        left.push_back(&left_)   ;
        right.push_back(&right_) ;
        size = 1 ;
    }
    // Constructor with a vector of boundaries
    SiteProblem(MPOTensor<Matrix, SymmGroup> const & mpo_,
                std::size_t const & idx1,
                std::size_t const & idx2,
                bound_database & database)
            : mpo(mpo_) {
        size = 0 ;
        std::size_t ov_dim = database.n_MPS_ ;
        left.resize(ov_dim) ;
        right.resize(ov_dim) ;
        vector_coefficients.resize(ov_dim) ;
        for (std::size_t i = 0; i < database.n_MPS_; i++) {
            size_t dim = database.get_num_bound(i) ;
            for (std::size_t j = 0; j < dim; j++) {
                left[i].push_back(&(*(database.get_boundaries_left_sp(i,j)))[idx1]);
                right[i].push_back(&(*(database.get_boundaries_right_sp(i,j)))[idx2]);
                vector_coefficients[i].push_back(database.get_coefficients(i,j)) ;
            }
            size += 1;
        }
    }
    // Methods
    float get_coefficient(const std::size_t& idx, const std::size_t& k) const
    {
        return vector_coefficients[idx][k] ;
    }
    //
    inline size_t get_size(void) const
    {
        return this->size ;
    }
public:
    // Attributes (public)
    boundary_matrix_ptr left ;
    boundary_matrix_ptr right ;
    std::size_t size ;
    MPOTensor<Matrix, SymmGroup> const & mpo ;
private:
    std::vector< std::vector <float> > vector_coefficients ;
};

namespace ietl {
    // +----+
    //  MULT
    // +----+
    // Overloading of the ietl::mult function, to mimic a matrix-vector multiplication
    template<class Matrix, class SymmGroup>
    void mult(SiteProblem<Matrix, SymmGroup> const &H,
              MPSTensor <Matrix, SymmGroup> const &x,
              MPSTensor <Matrix, SymmGroup> &y,
              std::size_t const& idx = 0)
    {
        assert( idx < H.get_size() && H.left[idx].size() == H.right[idx].size() ) ;
        std::size_t dimension = H.left[idx].size() ;
        float coeff = H.get_coefficient(idx,0) ;
        y = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left[idx][0]), *(H.right[idx][0]), H.mpo) * coeff ;
        for (std::size_t k = 1; k < dimension; k++) {
            coeff = H.get_coefficient(idx,k) ;
            y += contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left[idx][k]), *(H.right[idx][k]), H.mpo) * coeff ;
        }
        x.make_left_paired() ;
    }
    // +----------+
    //  GET_ENERGY
    // +----------+
    template<class Matrix, class SymmGroup>
    float get_energy(SiteProblem<Matrix, SymmGroup> const &H,
                     MPSTensor<Matrix, SymmGroup> const &x,
                     std::size_t const& idx)
    {
        MPSTensor<Matrix, SymmGroup> y = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left[idx][idx]), *(H.right[idx][idx]), H.mpo) ;
        float res = ietl::dot(x,y) ;
        x.make_left_paired() ;
        return res ;
    }
}

#endif
