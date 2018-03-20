/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *               2018 by Leon Freitag <lefreita@ethz.ch>
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
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef typename Matrix::value_type value_type;

    typedef Boundary<typename storage::constrained<Matrix>::type, SymmGroup>                      boundary_type ;
    typedef std::vector<boundary_type>                                                            boundary_vector ;
    typedef typename bound_database< MPS<Matrix, SymmGroup>, boundary_vector>::bound_database     bound_database ;
    typedef std::vector<boundary_vector>                                                          boundaries ;
    typedef std::vector<boundary_type*>                                                           boundary_vector_ptr ;
    typedef std::vector<boundary_vector_ptr>                                                      boundary_matrix_ptr ;
    typedef std::vector<typename contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>::schedule_t > contraction_schedule_vec;
    // Constructor with only one element
    SiteProblem(MPSTensor<Matrix, SymmGroup> const & initial,boundary_type & left_, boundary_type & right_, MPOTensor<Matrix, SymmGroup> const & mpo_)
            : mpo(mpo_)
    {
        left.push_back(&left_)   ;
        right.push_back(&right_) ;
        contraction_schedules.push_back(contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>::
                           right_contraction_schedule(initial, right, mpo));
        size = 1 ;
    }
    // Constructor with a vector of boundaries
    SiteProblem(
                std::vector<MPSTensor<Matrix, SymmGroup> >const & initial,
                std::size_t const & idx1,
                std::size_t const & idx2,
                MPOTensor<Matrix, SymmGroup> const & mpo_,
                bound_database & database)
            : mpo(mpo_)
    {
        size = 0 ;
        std::size_t ov_dim = database.n_MPS_ ;
        left.resize(ov_dim) ;
        right.resize(ov_dim) ;
        vector_coefficients.resize(ov_dim) ;
        for (std::size_t i = 0; i < database.n_MPS_; i++) {
            size_t dim = database.get_num_bound(i) ;
//             left[i].push_back(&(*(database.get_boundaries_left_sp(i,0)))[idx1]);
//             right[i].push_back(&(*(database.get_boundaries_right_sp(i,0)))[idx2]);
            //initialise the boundary with right dimension but all values set to 0. TODO: how to do it better?
            boundary_type avg_boundary = 0.0 * (*(&(*(database.get_boundaries_right_sp(0,0)))[idx2]));
            //average the boundaries
            for (std::size_t j = 0; j < dim; j++) {
                left[i].push_back(&(*(database.get_boundaries_left_sp(i,j)))[idx1]);
                right[i].push_back(&(*(database.get_boundaries_right_sp(i,j)))[idx2]);
                vector_coefficients[i].push_back(database.get_coefficients(i,j)) ;
                boundary_type add_boundary = ((scalar_type)( database.get_coefficients(i,j)*(1/dim)) ) * *((&(*(database.get_boundaries_right_sp(i,j)))[idx2]));
                avg_boundary += add_boundary;
            }
//             vector_coefficients[i].push_back(database.get_coefficients(i,0)) ;
//             boundary_type avg_boundary = *(&(*(database.get_boundaries_right_sp(0,0)))[idx2]);
            size += 1;
            contraction_schedules.push_back(contraction::Engine<Matrix, typename storage::constrained<Matrix>::type,SymmGroup>::right_contraction_schedule(initial[i], avg_boundary, mpo));

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
    contraction_schedule_vec contraction_schedules;
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
        y = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left[idx][0]), *(H.right[idx][0]), H.mpo, H.contraction_schedules[idx]) * coeff ;
        for (std::size_t k = 1; k < dimension; k++) {
            coeff = H.get_coefficient(idx,k) ;
            y += contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left[idx][k]), *(H.right[idx][k]), H.mpo, H.contraction_schedules[idx]) * coeff ;
        }
        x.make_left_paired() ;
    }
    // +----------+
    //  GET_ENERGY
    // +----------+
    template<class Matrix, class SymmGroup>
    typename maquis::traits::real_type<typename Matrix::value_type>::type get_energy(SiteProblem<Matrix, SymmGroup> const &H,
                     MPSTensor<Matrix, SymmGroup> const &x,
                     std::size_t const& idx)
    {
        MPSTensor<Matrix, SymmGroup> y = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left[idx][idx]), *(H.right[idx][idx]), H.mpo, H.contraction_schedules[idx]) ;
        typename MPS<Matrix, SymmGroup>::scalar_type res = ietl::dot(x,y) ;
        assert(check_real(res));
        typename maquis::traits::real_type<typename Matrix::value_type>::type res_dbl = maquis::real(res);
        x.make_left_paired() ;
        return res_dbl ;
    }
}

#endif
