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
    SiteProblem(MPOTensor<Matrix, SymmGroup> const & mpo_,
                MPOTensor<Matrix, SymmGroup> const & mpo_squared_,
                std::size_t const & idx1,
                std::size_t const & idx2,
                bound_database & database,
                const bool& add_squared)
            : mpo(mpo_), mpo_squared(mpo_squared_), is_squared_(add_squared)
    {
        size = 0 ;
        std::size_t ov_dim = database.n_MPS_ ;
        left.resize(ov_dim) ;
        left_squared.resize(ov_dim) ;
        right.resize(ov_dim) ;
        right_squared.resize(ov_dim) ;
        vector_coefficients.resize(ov_dim) ;
        for (std::size_t i = 0; i < database.n_MPS_; i++) {
            size_t dim = database.get_num_bound(i) ;
            boundary_type avg_boundary = 0.0 * (*(&(*(database.get_boundaries_right_sp(0,0, false)))[idx2]));
            for (std::size_t j = 0; j < dim; j++) {
                // Normal operator
                left[i].push_back(&(*(database.get_boundaries_left_sp(i,j, false)))[idx1]);
                right[i].push_back(&(*(database.get_boundaries_right_sp(i,j, false)))[idx2]);
                // Squared operator
                if (is_squared_) {
                    left_squared[i].push_back(&(*(database.get_boundaries_left_sp(i,j, true)))[idx1]);
                    right_squared[i].push_back(&(*(database.get_boundaries_right_sp(i,j, true)))[idx2]);
                }
                // Coefficients
                vector_coefficients[i].push_back(database.get_coefficients(i,j)) ;
                boundary_type add_boundary = ((scalar_type)( database.get_coefficients(i,j)*(1/dim)) ) * *((&(*(database.get_boundaries_right_sp(i,j, false)))[idx2]));
                avg_boundary += add_boundary;
            }
            size += 1;
            contraction_schedules.push_back(contraction::Engine<Matrix, typename storage::constrained<Matrix>::type,SymmGroup>::right_contraction_schedule(*(database.get_mps(i,idx1)), avg_boundary, mpo));
        }
    }
    // Methods
    double get_coefficient(const std::size_t& idx,
                           const std::size_t& k) const
    {
        return vector_coefficients[idx][k] ;
    }
    //
    inline size_t get_size() const
    {
        return this->size ;
    }
    //
    inline bool is_squared() const
    {
        return this->is_squared_ ;
    }
public:
    // Attributes (public)
    boundary_matrix_ptr left, left_squared ;
    boundary_matrix_ptr right, right_squared ;
    std::size_t size ;
    MPOTensor<Matrix, SymmGroup> const & mpo, mpo_squared ;
    contraction_schedule_vec contraction_schedules;
private:
    std::vector< std::vector <double> > vector_coefficients ;
    bool is_squared_ ;
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
              std::size_t const& idx = 0,
              bool const& do_squared = false)
    {
        assert( idx < H.get_size() && H.left[idx].size() == H.right[idx].size() ) ;
        std::size_t dimension = H.left[idx].size() ;
        double coeff = H.get_coefficient(idx,0) ;
        if (do_squared) {
            y = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left_squared[idx][0]),
                                                                               *(H.right_squared[idx][0]),
                                                                                 H.mpo_squared,
                                                                                 H.contraction_schedules[idx]) * coeff ;
        } else{
            y = contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left[idx][0]),
                                                                               *(H.right[idx][0]),
                                                                                 H.mpo,
                                                                                 H.contraction_schedules[idx]) * coeff ;
        }
        for (std::size_t k = 1; k < dimension; k++) {
            coeff = H.get_coefficient(idx,k) ;
            if (do_squared)
                if (H.is_squared())
                    y += contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left_squared[idx][k]),
                                                                                        *(H.right_squared[idx][k]),
                                                                                          H.mpo_squared,
                                                                                          H.contraction_schedules[idx]) * coeff ;
                else
                    throw std::runtime_error("Squared boundaries not available") ;
            else
                y += contraction::Engine<Matrix, Matrix, SymmGroup>::site_hamil2(x, *(H.left[idx][k]),
                                                                                    *(H.right[idx][k]),
                                                                                      H.mpo,
                                                                                      H.contraction_schedules[idx]) * coeff ;
        }
        x.make_left_paired() ;
    }
    // +----------+
    //  GET_ENERGY
    // +----------+
    template<class Matrix, class SymmGroup>
    typename maquis::traits::real_type<typename Matrix::value_type>::type get_energy(SiteProblem<Matrix, SymmGroup> const &H,
                     MPSTensor<Matrix, SymmGroup> const &x,
                     std::size_t const& idx,
                     const bool& do_squared)
    {
        MPSTensor<Matrix, SymmGroup> y ;
        ietl::mult(H, x, y, idx, do_squared) ;
        typename maquis::traits::real_type<typename Matrix::value_type>::type res_dbl = maquis::real(ietl::dot(x,y));
        x.make_left_paired() ;
        return res_dbl ;
    }
}

#endif
