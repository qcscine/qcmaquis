/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *
 * This software is part of the ALPS libraries, published under the ALPS
 * Library License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Library License along with
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

#include <boost/numeric/bindings/blas/detail/blas_names.h>
#include <boost/numeric/bindings/lapack/detail/lapack_names.h>
#include <boost/numeric/bindings/lapack/detail/lapack.h>
#include <iostream>
#include <cmath>
#include <numeric>
#include <cassert>

#include "varianceoptimizer.h"

// -- Default constructor --
template<class SymmGroup>
VarianceOptimizer<SymmGroup>::VarianceOptimizer() : H_(0,0), H_squared_(0,0), Hessian_(0,0), Hessian_inverse_(0,0),
                                         n_elements_(0), coefficients_(), coefficients_current_(), scaling_(1),
                                         lagrange_multiplier_(0.), penalty_(10.)
{} ;

// -- Useful constructor --
template<class SymmGroup>
VarianceOptimizer<SymmGroup>::VarianceOptimizer(matrix_type& H, matrix_type& H_squared) : H_(0,0), H_squared_(0,0), Hessian_(0,0),
                                                                               Hessian_inverse_(0,0), coefficients_(),
                                                                               coefficients_current_(), scaling_(1),
                                                                               lagrange_multiplier_(0.), penalty_(10.)
{
    // Simple checks
    assert (H.nrows() == H_squared.nrows() && H.ncols() == H_squared.ncols() && H.nrows() == H.ncols()) ;
    n_elements_ = H.ncols() ;
    n_elements_squared_ = n_elements_*n_elements_ ;
    // Scaling
    scaling_ = H(n_elements_-1, n_elements_-1) ;
    // Copy the data
    H_.resize(n_elements_, n_elements_) ;
    H_squared_.resize(n_elements_, n_elements_) ;
    std::copy(H.data(), H.data()+n_elements_squared_, H_.data()) ;
    std::copy(H_squared.data(), H_squared.data()+n_elements_squared_, H_squared_.data()) ;
    // Debug print
    if (verbose_) {
        std::cout << " -- Initial variance -- " << std::endl ;
        for (std::size_t idx = 0; idx < n_elements_; idx++)
            std::cout << " State " << idx << " variance " << sqrt(H_squared(idx, idx)-std::pow(H(idx, idx),2.))
                      << std::endl ;
    }
    std::for_each(H_.data(), H_.data()+n_elements_squared_, [&](double& a) { a /= scaling_ ; } ) ;
    std::for_each(H_squared_.data(), H_squared_.data()+n_elements_squared_, [&](double& a) { a /= (scaling_*scaling_) ; } ) ;
    // Prepares the Hessian matrices
    gradient_.resize(n_elements_, 0.) ;
    Hessian_.resize(n_elements_, n_elements_) ;
    Hessian_inverse_.resize(n_elements_, n_elements_) ;
    // Coefficients initialization
    coefficients_.resize(n_elements_, 0.) ;
    coefficients_[0] = 1. ;
    coefficients_current_.resize(n_elements_, 0.) ;
};

// -- Distance between two vectors --
template<class SymmGroup>
double VarianceOptimizer<SymmGroup>::compute_distance(const vector_ref& vec1, const vector_ref& vec2)
{
    return std::inner_product(vec1.begin(), vec1.end(), vec2.begin(), 0., std::plus<double>(),
                              [] (double a, double b) {return pow(a-b,2) ;} ) ;
}

// -- Apply the Hamiltonian to a vector --
template<class SymmGroup>
typename VarianceOptimizer<SymmGroup>::vector_ref
         VarianceOptimizer<SymmGroup>::apply_operator(const vector_ref& input,
                                           const bool& is_squared,
                                           const bool& is_inverted)
{
    // Initialization
    double* matrix_ptr  ;
    vector_ref output(n_elements_, 0.) ;
    if (is_squared)
        matrix_ptr = H_squared_.data() ;
    else if (is_inverted)
        matrix_ptr = Hessian_inverse_.data() ;
    else
        matrix_ptr = H_.data() ;
    //
    for (std::size_t idx = 0; idx < n_elements_; idx++)
        output[idx] = std::accumulate(input.begin(), input.end(), 0.,
                                      [&] (double a, double b) { return a+b*(*(matrix_ptr++)) ; } ) ;
    return output ;
}

// -- Get the coefficients --
template<class SymmGroup>
std::vector<double> VarianceOptimizer<SymmGroup>::get_coeff()
{
    return coefficients_ ;
}

// -- Performs the optimization --
template<class SymmGroup>
void VarianceOptimizer<SymmGroup>::PerformOptimization()
{
    // Initialization
    bool converged_outer = false ;
    std::size_t n_penalty_iter = 0 ;
    vector_ref step(n_elements_, 0.) ;
    if (verbose_)
        std::cout << std::endl << " -- Variance optimization -- " << std::endl ;
    // Actual optimization
    double dist ;
    while (n_penalty_iter <= max_iter_outer_ && !converged_outer) {
        // Reinitialization of the variable
        bool converged = false ;
        std::size_t n_iter = 0 ;
        std::copy(coefficients_.begin(), coefficients_.end(), coefficients_current_.begin()) ;
        // Inner iteration
        while (!converged && n_iter <= max_iter_inner_) {
            // Computes the gradient
            double norm_gradient = compute_derivatives() ;
            invert_Hessian() ;
            if (norm_gradient < thresh_inner_) {
                converged = true;
            } else {
                // Update of the coefficients with the Newton method
                std::fill_n(step.begin(), n_elements_, 0.) ;
                step = apply_operator(gradient_, false, true) ;
                vector_ref::iterator it = step.begin() ;
                std::for_each(coefficients_.begin(), coefficients_.end(), [&](double& a) { a -= *(it++) ; } ) ;
                ++n_iter;
            }
        }
        if (verbose_)
            std::cout << " Scaling parameter - " << penalty_  << " convergence in " << n_iter
                      << " iterations. Gradient norm " << compute_derivatives() << std::endl  ;
        dist = compute_distance(coefficients_, coefficients_current_) ;
        if (dist < thresh_outer_)
            converged_outer = true ;
        if (add_lagrange_)
            lagrange_multiplier_ -= penalty_*(compute_norm(true)-1.) ;
        ++n_penalty_iter ;
        penalty_ *= scaling_penalty_ ;
    }
    if (verbose_) {
        std::cout << " Overall convergence reached with penalty factor " << penalty_
                  << " difference " << dist << std::endl ;
        for (std::size_t idx = 0; idx < n_elements_; idx++)
            std::cout << " Coefficients number " << idx << " = " << coefficients_[idx] << std::endl;
    }
}

// -- Derivatives calculation --
template<class SymmGroup>
double VarianceOptimizer<SymmGroup>::compute_derivatives()
{
    assert (gradient_.size() == n_elements_) ;
    // Computes the variance of the current approximation and some preliminary quantities
    double norm = compute_norm(true), H_value = 0. , H_squared_value = 0. ;
    for (std::size_t idx1 = 0; idx1 < n_elements_; idx1++) {
        for (std::size_t idx2 = 0; idx2 < n_elements_; idx2++) {
            H_value += coefficients_[idx1]*coefficients_[idx2]*H_(idx1, idx2) ;
            H_squared_value += coefficients_[idx1]*coefficients_[idx2]*H_squared_(idx1, idx2) ;
        }
    }
    // Loop over the indexes of the gradient
    for (std::size_t idx = 0; idx < n_elements_; idx++) {
        // Loop over the "fake" indexes
        double jnk = 0. , jnk2 = 0. ;
        gradient_[idx] = 0. ;
        //
        for (std::size_t i = 0; i < n_elements_; i++) {
            jnk  += coefficients_[i] * H_(i, idx) ;
            jnk2 += coefficients_[i] * H_squared_(i, idx) ;
        }
        gradient_[idx]  = 2.*jnk2 - 4.*H_value*jnk ;
        gradient_[idx] += (2.*penalty_*(norm-1.) - lagrange_multiplier_)*2.*coefficients_[idx] ;
        // Hessian calculation
        for (std::size_t idx1 = 0; idx1 < n_elements_; idx1++) {
            double jnk3 = 0. ;
            for (std::size_t i = 0; i < n_elements_; i++)
                jnk3 += coefficients_[i] * H_(i, idx1) ;
            Hessian_(idx, idx1) = 2.*H_squared_(idx, idx1) - 4.*H_value*H_(idx, idx1) - 8.*jnk*jnk3
                                + 8.*penalty_*coefficients_[idx]*coefficients_[idx1] ;
            if (idx == idx1) {
                Hessian_(idx, idx1) += 4. * penalty_ * (norm - 1.) ;
                if (add_lagrange_)
                    Hessian_(idx, idx1) -= 2.*lagrange_multiplier_ ;
            }
        }
    }
    // Calculation of the norm of the gradient
    double norm_gradient = compute_norm(false) ;
    return sqrt(norm_gradient) ;
}

// -- Computes norm of relevant vectors --
template<class SymmGroup>
double VarianceOptimizer<SymmGroup>::compute_norm(const bool& is_coeff)
{
    vector_iterator it_begin, it_end ;
    if (is_coeff) {
        it_begin = coefficients_.begin() ;
        it_end   = coefficients_.end() ;
    } else {
        it_begin = gradient_.begin() ;
        it_end   = gradient_.end() ;
    }
    double res = std::inner_product(it_begin, it_end, it_begin, 0.) ;
    return res ;
}

// -- Invert the Hessian matrix --
template<class SymmGroup>
void VarianceOptimizer<SymmGroup>::invert_Hessian()
{
    // Set up dimensions
    fortran_int_t m = n_elements_, n = n_elements_ ;
    fortran_int_t dim_iwork = n_elements_ , dim_work = n_elements_squared_ ;
    std::copy(Hessian_.data(), Hessian_.data()+n_elements_squared_, Hessian_inverse_.data()) ;
    fortran_int_t *ipiv = new fortran_int_t[n_elements_] ;
    auto *work        = new double[n_elements_squared_] ;
    fortran_int_t info ;
    // LU factorization
    LAPACK_DGETRF(&m, &n, Hessian_inverse_.data(), &dim_iwork, ipiv, &info) ;
    // Inversion of the matrix
    LAPACK_DGETRI(&m, Hessian_inverse_.data(), &m, ipiv, work, &dim_work, &info) ;
    // Frees the memory
    delete [] ipiv ;
    delete [] work ;
}

