/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef MAQUIS_DMRG_BICGS_OPTIMIZER_H
#define MAQUIS_DMRG_BICGS_OPTIMIZER_H

#include "dmrg/optimize/MicroOptimizer/optimizationalgorithm.h"

template<class MATRIX, class VectorSpace, class CorrectionEquation>
class BICGS_optimizer : public OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>
{
private:
    // Types declaration
    typedef OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>   base ;
    typedef typename base::scalar_type                                       scalar_type ;
    typedef typename base::vector_type                                       vector_type ;
    typedef typename base::size_t                                            size_t ;
    typedef typename base::matrix_scalar                                     matrix_scalar ;
    typedef typename base::vector_scalar                                     vector_scalar ;
    typedef typename base::vector_space                                      vector_space ;
    // Attributes
    using base::abs_tol_ ;
    using base::activate_verbosity ;
    using base::b_ ;
    using base::max_iter_ ;
    using base::rel_tol_ ;
    using base::verbosity_ ;
    using base::correction_;
    using base::base;
public:
    // Overriding of the perform optimization algorithm
    virtual vector_type PerformOptimization(const vector_type& x0) ;
};

// Routine performing the optimization
template<class MATRIX, class VectorSpace, class CorrectionEquation>
typename BICGS_optimizer<MATRIX, VectorSpace, CorrectionEquation>::vector_type
         BICGS_optimizer<MATRIX, VectorSpace, CorrectionEquation>::PerformOptimization(const vector_type& x0)
{
    // Types definiton
    vector_type x = x0 ;        //  Approximation of the solution at the current iteration
    vector_type p , Ap , pinv ; //  Conjugate vector at the current iteration, and the results of the application of A to p
    vector_type s , As , sinv , Asinv ;
    vector_type r, r0 ;         //  Error vectors
    vector_type jnk ;           //  Scratch vector
    scalar_type alpha = 1. ;    //  Coefficient to use to update x
    scalar_type beta  ;
    scalar_type omega = 1. ;
    double normb = ietl::two_norm(b_) ;
    scalar_type rho_pr = 1. ;
    scalar_type rho_act ;
    // Initialization (0-th iteration)
    jnk = correction_.apply_correction(x) ;
    r0 = b_ - jnk ;
    r  = r0 ;
    p  = 0.*r ;
    Ap = 0.*r ;
    for (std::size_t i = 0 ; i < max_iter_ ; ++i)
    {
        rho_act = ietl::dot(r0, r) ;
        beta = (rho_act*alpha)/(rho_pr*omega) ;
        p    = r + beta*(p - omega*Ap) ;
        pinv = p ;
        correction_.apply_precondition(pinv) ;
        Ap      = correction_.apply_correction(pinv) ;
        alpha   = rho_act/ietl::dot(Ap, r0) ;
        s       = r - alpha*Ap ;
        if ( verbosity_ )
            std::cout << "Conj. Grad. iteration " << i << ", Abs. Err. = " << ietl::two_norm(s)
                      << ", Rel. Err. = " << ietl::two_norm(s)/normb
                      << std::endl ;
        if ( ietl::two_norm(s)/normb < rel_tol_ || ietl::two_norm(s) < abs_tol_ ) {
            x += alpha*pinv ;
            return x ;
        }
        sinv    = s ;
        correction_.apply_precondition(sinv) ;
        As      = correction_.apply_correction(sinv) ;
        Asinv   = As ;
        correction_.apply_precondition(Asinv) ;
        omega   = ietl::dot(Asinv,sinv) / ietl::dot(Asinv, Asinv) ;
        x      += alpha*pinv + omega*sinv ;
        r       = s - omega * As ;
        if ( verbosity_ )
            std::cout << "Conj. Grad. iteration " << i << ", Abs. Err. = " << ietl::two_norm(r)
                      << ", Rel. Err. = " << ietl::two_norm(r)/normb
                      << std::endl ;
        if ( ietl::two_norm(r)/normb < rel_tol_ || ietl::two_norm(r) < abs_tol_ )
            return x ;
        rho_pr  = rho_act ;
        rho_act = ietl::dot(r0, r) ;
        if (std::abs(rho_act) < 1.0E-6) {
            r0 = r ;
            p  = r ;
        }
    }
    return x ;
}

#endif

