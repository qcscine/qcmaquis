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

#ifndef MAQUIS_DMRG_CG_OPTIMIZER_H
#define MAQUIS_DMRG_CG_OPTIMIZER_H

#include "dmrg/optimize/MicroOptimizer/optimizationalgorithm.h"

template<class MATRIX, class VectorSpace, class CorrectionEquation>
class CG_optimizer : public OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>
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
typename CG_optimizer<MATRIX, VectorSpace, CorrectionEquation>::vector_type
         CG_optimizer<MATRIX, VectorSpace, CorrectionEquation>::PerformOptimization(const vector_type& x0)
{
    // Types definiton
    vector_type x = x0 ;     //  Approximation of the solution at the current iteration
    vector_type p , Ap;      //  Conjugate vector at the current iteration, and the results of the application of A to p
    vector_type r, r_back ;  //  Error vector
    vector_type z ;          //  Precondition vector
    vector_type jnk ;        //  Scratch vector
    scalar_type alpha ;      //  Coefficient to use to update x
    scalar_type beta ;       //  Coefficient to use to update r
    scalar_type aerr ;       //  Absolute error
    scalar_type rerr ;       //  Relative error
    scalar_type sprod ;
    double normb = ietl::two_norm(b_) , normr ;
    // Initialization (0-th iteration)
    vector_type init = x0 ;
    jnk = correction_.apply_correction(init) ;
    r = b_ - jnk ;
    z = r ;
    correction_.apply_precondition(z) ;
    p = z ;
    Ap = correction_.apply_correction(p) ;
    for (std::size_t i=0 ; i < max_iter_ ; ++i)
    {
        r_back = r ;
        sprod = ietl::dot(r, z) ;
        alpha = sprod / ietl::dot(p, Ap) ;
        x += alpha*p ;
        r -= alpha*Ap ;
        normr = ietl::two_norm(r) ;
        if (verbosity_)
            std::cout << "Conj. Grad. iteration " << i << ", Abs. Err. = " << ietl::two_norm(r)
                                                       << ", Rel. Err. = " << ietl::two_norm(r)/normb
                                                       << std::endl ;
        if (normr/normb < rel_tol_ || normr < abs_tol_) {
            break ;
        } else {
            z = r ;
            correction_.apply_precondition(z) ;
            beta = ietl::dot(z,(r-r_back)) / sprod ;
            jnk = beta*p ;
            p   = z + jnk ;
            Ap  = correction_.apply_correction(p) ;
        }
    }
    return x ;
}

#endif

