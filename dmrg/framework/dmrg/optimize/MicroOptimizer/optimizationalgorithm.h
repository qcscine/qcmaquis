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

#ifndef MAQUIS_DMRG_OPTIMIZATIONALGORITHM_H
#define MAQUIS_DMRG_OPTIMIZATIONALGORITHM_H

// Forward declaration

template<class MATRIX, class VectorSpace, class CorrectionEquation>
class MicroOptimizer ;

// +------------------------------+
//  GENERAL OPTIMIZATION ALGORITHM
// +------------------------------+

template<class MATRIX, class VectorSpace, class CorrectionEquation>
class OptimizationAlgorithm
{
protected:
    // Types declaration
    typedef typename ietl::vectorspace_traits< VectorSpace >::scalar_type           scalar_type ;
    typedef typename ietl::vectorspace_traits< VectorSpace >::vector_type           vector_type ;
    typedef typename std::size_t                                                    size_t ;
    typedef typename boost::numeric::ublas::matrix< scalar_type >                   matrix_scalar ;
    typedef typename std::vector< scalar_type >                                     vector_scalar ;
    typedef typename std::vector< vector_type >                                     vector_space ;
public:
    // Constructors
    OptimizationAlgorithm() : verbosity_(false) {} ;
    OptimizationAlgorithm(const float& abs_tol, const float& rel_tol, const size_t& max_iter) ;
    OptimizationAlgorithm(const float& abs_tol, const float& rel_tol, const size_t& max_iter, const vector_type& r) ;
    virtual ~OptimizationAlgorithm() {} ;
    // Set
    void activate_verbosity() ;
    void deactivate_verbosity() ;
    void set_max_iter(const size_t& max_iter) ;
    void set_n_restart(const size_t& n_restart) ;
    void set_abs_tol(const float& abs_tol) ;
    void set_rel_tol(const float& rel_tol) ;
    void set_error(vector_type& error) ;
    // Virtual function, specific for each optimization
    virtual vector_type perform_optimization(MicroOptimizer<MATRIX, VectorSpace, CorrectionEquation>* optimizer,
                                             const vector_type& x0) = 0 ;
protected:
    // Attributes
    bool verbosity_ ;
    vector_type b_ ;
    float abs_tol_, rel_tol_ ;
    std::size_t max_iter_, n_restart_ ;
} ;

#include "dmrg/optimize/MicroOptimizer/optimizationalgorithm.cpp"
#include "dmrg/optimize/MicroOptimizer/GMRES_optimizer.hpp"
#include "dmrg/optimize/MicroOptimizer/conjugategradient_optimizer.hpp"
#include "dmrg/optimize/MicroOptimizer/biconjugategradient_optimizer.hpp"
#include "dmrg/optimize/MicroOptimizer/onlyprec_optimizer.hpp"

#endif
