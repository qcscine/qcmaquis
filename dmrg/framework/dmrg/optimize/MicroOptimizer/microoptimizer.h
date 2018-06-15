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

#ifndef MAQUIS_DMRG_MICROOPTIMIZER_H
#define MAQUIS_DMRG_MICROOPTIMIZER_H

#include "dmrg/optimize/MicroOptimizer/optimizationalgorithm.h"

template<class MATRIX, class VecSpace, class CorrectionEquation>
class MicroOptimizer {
    // Types declaraton
    typedef typename ietl::vectorspace_traits< VecSpace >::scalar_type  scalar_type ;
    typedef typename ietl::vectorspace_traits< VecSpace >::vector_type  vector_type ;
public:
    // Constructor
    MicroOptimizer() = default ;
    MicroOptimizer(CorrectionEquation* correction, const float& abs_error, const float& rel_error,
                   const std::size_t& max_iter, const std::size_t& n_restart) ;
    ~MicroOptimizer() = default ;
    // Setter
    void activate_verbosity() ;
    void deactivate_verbosity() ;
    void set_opt_alg(OptimizationAlgorithm<MATRIX, VecSpace, CorrectionEquation>* opt_alg) ;
    void set_max_iter(const std::size_t& max_iter) ;
    void set_n_restart(const std::size_t& n_restart) ;
    void set_abs_tol(const float& abs_tol) ;
    void set_rel_tol(const float& rel_tol) ;
    void set_error(vector_type& error) ;
    // Getter
    CorrectionEquation get_correction() ;
    // Actual routines
    vector_type PerformOptimization(const vector_type& x0) ;
    // Attributes
protected:
    OptimizationAlgorithm<MATRIX, VecSpace, CorrectionEquation>* opt_alg_ ;
    CorrectionEquation* correction_ ;
    float abs_error_, rel_error_ ;
    std::size_t max_iter_, n_restart_ ;
};

#include "dmrg/optimize/MicroOptimizer/microoptimizer.cpp"

#endif
