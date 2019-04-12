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

#ifndef MAQUIS_DMRG_OP_OPTIMIZER_H
#define MAQUIS_DMRG_OP_OPTIMIZER_H

#include "dmrg/optimize/MicroOptimizer/optimizationalgorithm.h"

template<class MATRIX, class VectorSpace, class CorrectionEquation>
class OP_optimizer : public OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>
{
private:
    // Types declaration
    typedef OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>   base ;
    typedef typename base::scalar_type                                       scalar_type ;
    typedef typename base::vector_type                                       vector_type ;
    typedef typename base::size_t                                            size_t ;
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
typename OP_optimizer<MATRIX, VectorSpace, CorrectionEquation>::vector_type
         OP_optimizer<MATRIX, VectorSpace, CorrectionEquation>::PerformOptimization(const vector_type& x0)
{
    vector_type input = b_ ;
    correction_.apply_precondition(input) ;
    return input ;
}

#endif

