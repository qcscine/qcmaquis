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

#include "dmrg/optimize/MicroOptimizer/optimizationalgorithm.h"

// -- Constructor --

template<class MATRIX, class VectorSpace, class CorrectionEquation>
OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::OptimizationAlgorithm(const float &abs_tol,
                                                                                      const float &rel_tol,
                                                                                      const size_t &max_iter,
                                                                                      const vector_type &r)
    : abs_tol_(abs_tol) , rel_tol_(rel_tol) , max_iter_(max_iter) , b_(r), verbosity_(false)
{ } ;

template<class MATRIX, class VectorSpace, class CorrectionEquation>
OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::OptimizationAlgorithm(const float &abs_tol,
                                                                                      const float &rel_tol,
                                                                                      const size_t &max_iter)
        : abs_tol_(abs_tol) , rel_tol_(rel_tol) , max_iter_(max_iter), verbosity_(false)
{
    b_ = 0 ;
} ;

// -- Set --

template<class MATRIX, class VectorSpace, class CorrectionEquation>
void OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::activate_verbosity()
{
    verbosity_ = true ;
}

template<class MATRIX, class VectorSpace, class CorrectionEquation>
void OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::deactivate_verbosity()
{
    verbosity_ = false ;
}

template<class MATRIX, class VectorSpace, class CorrectionEquation>
void OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::set_abs_tol(const float &abs_tol)
{
    abs_tol_ = abs_tol ;
}

template<class MATRIX, class VectorSpace, class CorrectionEquation>
void OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::set_rel_tol(const float &rel_tol)
{
    rel_tol_ = rel_tol ;

}

template<class MATRIX, class VectorSpace, class CorrectionEquation>
void OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::set_max_iter(const size_t &max_iter)
{
    max_iter_ = max_iter ;
}

template<class MATRIX, class VectorSpace, class CorrectionEquation>
void OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::set_n_restart(const size_t &n_restart)
{
    n_restart_ = n_restart ;
}

template<class MATRIX, class VectorSpace, class CorrectionEquation>
void OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::set_error(vector_type &error)
{
    b_ = error ;
}


