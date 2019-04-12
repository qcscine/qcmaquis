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
OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::OptimizationAlgorithm(
                                                                                      CorrectionEquation& correction,
                                                                                      double abs_tol,
                                                                                      double rel_tol,
                                                                                      size_t max_iter,
                                                                                      size_t n_restart,
                                                                                      const vector_type &r)
    : correction_(correction),
      b_(correction_.do_refinement() ? r*correction_.get_rayleigh() : -r),
      abs_tol_(abs_tol) , rel_tol_(rel_tol) , max_iter_(max_iter) , n_restart_(n_restart), verbosity_(false) {};

template<class MATRIX, class VectorSpace, class CorrectionEquation>
OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::OptimizationAlgorithm(
                                                                                      CorrectionEquation& correction,
                                                                                      double abs_tol,
                                                                                      double rel_tol,
                                                                                      size_t max_iter,
                                                                                      size_t n_restart
                                                                                      )
    : correction_(correction),  b_(0), abs_tol_(abs_tol) , rel_tol_(rel_tol) , max_iter_(max_iter), n_restart_(n_restart), verbosity_(false) {};


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
void OptimizationAlgorithm<MATRIX, VectorSpace, CorrectionEquation>::set_error(const vector_type &error)
{
    b_ = correction_.do_refinement() ? error*correction_.get_rayleigh() : -error;
}


