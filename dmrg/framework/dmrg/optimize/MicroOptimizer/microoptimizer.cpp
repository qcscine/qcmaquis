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

#include "dmrg/optimize/MicroOptimizer/microoptimizer.h"

// -- CONSTRUCTOR --


template<class MATRIX, class VecSpace, class CorrectionEquation>
MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::MicroOptimizer(CorrectionEquation* correction,
                                                                     const float &abs_error,
                                                                     const float &rel_error,
                                                                     const std::size_t& max_iter,
                                                                     const std::size_t& n_restart)
        : correction_(correction),
          abs_error_(abs_error),
          rel_error_(rel_error),
          max_iter_(max_iter),
          n_restart_(n_restart)
{
    set_opt_alg(new GMRES_optimizer<MATRIX, VecSpace, CorrectionEquation>()) ;
    set_abs_tol(abs_error) ;
    set_rel_tol(rel_error) ;
    set_max_iter(max_iter) ;
    set_n_restart(n_restart) ;
}

// -- SET --

template<class MATRIX, class VecSpace, class CorrectionEquation>
void MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::activate_verbosity()
{
    opt_alg_->activate_verbosity() ;
}

template<class MATRIX, class VecSpace, class CorrectionEquation>
void MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::deactivate_verbosity()
{
    opt_alg_->deactivate_verbosity() ;
}

template<class MATRIX, class VecSpace, class CorrectionEquation>
void MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::set_opt_alg
        (OptimizationAlgorithm<MATRIX, VecSpace, CorrectionEquation>* opt_alg)
{
    opt_alg_ = opt_alg ;
    opt_alg_->set_abs_tol(abs_error_) ;
    opt_alg_->set_rel_tol(rel_error_) ;
    opt_alg_->set_max_iter(max_iter_) ;
    opt_alg_->set_n_restart(n_restart_) ;
}

template<class MATRIX, class VecSpace, class CorrectionEquation>
void MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::set_abs_tol(const float& abs_tol)
{
    opt_alg_->set_abs_tol(abs_tol) ;
}

template<class MATRIX, class VecSpace, class CorrectionEquation>
void MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::set_rel_tol(const float& rel_tol)
{
    opt_alg_->set_rel_tol(rel_tol) ;
}

template<class MATRIX, class VecSpace, class CorrectionEquation>
void MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::set_max_iter(const std::size_t &max_iter)
{
    opt_alg_->set_max_iter(max_iter) ;
}


template<class MATRIX, class VecSpace, class CorrectionEquation>
void MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::set_n_restart(const std::size_t &n_restart)
{
    opt_alg_->set_n_restart(n_restart) ;
}

template<class MATRIX, class VecSpace, class CorrectionEquation>
void MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::set_error(vector_type &error)
{
    vector_type jnk ;
    if (correction_->do_refinement())
        jnk = error*correction_->get_rayleigh() ;
    else
        jnk = -error ;
    opt_alg_->set_error(jnk) ;
}

// -- GET --

template<class MATRIX, class VecSpace, class CorrectionEquation>
CorrectionEquation MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::get_correction()
{
    return *correction_ ;
}

template<class MATRIX, class VecSpace, class CorrectionEquation>
typename MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::vector_type
         MicroOptimizer<MATRIX, VecSpace, CorrectionEquation>::PerformOptimization(const vector_type &x0)
{
    vector_type result ;
    result = opt_alg_->perform_optimization(this, x0) ;
    return result ;
};

