/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#include <dmrg/mp_tensors/contractions/abelian/h_diagonal.hpp>
#include "dmrg/optimize/CorrectionEquation/correctionequation.h"

// Class Constructors and Destructors

template<class MATRIX, class VecSpace>
CorrectionEquation<MATRIX, VecSpace>::CorrectionEquation() :
        do_omega_(false),
        do_precondition_(false),
        do_refine_error_(false),
        fully_initialized_(false)
{
    this->set_standard() ;
}

template<class MATRIX, class VecSpace>
CorrectionEquation<MATRIX, VecSpace>::CorrectionEquation(const MATRIX &Hamiltonian,
                                                         const scalar_type &theta,
                                                         const std::size_t &n_root,
                                                         const vector_type &error,
                                                         const vector_type &u,
                                                         VecSpace &vs) :
        Hamiltonian_(&Hamiltonian),
        theta_(theta),
        n_root_(n_root),
        error_(error),
        u_(u),
        vs_(&vs),
        do_omega_(false),
        do_precondition_(false),
        do_refine_error_(false),
        fully_initialized_(false)
{
    corrector_ = new StandardCorrection<MATRIX, VecSpace>() ;

}

template<class MATRIX, class VecSpace>
CorrectionEquation<MATRIX, VecSpace>::CorrectionEquation(const MATRIX &Hamiltonian,
                                                         const scalar_type &theta,
                                                         const std::size_t& n_root,
                                                         const vector_type &error,
                                                         const vector_type &u,
                                                         VecSpace &vs,
                                                         const scalar_type &omega) :
        Hamiltonian_(&Hamiltonian),
        theta_(theta),
        n_root_(n_root),
        error_(error),
        u_(u),
        vs_(&vs),
        do_omega_(true),
        do_precondition_(false),
        do_refine_error_(false),
        fully_initialized_(false),
        omega_(omega)
{
    corrector_ = new StandardCorrection<MATRIX, VecSpace>() ;
}

// Update the corrector algorithm
template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_corrector(Corrector<MATRIX, VecSpace> *pnt_crr)
{
    corrector_ = pnt_crr ;
}

// Get the Hamiltonian

template<class MATRIX, class VecSpace>
bool CorrectionEquation<MATRIX, VecSpace>::do_omega()
{
    return do_omega_ ;
}

template<class MATRIX, class VecSpace>
bool CorrectionEquation<MATRIX, VecSpace>::do_refinement()
{
    return do_refine_error_ ;
}
template<class MATRIX, class VecSpace>
MATRIX CorrectionEquation<MATRIX, VecSpace>::get_hamiltonian()
{
    return *Hamiltonian_ ;
}

template<class MATRIX, class VecSpace>
typename CorrectionEquation<MATRIX, VecSpace>::preconditioner_type
         CorrectionEquation<MATRIX, VecSpace>::get_preconditioner()
{
    return precond_ ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_hamiltonian(MATRIX &hamiltonian)
{
    Hamiltonian_ = &hamiltonian ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_preconditioner()
{
    if (fully_initialized_) {
        assert (Hamiltonian_->left[n_root_].size() == Hamiltonian_->right[n_root_].size());
        std::size_t dim = Hamiltonian_->left[n_root_].size();
        precond_ = contraction::diagonal_hamiltonian(*(Hamiltonian_->left[n_root_][0]),
                                                     *(Hamiltonian_->right[n_root_][0]),
                                                       Hamiltonian_->mpo, u_) *
                                                       Hamiltonian_->get_coefficient(n_root_, 0);
        for (std::size_t idx = 1; idx < dim; idx++)
            precond_ += contraction::diagonal_hamiltonian(*(Hamiltonian_->left[n_root_][idx]),
                                                          *(Hamiltonian_->right[n_root_][idx]),
                                                            Hamiltonian_->mpo, u_) *
                                                            Hamiltonian_->get_coefficient(n_root_, idx);
    } else {
        throw std::runtime_error("Number of root not initialized") ;
    }
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_error(const vector_type &error)
{
    error_ = error ;
}

template<class MATRIX, class VecSpace>
typename CorrectionEquation<MATRIX, VecSpace>::vector_type CorrectionEquation<MATRIX, VecSpace>::get_error()
{
    return error_ ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_omega(const scalar_type &omega)
{
    omega_ = omega ;
}

template<class MATRIX, class VecSpace>
typename CorrectionEquation<MATRIX, VecSpace>::scalar_type CorrectionEquation<MATRIX, VecSpace>::get_omega()
{
    return omega_ ;
}

template<class MATRIX, class VecSpace>
typename CorrectionEquation<MATRIX, VecSpace>::scalar_type CorrectionEquation<MATRIX, VecSpace>::get_rayleigh()
{
    return rayleigh_ ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_theta(const scalar_type &theta)
{
    theta_ = theta ;
}

template<class MATRIX, class VecSpace>
typename CorrectionEquation<MATRIX, VecSpace>::scalar_type CorrectionEquation<MATRIX, VecSpace>::get_theta()
{
    return theta_ ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_rayleigh()
{
    rayleigh_ = ietl::dot(u_,Au_) / ietl::dot(u_,u_) ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_u(const vector_type &u)
{
    u_ = u ;
}

template<class MATRIX, class VecSpace>
typename CorrectionEquation<MATRIX, VecSpace>::vector_type CorrectionEquation<MATRIX, VecSpace>::get_u()
{
    return u_ ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_Au(const vector_type &Au)
{
    Au_ = Au ;
}

template<class MATRIX, class VecSpace>
typename CorrectionEquation<MATRIX, VecSpace>::vector_type CorrectionEquation<MATRIX, VecSpace>::get_Au()
{
    return Au_ ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_vecspace(const VecSpace& vs)
{
    vs_ = &vs ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_vecspace(VecSpace& vs)
{
    vs_ = &vs ;
}

template<class MATRIX, class VecSpace>
VecSpace CorrectionEquation<MATRIX, VecSpace>::get_vecspace()
{
    return *vs_ ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::update_n_root(const size_t &n_root)
{
    n_root_ = n_root ;
    fully_initialized_ = true ;
    update_preconditioner() ;
}

template<class MATRIX, class VecSpace>
std::size_t CorrectionEquation<MATRIX, VecSpace>::get_n_root()
{
    return n_root_ ;
}

// -- Set different correction equation algorithm --

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::set_standard()
{
    corrector_ = new StandardCorrection<MATRIX, VecSpace>() ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::set_folded()
{
    corrector_ = new FoldedCorrection<MATRIX, VecSpace>() ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::set_skew()
{
    corrector_ = new SkewCorrection<MATRIX, VecSpace>() ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::set_modified()
{
    corrector_ = new ModifiedCorrection<MATRIX, VecSpace>() ;
    do_refine_error_ = true ;
}

// -- Method used to activate preconditioning --

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::activate_preconditioner()
{
    do_precondition_ = true ;
}


// -- Method used to activate omega --

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::activate_omega()
{
    do_omega_ = true ;
}

// -- Actual methods --

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::orthogonalize_simple(vector_type &input)
{
    vs_->project(input) ;
}

template<class MATRIX, class VecSpace>
typename CorrectionEquation<MATRIX, VecSpace>::vector_type
         CorrectionEquation<MATRIX, VecSpace>::apply_correction(vector_type &input)
{
    vector_type result ;
    result = corrector_->apply(this, input) ;
    return result ;
}

template<class MATRIX, class VecSpace>
void CorrectionEquation<MATRIX, VecSpace>::apply_precondition(vector_type &input)
{
    if (do_precondition_)
        corrector_->precondition(this, input) ;
}
