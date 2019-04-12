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

#ifndef MAQUIS_CORRECTIONEQUATION_H
#define MAQUIS_CORRECTIONEQUATION_H

#include "dmrg/optimize/singlesitevs.h"
#include "dmrg/optimize/CorrectionEquation/corrector.cpp"

template<class MATRIX, class VecSpace>
class CorrectionEquation {
public:
    // Types declaration
    typedef typename ietl::vectorspace_traits< VecSpace >::scalar_type     scalar_type ;
    typedef typename ietl::vectorspace_traits< VecSpace >::vector_type     vector_type ;
    typedef typename vector_type::bm_type                                  preconditioner_type ;
    typedef typename std::size_t                                           size_t ;
    // Methods
    CorrectionEquation() ;
    CorrectionEquation(const MATRIX& Hamiltonian, scalar_type theta, std::size_t n_root,
                       const vector_type& error, const vector_type& u, VecSpace& vs, scalar_type omega = 0.0) ;
    // Updater
    void activate_preconditioner() ;
    void update_corrector(Corrector<MATRIX, VecSpace>* pnt_crr) ;
    void update_error(const vector_type& error) ;
    void update_hamiltonian(MATRIX& hamiltonian) ;
    void update_n_root(const size_t& n_root) ;
    void update_omega(const scalar_type& omega) ;
    void update_rayleigh() ;
    void update_theta(const scalar_type& theta) ;
    void update_u(const vector_type& u) ;
    void update_Au(const vector_type& Au) ;
    void update_vecspace(const VecSpace& vs) ;
    void update_vecspace(VecSpace& vs) ;
    // Get elements
    bool do_refinement() ;
    bool do_omega() ;
    const MATRIX& get_hamiltonian() ;
    const preconditioner_type& get_preconditioner() ;
    scalar_type get_omega() ;
    scalar_type get_rayleigh() ;
    scalar_type get_theta() ;
    size_t get_n_root() ;
    vector_type get_error() ;
    vector_type get_u() ;
    vector_type get_Au() ;
    VecSpace get_vecspace() ;
    // Set correction equation
	void set_folded() ;
    void set_skew() ;
    void set_standard() ;
    void set_modified() ;
    // Methods
    vector_type apply_correction(vector_type& input) ;
    void apply_precondition(vector_type& input) ;
    void orthogonalize_simple(vector_type& input) ;
protected:
    // Methods
    void update_preconditioner() ;
    // Attributes
    bool do_precondition_, do_omega_, do_refine_error_, fully_initialized_ ;
    preconditioner_type precond_ ;
    MATRIX* Hamiltonian_ ;
    std::shared_ptr<Corrector<MATRIX, VecSpace> > corrector_ ;
    scalar_type omega_, theta_, rayleigh_ ;
    std::size_t n_root_ ;
    vector_type error_, u_, Au_ ;
    VecSpace* vs_ ;
};

#include "dmrg/optimize/CorrectionEquation/correctionequation.cpp"

#endif
