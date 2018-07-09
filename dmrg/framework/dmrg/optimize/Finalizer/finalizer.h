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

#ifndef MAQUIS_DMRG_FINALIZER_H
#define MAQUIS_DMRG_FINALIZER_H

#include <dmrg/optimize/JacobiDavidson/jacobi.h>
#include "dmrg/optimize/Finalizer/energy_computer.hpp"
#include "dmrg/optimize/Finalizer/error_computer.hpp"
#include "dmrg/optimize/singlesitevs.h"
#include "dmrg/optimize/utils/state_prop.h"

template<class MATRIX, class VecSpace>
class Finalizer {
public:
    // Types declaraton
    typedef typename ietl::vectorspace_traits<VecSpace>::real_type        real_type ;
    typedef typename ietl::vectorspace_traits<VecSpace>::scalar_type      scalar_type ;
    typedef typename ietl::vectorspace_traits<VecSpace>::vector_type      vector_type ;
    typedef typename ietl::number_traits<scalar_type>::magnitude_type     magnitude_type;
    typedef typename ietl::state_prop<VecSpace>                           state_prop ;
    typedef typename std::vector< state_prop >                            vec_prop ;
    // Constructor
    Finalizer() ;
    ~Finalizer() {} ;
    // Setter
    void set_Hamiltonian(MATRIX& ham) ;
    void set_omega(real_type omega) ;
    void set_candidate(vec_prop const& candidates) ;
    // Set the energy calculator
    void set_energy_standard() ;
    void set_energy_si_standard() ;
    void set_energy_modified() ;
    void set_energy_skew() ;
	void set_energy_folded() ;
    // Set the error calculator
    void set_error_standard() ;
    void set_error_si_standard() ;
    void set_error_modified() ;
    void set_error_skew() ;
	void set_error_folded() ;
    // Getter
    bool get_is_si() ;
    magnitude_type get_omega() ;
    magnitude_type theta_converter(magnitude_type const& theta) ;
    magnitude_type get_eigen(size_t const& idx) ;
    vector_type get_u(size_t const& idx) ;
    vector_type get_uA(size_t const& idx) ;
    MATRIX* get_Hamiltonian() ;
    // Interfaces to the call
    magnitude_type compute_energy(size_t const& i_state, size_t const& idx) ;
    vector_type compute_error(size_t const& i_state, size_t const& idx) ;
private:
    // Attributes
    MATRIX* Hamiltonian_ ;
    bool is_si_ ;
    const vec_prop* candidates_ ;
    magnitude_type omega_ ;
    EnergyComputer<MATRIX, VecSpace>* energy_computer_ ;
    ErrorComputer<MATRIX, VecSpace>* error_computer_ ;
} ;

#include "dmrg/optimize/Finalizer/finalizer.cpp"

#endif
