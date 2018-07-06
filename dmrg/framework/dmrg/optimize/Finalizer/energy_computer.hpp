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

#ifndef MAQUIS_ENERGY_COMPUTER_H
#define MAQUIS_ENERGY_COMPUTER_H

// Forward declaration

#include "finalizer.h"

template<class MATRIX, class VecSpace>
class Finalizer ;

// +------------------------+
//  General corrector object
// +------------------------+

template<class MATRIX, class VecSpace>
class EnergyComputer {
public:
    // Types definition
    typedef Finalizer<MATRIX, VecSpace>       base ;
    typedef typename base::real_type          real_type ;
    typedef typename base::scalar_type        scalar_type ;
    typedef typename base::vector_type        vector_type ;
    // Constructor and desctructor
    virtual ~EnergyComputer() {} ;
    // Actual method
    virtual scalar_type compute_energy(Finalizer<MATRIX, VecSpace>* finalizer,
                                       size_t const& i_state,
                                       size_t const& idx ) = 0 ;
    virtual real_type theta_converter(real_type const& theta) = 0  ;
} ;

#include "dmrg/optimize/Finalizer/standard_energy_computer.hpp"
#include "dmrg/optimize/Finalizer/standard_si_energy_computer.hpp"
#include "dmrg/optimize/Finalizer/modified_energy_computer.hpp"
#include "dmrg/optimize/Finalizer/skew_energy_computer.hpp"
#include "dmrg/optimize/Finalizer/folded_energy_computer.hpp"

#endif
