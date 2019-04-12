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

#include "dmrg/optimize/Finalizer/finalizer.h"

// Constructor

template<class MATRIX, class VecSpace>
Finalizer<MATRIX, VecSpace>::Finalizer()
    : is_si_(false)
{}

template<class MATRIX, class VecSpace>
void Finalizer<MATRIX, VecSpace>::set_candidate(vec_prop const& candidates)
{
    candidates_ = &candidates ;
}

template<class MATRIX, class VecSpace>
void Finalizer<MATRIX, VecSpace>::set_omega(real_type omega)
{
    omega_ = omega ;
}

template<class MATRIX, class VecSpace>
void Finalizer<MATRIX, VecSpace>::set_Hamiltonian(MATRIX &ham)
{
    Hamiltonian_ = &ham ;
}

// Getter

template<class MATRIX, class VecSpace>
bool Finalizer<MATRIX, VecSpace>::get_is_si()
{
    return is_si_ ;
}
template<class MATRIX, class VecSpace>
typename Finalizer<MATRIX, VecSpace>::magnitude_type
         Finalizer<MATRIX, VecSpace>::get_eigen(size_t idx)
{
    return (*candidates_)[idx].theta_ ;
}

template<class MATRIX, class VecSpace>
typename Finalizer<MATRIX, VecSpace>::magnitude_type
         Finalizer<MATRIX, VecSpace>::get_omega()
{
    return omega_ ;
}

template<class MATRIX, class VecSpace>
typename Finalizer<MATRIX, VecSpace>::vector_type
         Finalizer<MATRIX, VecSpace>::get_u(size_t idx)
{
    return ((*candidates_)[idx].u_) ;
}

template<class MATRIX, class VecSpace>
typename Finalizer<MATRIX, VecSpace>::vector_type
         Finalizer<MATRIX, VecSpace>::get_uA(size_t idx)
{
    return ((*candidates_)[idx].uA_) ;
}

template<class MATRIX, class VecSpace>
MATRIX* Finalizer<MATRIX, VecSpace>::get_Hamiltonian()
{
    return Hamiltonian_ ;
}
