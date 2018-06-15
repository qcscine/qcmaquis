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

#include "dmrg/optimize/Orthogonalizer/orthogonalizer.h"

// -- CONSTRUCTOR --

template<class VecSpace>
Orthogonalizer<VecSpace>::Orthogonalizer()
        : has_additional_(false),
          thresh_refinement_(0.25)
{ }

// -- SETTER --

template<class VecSpace>
void Orthogonalizer<VecSpace>::set_vecspace(vector_space& input_space_u)
{
    vspace_reference_ = &input_space_u ;
}

template<class VecSpace>
void Orthogonalizer<VecSpace>::set_addspace(vector_space& input_space_additional)
{
    vspace_additional_ = &input_space_additional ;
    has_additional_ = true ;
}

template<class VecSpace>
void Orthogonalizer<VecSpace>::set_diagonal(vector_scalar& diag_element)
{
    diagonal_elements_ = &diag_element ;
}

template<class VecSpace>
void Orthogonalizer<VecSpace>::activate_refinement(const double& thresh)
{
    thresh_refinement_ = thresh ;
    do_refinement_ = true ;
}

// -- GETTER --

template<class VecSpace>
typename Orthogonalizer<VecSpace>::vector_space*
         Orthogonalizer<VecSpace>::get_vecspace()
{
    return vspace_reference_ ;
}

template<class VecSpace>
typename Orthogonalizer<VecSpace>::vector_space*
         Orthogonalizer<VecSpace>::get_addspace()
{
    return vspace_additional_ ;
}

template<class VecSpace>
bool Orthogonalizer<VecSpace>::has_refinement()
{
    return  do_refinement_ ;
}

template<class VecSpace>
void Orthogonalizer<VecSpace>::normalize(vector_type& t)
{
    t /= ietl::two_norm(t) ;
}

