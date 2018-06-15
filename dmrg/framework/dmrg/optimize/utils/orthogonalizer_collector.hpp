/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#include <ietl/vectorspace.h>
#include "orthogonalizer_collector.h"

// -- STANDARD CONSTRUCTOR --

template<class MATRIX>
orthogonalizer_collector<MATRIX>::orthogonalizer_collector() :
    orthovec_constant(), orthovec_otherstate(), orthovec_withinstate(),
    dim_const(0), dim_otherstate(0), dim_overall(0), dim_withinstate(0)
{
    orthovec_constant.resize(0) ;
    orthovec_otherstate.resize(0) ;
    orthovec_withinstate.resize(0) ;
} ;

// -- CONSTRUCTOR WITH CONSTANT VECTORS GIVEN AS INPUT --

template<class MATRIX>
orthogonalizer_collector<MATRIX>::orthogonalizer_collector(vecspace &const_vec)
        : orthovec_withinstate(),
          orthovec_otherstate(),
          orthovec_constant(const_vec),
          dim_const(0), dim_otherstate(0), dim_overall(0), dim_withinstate(0)
{
    // Set dimensions
    dim_const       = orthovec_constant.size() ;
    dim_overall     = dim_const ;
}

// -- COPY CONSTRUCTOR --

template<class MATRIX>
orthogonalizer_collector<MATRIX>::orthogonalizer_collector(orthogonalizer_collector& rhs)
{
    // Dimensions
    dim_const        = rhs.dim_const ;
    dim_overall      = rhs.dim_overall ;
    dim_withinstate  = rhs.dim_withinstate ;
    dim_otherstate   = rhs.dim_otherstate ;
    // Vector spaces
    orthovec_withinstate = rhs.orthovec_withinstate ;
    orthovec_constant    = rhs.orthovec_constant ;
    orthovec_otherstate  = rhs.orthovec_otherstate ;
}

// -- CONSTRUCTOR WITH VECTOR GIVEN AS INPUT --

template<class MATRIX>
orthogonalizer_collector<MATRIX>::orthogonalizer_collector(vecspace &const_vec,
                                                           vecspace &within_vec,
                                                           vecspace &other_vec)
    : orthovec_withinstate(within_vec),
      orthovec_otherstate(other_vec),
      orthovec_constant(const_vec),
      dim_const(0), dim_otherstate(0), dim_overall(0), dim_withinstate(0)
{
    // Set dimensions
    dim_const       = orthovec_constant.size() ;
    dim_otherstate  = orthovec_otherstate.size() ;
    dim_withinstate = orthovec_withinstate.size() ;
    dim_overall     = dim_const + dim_otherstate + dim_withinstate ;
}

// -- METHODS --

template <class MATRIX>
void orthogonalizer_collector<MATRIX>::clear_all()
{
    // Reset dimensions
    dim_const       = 0 ;
    dim_otherstate  = 0 ;
    dim_withinstate = 0 ;
    dim_overall     = 0 ;
    // Reset vector spaces
    orthovec_withinstate.clear() ;
    orthovec_otherstate.clear() ;
    orthovec_constant.clear() ;
}

template<class MATRIX>
void orthogonalizer_collector<MATRIX>::clear_local()
{
    // Reset dimensions
    dim_overall     = dim_const + dim_otherstate ;
    dim_withinstate = 0 ;
    // Reset vector spaces
    orthovec_withinstate.clear() ;
}

template<class MATRIX>
void orthogonalizer_collector<MATRIX>::add_other_vec(const MATRIX &to_add)
{
    MATRIX jnk = to_add ;
    this->orthogonalize(jnk) ;
    bool add = this->normalize(jnk) ;
    if (add) {
        orthovec_otherstate.push_back(to_add);
        dim_otherstate++;
    }
}

template<class MATRIX>
void orthogonalizer_collector<MATRIX>::add_within_vec(const MATRIX &to_add)
{
    MATRIX jnk = to_add ;
    this->orthogonalize(jnk) ;
    bool add = this->normalize(jnk) ;
    if (add) {
        orthovec_withinstate.push_back(jnk);
        dim_withinstate++;
    }
}

// -- ORTHOGONALIZATION ROUTINE --

template<class MATRIX>
void orthogonalizer_collector<MATRIX>::orthogonalize(MATRIX &to_ortho) const
{
    // Initialization
    ref_iterator it ;
    double ref = ietl::two_norm(to_ortho) ;
    // Loop over the constant vector
    for (it = orthovec_constant.begin(); it != orthovec_constant.end(); it++)
        to_ortho -= *it * ietl::dot(*it, to_ortho) / ietl::dot(*it, *it) ;
    if (ietl::two_norm(to_ortho)/ref < refinement_threshold)
        for (it = orthovec_constant.begin(); it != orthovec_constant.end(); it++)
            to_ortho -= *it * ietl::dot(*it, to_ortho) / ietl::dot(*it, *it);
    // Loop over the within vector
    for (it = orthovec_withinstate.begin(); it != orthovec_withinstate.end(); it++)
        to_ortho -= *it * ietl::dot(*it, to_ortho) / ietl::dot(*it, *it) ;
    if (ietl::two_norm(to_ortho)/ref < refinement_threshold)
        for (it = orthovec_withinstate.begin(); it != orthovec_withinstate.end(); it++)
            to_ortho -= *it * ietl::dot(*it, to_ortho) / ietl::dot(*it, *it) ;
    // Loop over the within vector
    for (it = orthovec_otherstate.begin(); it != orthovec_otherstate.end(); it++)
        to_ortho -= *it * ietl::dot(*it, to_ortho) / ietl::dot(*it, *it) ;
    if (ietl::two_norm(to_ortho)/ref < refinement_threshold)
        for (it = orthovec_otherstate.begin(); it != orthovec_otherstate.end(); it++)
            to_ortho -= *it * ietl::dot(*it, to_ortho) / ietl::dot(*it, *it) ;
}

// -= NORMALIZATION ROUTINE --

template<class MATRIX>
bool orthogonalizer_collector<MATRIX>::normalize(MATRIX& to_norm)
{
    bool ret = true ;
    if (ietl::two_norm(to_norm) > this->normalize_threshold) {
        to_norm /= ietl::two_norm(to_norm);
    } else {
        ret = false ;
    }
    return ret ;
}
