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

#ifndef MAQUIS_DMRG_GSORTHO_MODIFIED_H
#define MAQUIS_DMRG_GSORTHO_MODIFIED_H

#include <cstdio>
#include "dmrg/optimize/Orthogonalizer/orthogonalizer.h"

template<class VecSpace>
class GS_ortho_mod : public Orthogonalizer<VecSpace>
{
public:
    // Types definition
    typedef Orthogonalizer<VecSpace>          base ;
    typedef typename base::magnitude_type     magnitude_type ;
    typedef typename base::scalar_type        scalar_type ;
    typedef typename base::vector_type        vector_type ;
    typedef typename base::vector_space       vector_space ;
    // Inheritance
    using base::diagonal_elements_ ;
    using base::has_additional_ ;
    using base::has_refinement ;
    using base::normalize ;
    using base::thresh_refinement_ ;
    using base::vspace_reference_ ;
    using base::vspace_additional_ ;
    using base::do_refinement_ ;
    // Constructor
    GS_ortho_mod(bool refinement = false) : base::Orthogonalizer() { do_refinement_ = refinement; } ;
    // Implementation of virtual functions
    magnitude_type get_hamiltonian(const vector_space& t, const vector_space& tA, const std::size_t& i,
                                   const std::size_t& j) ;
    void normalize(vector_type& t, vector_type& tA) ;
    void update_diagonal(vector_type& t, vector_type& tA) ;
    void orthogonalize(vector_type& t, vector_type& tA) ;
} ;

template<class VecSpace>
void GS_ortho_mod<VecSpace>::orthogonalize(vector_type& t, vector_type& tA)
{
    // Orthogonalization
    double ref = ietl::two_norm(tA) ;
    for (std::size_t jcont = 0; jcont < (*vspace_reference_).size(); jcont++) {
        t  -= (*vspace_reference_)[jcont] * ietl::dot((*vspace_additional_)[jcont],tA) ;
        tA -= (*vspace_additional_)[jcont] * ietl::dot((*vspace_additional_)[jcont],tA) ;
    }
    if (std::abs(ietl::two_norm(tA)/ref) < thresh_refinement_ && this->has_refinement()) {
        for (std::size_t jcont = 0; jcont < (*vspace_reference_).size(); jcont++) {
            t  -= (*vspace_reference_)[jcont] * ietl::dot((*vspace_additional_)[jcont],tA) ;
            tA -= (*vspace_additional_)[jcont] * ietl::dot((*vspace_additional_)[jcont],tA) ;
        }
    }
}

template<class VecSpace>
void GS_ortho_mod<VecSpace>::update_diagonal(vector_type& t, vector_type& tA)
{
    (*diagonal_elements_).push_back(1.) ;
}

template<class VecSpace>
void GS_ortho_mod<VecSpace>::normalize(vector_type& t, vector_type& tA)
{
    t  /= ietl::two_norm(tA) ;
    tA /= ietl::two_norm(tA) ;
}

template<class VecSpace>
typename GS_ortho_mod<VecSpace>::magnitude_type
         GS_ortho_mod<VecSpace>::get_hamiltonian(const vector_space& t,
                                                 const vector_space& tA,
                                                 const std::size_t& i,
                                                 const std::size_t& j)
{
    return std::real(ietl::dot(t[i], tA[j])) ;
}


#endif


