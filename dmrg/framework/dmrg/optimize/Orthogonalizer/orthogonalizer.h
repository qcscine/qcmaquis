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

#ifndef MAQUIS_DMRG_ORTHOGONALIZER_H
#define MAQUIS_DMRG_ORTHOGONALIZER_H

template<class VecSpace>
class Orthogonalizer {
protected:
    // Types declaraton
    typedef typename ietl::vectorspace_traits< VecSpace >::scalar_type    scalar_type ;
    typedef typename ietl::vectorspace_traits< VecSpace >::vector_type    vector_type ;
    typedef typename ietl::number_traits<scalar_type>::magnitude_type     magnitude_type;
    typedef typename std::vector< vector_type >                           vector_space ;
    typedef typename std::vector< scalar_type >                           vector_scalar ;
public:
    // Constructor
    Orthogonalizer() ;
    ~Orthogonalizer() {} ;
    // Setter
    void activate_refinement(const double& thresh) ;
    void set_addspace(vector_space& input_space_additional) ;
    void set_diagonal(vector_scalar& diag_element) ;
    void set_vecspace(vector_space& input_space_u) ;
    // Getter
    bool has_refinement() ;
    std::vector<vector_type>* get_vecspace() ;
    std::vector<vector_type>* get_addspace() ;
    // Standard orthonormalizer
    void normalize(vector_type& t) ;
    // Actual routines
    virtual magnitude_type get_hamiltonian(const vector_space& t, const vector_space& tA, const std::size_t& i,
                                           const std::size_t& j) = 0 ;
    virtual void normalize(vector_type& t, vector_type& tA) = 0 ;
    virtual void update_diagonal(vector_type& t, vector_type& tA) = 0 ;
    virtual void orthogonalize(vector_type& t, vector_type& tA) = 0 ;
    // Attributes
protected:
    // "Real" data
    vector_space* vspace_reference_ ;
    vector_space* vspace_additional_ ;
    vector_scalar* diagonal_elements_ ;
    // Internal data
    double thresh_refinement_ ;
    bool do_refinement_, has_additional_ ;
};

#include "dmrg/optimize/Orthogonalizer/orthogonalizer.cpp"

// -- Derived classes --
// Up to now, Gram-Schmidt and Biorthogonal basis construction is supported.

#include "dmrg/optimize/Orthogonalizer/GS_ortho.hpp"
#include "dmrg/optimize/Orthogonalizer/GS_ortho_modified.hpp"
#include "dmrg/optimize/Orthogonalizer/BI_ortho.hpp"

#endif

