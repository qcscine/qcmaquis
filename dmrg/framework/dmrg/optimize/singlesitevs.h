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

#ifndef SINGLESITEVS_H
#define SINGLESITEVS_H

#include <dmrg/mp_tensors/contractions/abelian/special.hpp>
#include "dmrg/optimize/utils/orthogonalizer_collector.h"

template<class Matrix, class SymmGroup>
class SingleSiteVS
{
    typedef typename block_matrix< typename storage::constrained<Matrix>::type, SymmGroup>::block_matrix block_matrix ;
    typedef typename std::vector< class std::vector< std::vector< block_matrix > > >                     vector_ortho_type ;
    typedef typename MPS<Matrix, SymmGroup>::MPS                                                         MPS ;
    typedef typename MPSTensor<Matrix, SymmGroup>::MPSTensor                                             MPSTensor ;
    typedef typename VectorSet<Matrix, SymmGroup>::VectorSet                                             VectorSet ;
    typedef typename orthogonalizer_collector< MPSTensor >::orthogonalizer_collector                     orthogonalizer_collector ;
public:
    // Note that here two possibilities are present:
    SingleSiteVS(VectorSet & vs, vector_ortho_type & vec_sa_left, vector_ortho_type & vec_sa_right,
                 orthogonalizer_collector& ortho_vectors) :
              vec_sa_left_(vec_sa_left)
            , vec_sa_right_(vec_sa_right)
            , ortho_vec_(ortho_vectors)
    {
        // Loads basis vectors
        N_root = vs.n_vec;
        for (std::size_t k = 0 ; k < N_root ;  k++)
            MPSTns_vec.push_back(&(vs.MPSTns_SA[k])) ;
    }
    // Function to access data
    friend std::size_t n_root(SingleSiteVS const & vs) { return vs.N_root ; }
    MPSTensor new_vector(const std::size_t & k) { return *(this->MPSTns_vec[k]) ; }
    // Function to perform orthogonalization and to modify orthogonal vectors
    void project(MPSTensor & t) const { ortho_vec_.orthogonalize(t) ; } ;
    void add_within_vec(const MPSTensor & t) { ortho_vec_.add_within_vec(t) ; } ;
    void add_other_vec(const MPSTensor & t) { ortho_vec_.add_other_vec(t) ; } ;
    void clear_projector() { ortho_vec_.clear_local() ; } ;
    //
    // -- RETURN_ORTHOVEC --
    // Computes the vector to be used in the constrained optimization
    MPSTensor return_orthovec(const MPSTensor & t, const std::size_t& idx, const std::size_t& i_2ortho,
                              const std::size_t& site1, const std::size_t& site2)
    {
        MPSTensor tmp ;
        tmp = contraction::site_ortho_boundaries(*(MPSTns_vec[idx]), t, vec_sa_left_[idx][i_2ortho][site1],
                                                 vec_sa_right_[idx][i_2ortho][site2+1]);
        return tmp ;
    };
private:
    // +----------+
    //  ATTRIBUTES
    // +----------+
    orthogonalizer_collector ortho_vec_ ;
    std::vector<MPSTensor* > MPSTns_vec ;
    vector_ortho_type vec_sa_left_, vec_sa_right_ ;
    std::size_t N_root ;
};

// +------------------+
//  VECTORSPACE_TRAITS
// +------------------+
// Structure for getting types given by a MPSTensor

namespace ietl
{
    template<class Matrix, class SymmGroup>
    struct vectorspace_traits<SingleSiteVS<Matrix, SymmGroup> >
    {
        typedef MPSTensor<Matrix, SymmGroup>                             vector_type;
        typedef typename MPSTensor<Matrix, SymmGroup>::value_type        scalar_type;
        typedef typename MPSTensor<Matrix, SymmGroup>::magnitude_type    magnitude_type;
        typedef typename MPSTensor<Matrix, SymmGroup>::real_type         real_type ;
        typedef std::size_t size_type;
    };
}

#endif
