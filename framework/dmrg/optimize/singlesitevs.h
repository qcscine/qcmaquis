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

#include "dmrg/optimize/utils/bound_database.h"

template<class Matrix, class SymmGroup>
class SingleSiteVS
{
    typedef typename std::vector< class std::vector< std::vector< class block_matrix< typename storage::constrained<Matrix>::type, SymmGroup> > > > vector_ortho_type ;
    typedef Boundary<typename storage::constrained<Matrix>::type, SymmGroup>                      boundary ; 
    typedef std::vector< boundary >                                                               boundaries_type ;
    typedef typename bound_database< MPS<Matrix, SymmGroup>, boundaries_type>::bound_database     bound_database ; 
public:
    // Initializer from a single MPSTns
    SingleSiteVS(MPSTensor<Matrix, SymmGroup> & m,
                 std::vector< class MPSTensor<Matrix, SymmGroup> > const & ortho_vecs)
            : ortho_vecs_(ortho_vecs)
    {
        N_root  = 1 ;
        N_ortho = 0 ;
        for (std::size_t k = 0; k < m.data().n_blocks(); ++k)
            N_ortho += num_rows(m.data()[k]) * num_cols(m.data()[k]) ;
        MPSTns_vec.push_back(&m) ;
    }
    // Initializer from a VectorSet object
    // Note that here two possibilities are present:
    SingleSiteVS(VectorSet<Matrix, SymmGroup> & vs,
                 std::vector< class MPSTensor<Matrix, SymmGroup> > const & ortho_vecs,
                 vector_ortho_type vec_sa_left,
                 vector_ortho_type vec_sa_right,
                 bound_database database)
            : ortho_vecs_(ortho_vecs)
            , vec_sa_left_(vec_sa_left)
            , vec_sa_right_(vec_sa_right)
            , database_(database)
    {
        // Loads basis vectors
        N_root = vs.n_vec;
        for (std::size_t k = 0 ; k < N_root ;  k++)
            MPSTns_vec.push_back(&(vs.MPSTns_SA[k])) ;
        // Loads orthogonal vectors
        N_ortho = 0 ;
        for (std::size_t k = 0; k < vs.MPSTns_SA[0].data().n_blocks(); ++k)
            N_ortho += num_rows(vs.MPSTns_SA[0].data()[k]) * num_cols(vs.MPSTns_SA[0].data()[k]) ;
        ortho_vecs_add_.resize(0) ;
    }
    // Function to access data
    friend std::size_t n_root(SingleSiteVS const & vs) { return vs.N_root ; }
    friend MPSTensor<Matrix, SymmGroup> new_vector(SingleSiteVS & vs) { return *(vs.MPSTns_vec[0]) ; }
    friend MPSTensor<Matrix, SymmGroup> new_vector(SingleSiteVS & vs, const std::size_t & k) { return *(vs.MPSTns_vec[k]) ; }
    // Function to perform orthogonalization
    void project(MPSTensor<Matrix, SymmGroup> & t) const
    {
        // Orthogonalization vs excited states
        for (typename std::vector<MPSTensor<Matrix, SymmGroup> >::const_iterator it = ortho_vecs_.begin(); it != ortho_vecs_.end(); ++it)
            t -= ietl::dot(*it,t)/ietl::dot(*it,*it)**it;
        // Orthogonalization vs already converged roots
        for (typename std::vector<MPSTensor<Matrix, SymmGroup> >::const_iterator it = ortho_vecs_add_.begin(); it != ortho_vecs_add_.end(); ++it) {
            t -= ietl::dot(*it, t)/ietl::dot(*it,*it) *(*it);
        }
    }
    //
    std::vector<double> coeff(const MPSTensor<Matrix, SymmGroup> & t)
    {
        std::vector<double> res ;
        for (typename std::vector<MPSTensor<Matrix, SymmGroup> >::const_iterator it = ortho_vecs_add_.begin(); it != ortho_vecs_add_.end(); ++it)
            res.push_back(ietl::dot(*it,t));
        return res ;
    }
    //
    void add_orthovec(const MPSTensor<Matrix, SymmGroup> & t, 
                      const std::size_t& idx, 
                      const std::size_t& site)
    {
        ortho_vecs_add_.resize(0) ;
        MPSTensor<Matrix, SymmGroup> tmp ;
        for (int i = 0; i < idx; i++) {
            tmp = contraction::site_ortho_boundaries(*(MPSTns_vec[idx]), t, vec_sa_left_[idx][i][site], vec_sa_right_[i][idx][site+1]);
            ortho_vecs_add_.push_back(tmp) ;
        }
    }
    //
    MPSTensor<Matrix, SymmGroup> return_orthovec(const MPSTensor<Matrix, SymmGroup> & t,
                                                 const std::size_t& idx,
                                                 const std::size_t& i_2ortho,
                                                 const std::size_t& site1,
                                                 const std::size_t& site2)
    {
        MPSTensor<Matrix, SymmGroup> tmp ;
        tmp = contraction::site_ortho_boundaries(*(MPSTns_vec[idx]), t, vec_sa_left_[idx][i_2ortho][site1], vec_sa_right_[i_2ortho][idx][site2+1]);
        return tmp ;
    };
private:
    bound_database database_ ;
    std::vector<MPSTensor<Matrix, SymmGroup>* > MPSTns_vec ;
    std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs_ ;
    std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs_add_ ;
    vector_ortho_type vec_sa_left_, vec_sa_right_ ;
    std::size_t N_ortho, N_root ;
};

#endif
