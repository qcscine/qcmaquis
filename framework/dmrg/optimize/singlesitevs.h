/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *
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

template<class Matrix, class SymmGroup>
class SingleSiteVS
{
public:
    // Initializer from a single MPSTns
    SingleSiteVS(MPSTensor<Matrix, SymmGroup> const & m,
                 std::vector< class MPSTensor<Matrix, SymmGroup> > const & ortho_vecs)
            : MPSTns(m)
            , ortho_vecs_(ortho_vecs)
    {
        N_sa    = 0 ;
        N_ortho = 0 ;
        for (std::size_t k = 0; k < m.data().n_blocks(); ++k)
            N_ortho += num_rows(m.data()[k]) * num_cols(m.data()[k]);
        MPSTns_vec.resize(N_sa) ;
    }
    // Initializer from a VectorSet object
    SingleSiteVS(VectorSet<Matrix, SymmGroup> const & vs,
                 std::vector< class MPSTensor<Matrix, SymmGroup> > const & ortho_vecs)
            : MPSTns(vs.MPSTns_averaged)
            , ortho_vecs_(ortho_vecs)
    {
        N_sa    = vs.n_sa ;
        N_ortho = 0 ;
        for (std::size_t k = 0; k < vs.MPSTns_averaged.data().n_blocks(); ++k)
            N_ortho += num_rows(vs.MPSTns_averaged.data()[k]) * num_cols(vs.MPSTns_averaged.data()[k]);
        for (std::size_t k = 0 ; k < N_sa ;  k++)
            MPSTns_vec.push_back(vs.MPSTns_SA[k]) ;
    }
    // Function to access data
    friend MPSTensor<Matrix, SymmGroup> new_vector(SingleSiteVS const & vs) { return vs.MPSTns ; }
    friend std::size_t vec_dimension(SingleSiteVS const & vs)  { return vs.N ; }
    friend std::size_t n_sa_dimension(SingleSiteVS const & vs) { return vs.N_sa ; }
    friend MPSTensor<Matrix, SymmGroup> new_vector_sa(SingleSiteVS const & vs, const std::size_t & k) { return vs.MPSTns_vec[k] ; }
    // Function to perform orthogonalization
    void project(MPSTensor<Matrix, SymmGroup> & t) const
    {
        for (typename std::vector<MPSTensor<Matrix, SymmGroup> >::const_iterator it = ortho_vecs_.begin(); it != ortho_vecs_.end(); ++it)
            t -= ietl::dot(*it,t)/ietl::dot(*it,*it)**it;
    }
private:
    MPSTensor<Matrix, SymmGroup> MPSTns;
    std::vector< MPSTensor<Matrix, SymmGroup> > MPSTns_vec ;
    std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs_ ;
    std::size_t N_ortho, N_sa ;
};

#endif
