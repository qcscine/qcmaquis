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

#ifndef ORTHOGONALIZERS_H
#define ORTHOGONALIZERS_H

// +---------------------+
//  Mofidied Gram-Schmidt
// +---------------------+

template<class VECTOR, class MAGNITUDE>
std::size_t gram_schmidt_orthogonalizer(std::vector<VECTOR> &V,
                                        std::vector<VECTOR> &t,
                                        const MAGNITUDE& thresh_ortho = 1.0E-10 ) {
    // Do some check
    typedef std::size_t size_t ;
    std::vector<size_t> vector_size ;
    size_t size_ortho = t.size() , n_lin_ind = 0 , jcont = 0 ;
    for (typename std::vector<VECTOR>::iterator it_t = t.begin(); it_t != t.end() ; it_t++) {
        for (typename std::vector<VECTOR>::iterator it_V = V.begin(); it_V != V.end(); it_V++)
            *it_t -= ietl::dot(*it_V, *it_t) * *it_V;
        if (ietl::two_norm(*it_t) > thresh_ortho) {
            n_lin_ind += 1 ;
            V.push_back(t[jcont] / ietl::two_norm(t[jcont]));
        } else {
            vector_size.push_back(jcont) ;
        }
        jcont += 1 ;
    }
    for (typename std::vector<size_t>::iterator it = vector_size.end() ; it != vector_size.begin() ; it-- )
        t.erase(t.begin() + *it) ;
    return n_lin_ind ;
};

// +-------------------------------------+
//  Mofidied Gram-Schmidt with refinement
// +-------------------------------------+

template<class VECTOR, class MAGNITUDE>
std::size_t gram_schmidt_orthogonalizer_refinement(std::vector<VECTOR> &V,
                                                   std::vector<VECTOR> &t,
                                                   const MAGNITUDE& thresh_refinement = 0.25,
                                                   const MAGNITUDE& thresh_ortho = 1.0E-10 ) {
    // Do some check
    typedef std::size_t size_t ;
    std::vector<size_t> vector_size ;
    MAGNITUDE tau ;
    size_t size_ortho = t.size() , n_lin_ind = 0 , jcont = 0 ;
    for (typename std::vector<VECTOR>::iterator it_t = t.begin(); it_t != t.end() ; it_t++) {
        tau = ietl::two_norm(*it_t);
        for (typename std::vector<VECTOR>::iterator it_V = V.begin(); it_V != V.end(); it_V++)
            *it_t -= ietl::dot(*it_V, *it_t) * *it_V;
        if (ietl::two_norm(*it_t) < thresh_refinement*tau )
            for (typename std::vector<VECTOR>::iterator it_V = V.begin(); it_V != V.end(); it_V++)
                *it_t -= ietl::dot(*it_V, *it_t) * *it_V;
        if (ietl::two_norm(*it_t) > thresh_ortho) {
            n_lin_ind += 1 ;
            V.push_back(t[jcont] / ietl::two_norm(t[jcont]));
        } else {
            vector_size.push_back(jcont) ;
        }
        jcont += 1 ;
    }
    for (typename std::vector<size_t>::iterator it = vector_size.end() ; it != vector_size.begin() ; it--)
        t.erase(t.begin() + *it) ;
    return n_lin_ind ;
};

// +----------------------------------------------------------------------------+
//  Mofidied Gram-Schmidt with refinement + additional vector space to transform
// +----------------------------------------------------------------------------+

template<class VECTOR, class MAGNITUDE>
std::size_t gram_schmidt_orthogonalizer_additional(std::vector<VECTOR> &V,
                                                   std::vector<VECTOR> &V_add,
                                                   std::vector<VECTOR> &t,
                                                   std::vector<VECTOR> &t_add,
                                                   const MAGNITUDE& thresh_refinement = 0.25,
                                                   const MAGNITUDE& thresh_ortho = 1.0E-10 ) {
    // Do some check
    typedef std::size_t size_t ;
    typedef typename std::vector<VECTOR>::iterator iter ;
    std::vector<size_t> vector_size ;
    size_t size_ortho = t.size() , n_lin_ind = 0 , jcont = 0 ;
    for (iter it_t = t.begin(), it_tA = t_add.begin(); it_t != t.end(), it_tA != t_add.end() ; it_t++, it_tA++) {
        for (iter it_V = V.begin(), it_VA = V_add.begin() ; it_V != V.end(), it_VA != V_add.end(); it_V++, it_VA++) {
            *it_t  -= ietl::dot(*it_V, *it_t) * *it_V  ;
            *it_tA -= ietl::dot(*it_V, *it_t) * *it_VA ;
        }
        if (ietl::two_norm(*it_t) > thresh_ortho) {
            n_lin_ind += 1 ;
            V.push_back(*it_t / ietl::two_norm(*it_t));
            V_add.push_back(*it_tA / ietl::two_norm(*it_t));
        } else {
            vector_size.push_back(jcont) ;
        }
        jcont += 1 ;
    }
    for (typename std::vector<size_t>::iterator it = vector_size.end() ; it != vector_size.begin() ; it--) {
        t.erase(t.begin() + *it);
        t_add.erase(t_add.begin() + *it);
    }
    return n_lin_ind ;
};

#endif
