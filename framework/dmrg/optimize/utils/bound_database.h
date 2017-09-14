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

// Object used for multi-state targeting algorithms

#ifndef BOUNDARY_DATABASE_H
#define BOUNDARY_DATABASE_H

template<class MPS, class Bound>
class bound_database
{
//
// Public attributes
// -----------------
public:
    // Types definition
    typedef typename std::size_t               size_t ;
    typedef typename std::vector<Bound>        Bound_vector ; 
    typedef typename std::vector<MPS>          MPS_vector ;
    typedef typename MPS::data_t::value_type   MPSTensor ;
    // Constructors
    bound_database(void) ;
    bound_database(MPS_vector& MPSVec, MPS& MPSAverage, Bound_vector& bound_left, Bound_vector& bound_right, 
                   const int& sa_algorithm);
    // Public methods
    Bound*     get_boundaries_left(const size_t& idx) ;
    Bound*     get_boundaries_right(const size_t& idx) ;
    MPS*       get_mps(const size_t& idx) ;
    MPSTensor* get_mps(const size_t& idx, const size_t& site) ;
//
// Private attributes
// ------------------
private:
    // Main attribures
    std::vector< Bound* > vec_bound_left_, vec_bound_right_ ; 
    std::vector< MPS* > vec_MPS_ ;
    // Dimensions
public:
    int n_MPS_ ;
} ;

// 
// Constructor
// -----------
template<class MPS, class Bound>
bound_database<MPS, Bound>::bound_database(void)
{
    n_MPS_ = 0 ;
}

template<class MPS, class Bound>
bound_database<MPS, Bound>::bound_database(MPS_vector& MPSVec,
                                           MPS& MPSAverage, 
                                           Bound_vector& bound_left,
                                           Bound_vector& bound_right, 
                                           const int& sa_algorithm)
{
    // The initialization depends on the value of the sa_algorithm parameter
    n_MPS_ = MPSVec.size() ;
    vec_bound_left_.resize(n_MPS_) ;
    vec_bound_right_.resize(n_MPS_) ;
    vec_MPS_.resize(n_MPS_) ;
    // Case 1 - boundaries built from a specific state
    if (sa_algorithm >= -1 ) {
        if (sa_algorithm >= n_MPS_) { 
            throw std::runtime_error("sa_algorithm parameter must be <= number of SA states") ;
        } else {
            for (size_t idx = 0; idx < n_MPS_; idx++) {
                vec_bound_left_[idx]  = &(bound_left[0]) ;
                vec_bound_right_[idx] = &(bound_right[0]) ;
                if (sa_algorithm == -1)
                    vec_MPS_[idx] = &(MPSAverage) ;
                else
                    vec_MPS_[idx] = &(MPSVec[sa_algorithm]) ;
            }
        }
    } else if (sa_algorithm == -2) {
        for (size_t idx = 0; idx < n_MPS_; idx++) {
            vec_bound_left_[idx]  = &(bound_left[idx]) ;
            vec_bound_right_[idx] = &(bound_right[idx]) ;
            vec_MPS_[idx]         = &(MPSVec[idx]) ;
        }
    }
}

//
// Public methods
// --------------
template<class MPS, class Bound>
Bound* bound_database<MPS, Bound>::get_boundaries_left(const size_t& idx)
{
    return vec_bound_left_[idx] ;
}

template<class MPS, class Bound>
Bound* bound_database<MPS, Bound>::get_boundaries_right(const size_t& idx)
{
    return vec_bound_right_[idx] ;
}

template<class MPS, class Bound>
MPS* bound_database<MPS, Bound>::get_mps(const size_t& idx)
{
    return vec_MPS_[idx] ;
}

template<class MPS, class Bound>
typename bound_database<MPS, Bound>::MPSTensor* bound_database<MPS, Bound>::get_mps(const size_t& idx, 
                                                                                    const size_t& site)
{
    return &((*(vec_MPS_[idx]))[site]) ;
}

#endif

