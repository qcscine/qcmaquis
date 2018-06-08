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
    Bound*     get_boundaries_left_sp(const size_t& idx, const size_t& k) ;
    Bound*     get_boundaries_right_sp(const size_t& idx, const size_t& k) ;
    float      get_coefficients(const size_t& idx, const size_t& k) ;
    MPS*       get_mps(const size_t& idx) ;
    MPSTensor* get_mps(const size_t& idx, const size_t& site) ;
    size_t     get_num_bound(const size_t& idx) ;
    void       print();
//
// Private attributes
// ------------------
private:
    // Main attribures
    std::vector< Bound* >                vec_bound_left_, vec_bound_right_  ;
    std::vector< std::vector< Bound*> >  vec_left_sp_, vec_right_sp_ ;
    std::vector< MPS* >                  vec_MPS_ ;
    std::vector< std::vector<float> >    matrix_coefficients_ ;
    // Dimensions
public:
    int n_MPS_ ;
} ;

// +-----------+
//  Constructor
// +-----------+
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
    vec_left_sp_.resize(n_MPS_) ;
    vec_right_sp_.resize(n_MPS_) ;
    vec_MPS_.resize(n_MPS_) ;
    matrix_coefficients_.resize(n_MPS_) ;
    // boundaries built from a specific state
    if (sa_algorithm >= 0 ) {
        if (sa_algorithm >= n_MPS_) {
            throw std::runtime_error("sa_algorithm parameter must be <= number of SA states") ;
        } else {
            for (size_t idx = 0; idx < n_MPS_; idx++) {
                vec_bound_left_[idx] = &(bound_left[0]);
                vec_bound_right_[idx] = &(bound_right[0]);
                vec_left_sp_[idx].push_back(&(bound_left[0])) ;
                vec_right_sp_[idx].push_back(&(bound_right[0])) ;
                vec_MPS_[idx] = &(MPSVec[idx]) ;
                matrix_coefficients_[idx].push_back(1.0) ;
            }
        }
    // State-specific boundaries for sa_algorithm =-2 (full state-specific optimisation)
    // and sa_algorithm =-1 (in this version, state-specific optimisation with truncation based on average MPStensors)
    } else if (sa_algorithm == -2 || sa_algorithm == -1) {
        for (size_t idx = 0; idx < n_MPS_; idx++) {
            vec_bound_left_[idx]  = &(bound_left[idx]) ;
            vec_bound_right_[idx] = &(bound_right[idx]) ;
            vec_left_sp_[idx].push_back(&(bound_left[idx])) ;
            vec_right_sp_[idx].push_back(&(bound_right[idx])) ;
            vec_MPS_[idx] = &(MPSVec[idx]) ;
            matrix_coefficients_[idx].push_back(1.0) ;
        }
    } else {
        throw std::runtime_error("sa_algorithm parameter not recognized") ;
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
Bound* bound_database<MPS, Bound>::get_boundaries_left_sp(const size_t& idx, const size_t& k)
{
    Bound* res ;
    if (idx >= n_MPS_) {
        throw std::runtime_error("error in the idx index");
    } else {
        if (k >= vec_left_sp_[idx].size()) {
            throw std::runtime_error("error in the k index");
        } else {
            res = vec_left_sp_[idx][k] ;
        }
    }
    return res ;
}

template<class MPS, class Bound>
Bound* bound_database<MPS, Bound>::get_boundaries_right(const size_t& idx)
{
    return vec_bound_right_[idx] ;
}

template<class MPS, class Bound>
Bound* bound_database<MPS, Bound>::get_boundaries_right_sp(const size_t& idx, const size_t& k)
{
    Bound* res ;
    if (idx >= n_MPS_) {
        throw std::runtime_error("error in the idx index");
    } else {
        if (k >= vec_right_sp_[idx].size()) {
            throw std::runtime_error("error in the k index");
        } else {
            res = vec_right_sp_[idx][k] ;
        }
    }
    return res ;
}

template<class MPS, class Bound>
float bound_database<MPS, Bound>::get_coefficients(const size_t &idx, const size_t &k)
{
    float res ;
    if (idx >= n_MPS_) {
        throw std::runtime_error("error in the idx index");
    } else {
        if (k >= vec_right_sp_[idx].size()) {
            throw std::runtime_error("error in the k index");
        } else {
            res = matrix_coefficients_[idx][k] ;
        }
    }
    return res ;
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

template<class MPS, class Bound>
typename bound_database<MPS, Bound>::size_t bound_database<MPS, Bound>::get_num_bound(const size_t& idx)
{
    size_t res ;
    if (idx >= this->n_MPS_ || idx < 0) {
        throw std::runtime_error("wrong index in bound_database") ;
    } else {
        res = matrix_coefficients_[idx].size() ;
    }
    return res ;
}

template<class MPS, class Bound>
void bound_database<MPS, Bound>::print()
{
//   using namespace std;
//   cout << "vec_bound_left_: " << endl;
//
//   for (auto x : vec_bound_left_)
//     for (auto i = 0; i < x->size(); i++)
//       (*x)[i].print();
//
//   cout << "vec_bound_right_: " << endl;
//   for (auto x : vec_bound_right_)
//     for (auto i = 0; i < x->size(); i++)
//       (*x)[i].print();
//
//   cout << "vec_left_sp_:" << endl;
//   for (auto x : vec_left_sp_)
//     for (auto y : x)
//       for (auto i = 0; i < y->size(); i++)
//         (*y)[i].print();
//
//   cout << "vec_right_sp_:" << endl;
//   for (auto x : vec_right_sp_)
//     for (auto y : x)
//       for (auto i = 0; i < y->size(); i++)
//         (*y)[i].print();

//   cout << "vec_MPS_:" << endl;
//   for (auto x: vec_MPS_)
//     for (auto i = 0; i < x->size(); i++)
//       cout << i << " " << (*x)[i] << endl;

//   cout << "Boundary coefficients:" << endl;
//   for (auto x : matrix_coefficients_)
//     for (auto y : x)
//       cout << y << " ";
//   cout << endl;
}
#endif

