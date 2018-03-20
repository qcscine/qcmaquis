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

#ifndef VECTORSET_H
#define VECTORSET_H

//
// Vector set structure
// --------------------
//
// This structure contains the MPS and the state-averaged states (for SA calculations).
// Generalizes the MPS object that was previously passed to the optimizer
//

template<class Matrix, class SymmGroup>
struct VectorSet
{
    typedef typename MPS<Matrix, SymmGroup>::MPS MPSTyp ;
    typedef typename MPSTensor<Matrix, SymmGroup>::MPSTensor MPSTns ;
    typedef typename std::vector< MPSTyp > MPSVec ;
    typedef typename std::vector< MPSTns > MPSTnsVec ;
    typedef typename std::size_t size_t  ;
    // Constructor
    VectorSet(const MPSVec& MPS_SA_, const size_t& i)
    {
        n_vec = MPS_SA_.size()   ;
        for (int j = 0 ; j < n_vec ; j++)
            MPSTns_SA.push_back(MPS_SA_[j][i]) ;
    } ;
    VectorSet(const MPSTnsVec& MPSTns_SA_)
    {
        n_vec = MPSTns_SA_.size()   ;
        for (int j = 0 ; j < n_vec ; j++)
            MPSTns_SA.push_back(MPSTns_SA_[j]) ;
    } ;
    // Attribute
    MPSTnsVec MPSTns_SA ;
    size_t    n_vec ;
};

#endif