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

#ifndef PARTIAL_OVERLAP_H
#define PARTIAL_OVERLAP_H

#include <sstream>
#include <algorithm>
#include <numeric>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>

#include "utils/traits.hpp"

#include <boost/ptr_container/ptr_vector.hpp>

template<class Matrix, class SymmGroup>
class partial_overlap
{
public:
    // +----------------+
    //  Types definition
    // +----------------+
    typedef contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>   contr ;
    typedef typename block_matrix<Matrix,SymmGroup>::block_matrix                                 bmatrix ;
    typedef typename std::size_t                                                                  dim_type ;
    typedef typename MPSTensor<Matrix,SymmGroup>::MPSTensor                                       MPSTensor ;
    typedef typename TwoSiteTensor<Matrix, SymmGroup>::TwoSiteTensor                              TSTensor ;
    typedef typename MPS<Matrix,SymmGroup>::MPS                                                   MPSWave ;
    typedef typename Matrix::value_type                                                           value_type ;
    typedef typename std::vector< bmatrix >                                                       vector_overlap ;
    // +------------+
    //  Constructors
    // +------------+
    partial_overlap();
    partial_overlap(const MPSWave& MPS, const MPSWave& MPS_reference);
    // +-------+
    //  Methods
    // +-------+
    // General methods, to get attributes of the partial overlap objects
    void update(const MPSWave& MPS, const dim_type& l, const int& direction) ;
    void update(const MPSWave& MPS, const dim_type& l1, const dim_type& l2, const int& direction) ;
    // Calculation of the overlap. Two class of methods are available
    value_type overlap(const dim_type& i);
    value_type overlap(const MPSTensor& MPSTns) ;
    // Prepares the object before the optimization
    void prepare (const MPSTensor& MPS, const dim_type& idx) ;
    void prepare (const MPSTensor& MPS, const dim_type& idx1, const dim_type& idx2) ;
private:
    // Private attributes
    vector_overlap data_left_ , data_right_ ;
    dim_type lattice_L_ ;
    MPSWave MPS_reference_ ;
    MPSTensor ortho_MPS_ ;
};

#include "partial_overlap.cpp"

#endif
