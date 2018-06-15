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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef ORTHOGONALIZER_COLLECTOR_H
#define ORTHOGONALIZER_COLLECTOR_H

// ==========================
//  ORTHOGONALIZER_COLLECTOR
// ==========================
//
// The orthogonalizer collector class collects vectors against which orthogonalziation has to be
// performed.

#include <vector>

template<class MATRIX>
class orthogonalizer_collector {
public:
    // Types definition
    typedef typename std::vector<MATRIX>   vecspace ;
    typedef typename vecspace::const_iterator    ref_iterator ;
    // Constructor
    orthogonalizer_collector() ;
    explicit orthogonalizer_collector(vecspace& const_vec) ;
    orthogonalizer_collector(orthogonalizer_collector& rhs) ;
    orthogonalizer_collector(vecspace& const_vec, vecspace& within_vec, vecspace& other_vec) ;
    // Methods
    void add_within_vec(const MATRIX& to_add) ;
    void add_other_vec(const MATRIX& to_add) ;
    void clear_local() ;
    void clear_all() ;
    // Orthogonalization routine
    void orthogonalize(MATRIX& to_ortho) const ;
private:
    // -- ATTRIBUTES --
    vecspace orthovec_constant, orthovec_withinstate, orthovec_otherstate ;
    std::size_t dim_overall, dim_const, dim_withinstate, dim_otherstate ;
    // Constant parameter
    const double refinement_threshold = 0.25 ;
    const double normalize_threshold  = 1.0E-10 ;
    // Private functions
    bool normalize(MATRIX& to_norm) ;
};

#include "orthogonalizer_collector.hpp"

#endif // ORTHOGONALIZER_COLLECTOR_H
