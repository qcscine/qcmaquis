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

#ifndef MAQUIS_CORRECTOR_H
#define MAQUIS_CORRECTOR_H

// Forward declaration

#include "correctionequation.h"

template<class MATRIX, class VecSpace>
class CorrectionEquation ;

// +------------------------+
//  General corrector object
// +------------------------+

template<class MATRIX, class VecSpace>
class Corrector {
public:
    // Types definition
    typedef CorrectionEquation<MATRIX, VecSpace>  base ;
    typedef typename base::scalar_type            scalar_type ;
    typedef typename base::vector_type            vector_type ;
    // Constructor and Desctructor
    virtual ~Corrector() {} ;
    // Actual method
    virtual vector_type apply(CorrectionEquation<MATRIX, VecSpace>* corr_eq, const vector_type& input) = 0 ;
    virtual void precondition(CorrectionEquation<MATRIX, VecSpace>* corr_eq, vector_type& input) = 0 ;
private:
    virtual void multiply_diagonal(CorrectionEquation<MATRIX, VecSpace>* corr_eq, vector_type& input) = 0 ;
} ;

#include "dmrg/optimize/CorrectionEquation/standardcorrection.hpp"
#include "dmrg/optimize/CorrectionEquation/skewcorrection.hpp"
#include "dmrg/optimize/CorrectionEquation/foldedcorrection.hpp"
#include "dmrg/optimize/CorrectionEquation/modifiedcorrection.hpp"

#endif
