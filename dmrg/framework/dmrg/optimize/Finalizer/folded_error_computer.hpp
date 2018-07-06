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

#ifndef MAQUIS_DMRG_FOLDED_ERROR_HPP
#define MAQUIS_DMRG_FOLDED_ERROR_HPP

#include "dmrg/optimize/Finalizer/error_computer.hpp"

// +-----------------------+
//  Folded corrector object
// +-----------------------+

template<class MATRIX, class VecSpace>
class FoldedError : public ErrorComputer<MATRIX, VecSpace> {
public:
    // Types definition
    typedef ErrorComputer<MATRIX, VecSpace>      base ;
    typedef typename base::scalar_type           scalar_type ;
    typedef typename base::vector_type           vector_type ;
    // Destructor
    FoldedError() {} ;
    ~FoldedError() {} ;
    // Routine used to apply the correction equation
    vector_type compute_error(Finalizer<MATRIX, VecSpace>* finalizer,
                              size_t const& i_state,
                              size_t const& idx)
    {
        vector_type junk = finalizer->get_uA(idx) ;
        junk -= finalizer->get_eigen(idx) * (finalizer->get_u(idx)) ;
        return junk ;
    }
};

#endif
