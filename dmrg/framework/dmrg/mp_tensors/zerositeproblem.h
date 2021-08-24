/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2018-2020 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MAQUIS_DMRG_ZEROSITEPROBLEM_H
#define MAQUIS_DMRG_ZEROSITEPROBLEM_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/utils/parallel.hpp"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include <sys/time.h>

/**
 * @brief ZeroSiteProblem class
 * 
 * Similar to SiteProblem, but used only in the back-propagation step of TD-DMRG.
 */

template<class Matrix, class SymmGroup>
class ZeroSiteProblem
{
public:
    // Types definition
    using BlockMatrixType = block_matrix<Matrix, SymmGroup>;
    using BoundaryType = Boundary<typename storage::constrained<Matrix>::type, SymmGroup>;
    using MPOTensorType = MPOTensor<Matrix, SymmGroup>;
    /** @brief Class constructor */
    ZeroSiteProblem(const MPOTensorType& mpo_ten_left, const MPOTensorType& mpo_ten_right,
                    const BoundaryType& left, const BoundaryType& right) 
      : MPOTen_left_(mpo_ten_left), MPOTen_right_(mpo_ten_right), left_(left), right_(right) 
    {}

    /** @brief Method to apply an operator */
    BlockMatrixType apply(const BlockMatrixType& input_MPS) const {
      return contraction::Engine<Matrix, Matrix, SymmGroup>::zerosite_hamil2(input_MPS, left_, right_, MPOTen_left_, MPOTen_right_);
    }

    /** @brief Energy getter */
    auto get_energy(const BlockMatrixType& x)
    {
        auto y = this->apply(x);
        auto res = ietl::dot(x, y)/ietl::dot(x, x);
        return maquis::real(res);
    }

private:
    // -- Private attributes --
    const MPOTensorType& MPOTen_left_, MPOTen_right_;
    const BoundaryType& left_, right_;
};

#endif //MAQUIS_DMRG_ZEROSITEPROBLEM_H