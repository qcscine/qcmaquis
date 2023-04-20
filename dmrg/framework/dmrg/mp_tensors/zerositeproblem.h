/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
