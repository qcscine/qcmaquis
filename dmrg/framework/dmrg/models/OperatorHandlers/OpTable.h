/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef OPTABLE_H
#define OPTABLE_H

#include <vector>
#include <utility>
#include <stdexcept>
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/models/tag_detail.h"

template <class Matrix, class SymmGroup>
class OPTable : public std::vector<typename operator_selector<Matrix, SymmGroup>::type>
{
public:
    using tag_type = tag_detail::tag_type;
    using op_t = typename operator_selector<Matrix, SymmGroup>::type;

private:
    using mvalue_type = typename Matrix::value_type;

public:
    tag_type register_op(op_t const & op_);
    std::pair<tag_type, mvalue_type> checked_register(const op_t& sample);
    bool hasRegistered(const op_t& sample);
};

#include "OpTable.hpp"

#endif
