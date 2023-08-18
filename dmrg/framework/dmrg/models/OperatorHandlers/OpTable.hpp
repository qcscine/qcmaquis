/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef OPTABLE_HPP
#define OPTABLE_HPP

#include "dmrg/models/OperatorHandlers/OpTable.h"
#include <algorithm>

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type
OPTable<Matrix, SymmGroup>::register_op(op_t const & op_)
{
    tag_type ret = this->size();
    this->push_back(op_);
    return ret;
}

/**
 * @brief Checkes if `sample` operator already exists in the table.
 *        This is achieved by comparing if any of the operators in the table
 *        are simply scaled versions of `sample`.
 *
 * @return Returns a pair containing the index of the operator in the table and
 *         the scale factor.
 */
template <class Matrix, class SymmGroup>
std::pair<typename OPTable<Matrix, SymmGroup>::tag_type, typename OPTable<Matrix, SymmGroup>::mvalue_type>
OPTable<Matrix, SymmGroup>::checked_register(op_t const& sample)
{
    mvalue_type scale_factor;
    auto get_scale = [&sample, &scale_factor](const op_t& op) {
      std::pair<bool, mvalue_type> cmp_result = tag_detail::is_scaled(op, sample);
      scale_factor = cmp_result.second;
      return cmp_result.first;
    };

    if (auto op_it = find_if(this->begin(), this->end(), get_scale); op_it == this->end()) {
        // if operator not found, register new operator
        return std::make_pair(this->register_op(sample), 1.0);
    } else {
        // if operator found, register the existing operator withh the new scale factor
        return std::make_pair(std::distance(this->begin(), op_it), scale_factor);
    }
}

/**
 * @brief Checks if `sample` operator already exists in the table.
 */
template <class Matrix, class SymmGroup>
bool OPTable<Matrix, SymmGroup>::hasRegistered(const op_t& sample)
{
    auto op_it = find_if(this->begin(), this->end(), [&sample](const op_t& op) {
        return tag_detail::is_scaled(op, sample).first;
    });
    bool is_found = op_it != this->end();
    return is_found;
}

#endif
