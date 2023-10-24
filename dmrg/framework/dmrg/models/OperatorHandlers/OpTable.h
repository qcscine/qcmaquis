/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef OPTABLE_H
#define OPTABLE_H

#include <vector>
#include <utility>
#include <stdexcept>
#include "dmrg/models/tag_detail.h"

template <class Matrix, class SymmGroup>
class OPTable : public std::vector<typename operator_selector<Matrix, SymmGroup>::type>
{
public:
    typedef tag_detail::tag_type tag_type;
    typedef typename operator_selector<Matrix, SymmGroup>::type op_t;

private:
    typedef typename Matrix::value_type mvalue_type;

public:
    tag_type register_op(op_t const & op_);
    std::pair<tag_type, mvalue_type> checked_register(const op_t& sample);
    bool hasRegistered(const op_t& sample);
};

#include "OpTable.hpp"

#endif
