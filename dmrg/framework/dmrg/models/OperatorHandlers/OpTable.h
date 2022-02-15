/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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
