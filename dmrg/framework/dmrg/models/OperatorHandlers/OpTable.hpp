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

#ifndef OPTABLE_HPP
#define OPTABLE_HPP

template <class Matrix, class SymmGroup>
typename OPTable<Matrix, SymmGroup>::tag_type
OPTable<Matrix, SymmGroup>::register_op(op_t const & op_)
{
    tag_type ret = this->size();
    this->push_back(op_);
    return ret;
}

template <class Matrix, class SymmGroup>
std::pair<typename OPTable<Matrix, SymmGroup>::tag_type, typename OPTable<Matrix, SymmGroup>::mvalue_type>
OPTable<Matrix, SymmGroup>::checked_register(op_t const& sample)
{
    std::pair<bool, mvalue_type> cmp_result;
    typename std::vector<op_t>::iterator it_pt = this->begin();
    for (; it_pt != this->end(); ++it_pt) {
        cmp_result = tag_detail::equal(*it_pt, sample);
        if (cmp_result.first)
            break;
    }
    if (it_pt == this->end()) {
        return std::make_pair(this->register_op(sample), 1.0);
    } else
        return std::make_pair(it_pt - this->begin(), cmp_result.second);

}

template <class Matrix, class SymmGroup>
bool OPTable<Matrix, SymmGroup>::hasRegistered(const op_t& sample)
{
    std::pair<bool, mvalue_type> cmp_result;
    typename std::vector<op_t>::iterator it_pt = this->begin();
    for (; it_pt != this->end(); ++it_pt) {
        cmp_result = tag_detail::equal(*it_pt, sample);
        if (cmp_result.first)
            break;
    }
    return !(it_pt == this->end());
}

#endif