/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
