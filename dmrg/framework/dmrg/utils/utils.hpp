/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <cstddef>
#include <complex>

#include "dmrg/utils/proc_statm.h"
#include "dmrg/utils/proc_status.h"

struct cmp_with_prefactor {
	static double prefactor;
	bool operator() (std::size_t i, std::size_t j) {
		bool ret = (i < j);
		if (ret) prefactor *= -1.;
		return ret;
	}
};

template<class T>
bool check_real(T x) { return true; }

template<class T>
bool check_real(std::complex<T> x)
{
    return std::imag(x)/std::real(x) < 1e-14 || std::imag(x) < 1e-14;
}


template <class InputIterator, class Predicate>
bool all_true (InputIterator first, InputIterator last, Predicate pred)
{
    bool allTrue = true;
    while (allTrue && first != last) 
        allTrue = pred(*first++);
    return allTrue;
}

template <class Pair>
struct compare_pair
{
    bool operator()(Pair const & i,
                    Pair const & j) const
    {
        if (i.first < j.first)
            return true;
        else if (i.first > j.first)
            return false;
        else
            return i.second < j.second;
    }
};

template <class Pair>
struct compare_pair_inverse
{
    bool operator()(Pair const & i,
                    Pair const & j) const
    {
        if (i.second < j.second)
            return true;
        else if (j.second < i.second)
            return false;
        else
            return i.first < j.first;
    }
};

#endif /* UTILS_HPP_ */
