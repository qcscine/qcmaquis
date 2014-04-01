/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <cstddef>
#include <complex>

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

#endif /* UTILS_HPP_ */
