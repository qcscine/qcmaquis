/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef SYMMETRY_U1_H
#define SYMMETRY_U1_H

#include <iostream>
#include <vector>

#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>


class U1
{
public:
	typedef int charge;
    typedef int subcharge;

	static const charge IdentityCharge = 0;
    static const bool finite = false;
	
	static charge fuse(charge a, charge b) { return a + b; }
	
	template<int R> static charge fuse(const boost::array<charge, R> &v)
	{
		charge ret = 0;
		for (int i = 0; i < R; i++)
			ret += v[i];
		return ret;
	}
};

template <class Archive>
inline void serialize(Archive & ar, U1::charge & c, const unsigned int version)
{
    ar & c;
}

#endif
