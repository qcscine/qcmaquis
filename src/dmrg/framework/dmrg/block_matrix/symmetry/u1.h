/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SYMMETRY_U1_H
#define SYMMETRY_U1_H

#include <iostream>
#include <vector>

#include <boost/functional/hash.hpp>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>


class U1
{
public:
	typedef int charge;

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
