/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#ifndef SYMMETRY_U1LPG_H
#define SYMMETRY_U1LPG_H

#include <iostream>
#include <vector>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>

#include <alps/numeric/matrix.hpp>

template<class S>
class U1LPG_tpl;

template<class S>
class U1ChargeLPG
{
	public:
	U1ChargeLPG(S init = 0)
	{
		(*this)[0] = init;
		(*this)[1] = 0;
	}

	U1ChargeLPG(boost::array<S, 2> const & rhs)
	{
		std::copy(rhs.begin(), rhs.end(), this->begin());
	}

	S * begin() { return &data_[0]; }
	S * end() { return &data_[2]; }

	S const * begin() const { return &data_[0]; }
	S const * end() const { return &data_[2]; }

	S & operator[](std::size_t p) { return data_[p]; }
	S const & operator[](std::size_t p) const { return data_[p]; }

	template<class Archive>
	void save(Archive & ar) const
	{
		ar[boost::lexical_cast<std::string>(0)] << (*this)[0];
		ar[boost::lexical_cast<std::string>(1)] << (*this)[1];
	}

	template<class Archive>
	void load(Archive & ar) 
	{
		ar[boost::lexical_cast<std::string>(0)] >> (*this)[0];
		ar[boost::lexical_cast<std::string>(1)] >> (*this)[1];
	}
	
	template<class Archive>
	void serialize(Archive & ar, const unsigned int version) 
	{
		ar & data_;
	}

	private:
	S data_[2];
};

namespace boost{
	template<class S>
	class hash<U1ChargeLPG<S> > {
		public:
			size_t operator()(U1ChargeLPG<S> const & Charge) const {
				return boost::hash_range(&Charge[0], &Charge[1]);
			}
	};

	template<class S>
	class hash<std::pair<U1ChargeLPG<S>, U1ChargeLPG<S> > >{
		public:
			size_t operator()(std::pair<U1ChargeLPG<S>, U1ChargeLPG<S> > const & Pair_of_charge) const {
				std::size_t seed = 0;
				boost::hash_combine(seed, (Pair_of_charge.first)[0]);
				boost::hash_combine(seed, (Pair_of_charge.second)[0]);
				boost::hash_combine(seed, (Pair_of_charge.first)[1]);
				boost::hash_combine(seed, (Pair_of_charge.second)[1]);
				return seed;
			}
	};
}

template<class S>
std::ostream& operator<<(std::ostream& os, U1ChargeLPG<S> const & c)
{
	os << "<" << c[0] << "," << c[1] << ">";
	return os;
}

struct tpl_ops_lpg_
{
	template<typename T>
	bool operator_lt(T const * a, T const * b) const
	{
		if (a[0] < b[0])
			return true;
		else if (a[0] > b[0])
			return false;
		else
			return (a[1] < b[1]);
	}

	template<typename T>
	bool operator_gt(T const * a, T const * b) const
	{
		if (a[0] > b[0])
			return true;
		else if (a[0] < b[0])
			return false;
		else
			return (a[1] > b[1]);
	}

	template<typename T>
	bool operator_eq(T const * a, T const * b) const
	{
		if (a[0] != b[0])
			return false;
		else
			return (a[1] == b[1]);
	}

};

template<class G>
struct tpl_arith_lpg_
{
	template<typename T>
	void operator_plus(T const * a, T const * b, T * c) const
	{
		c[0] = a[0] + b[0];
		c[1] = G::mult_table(a[1], b[1]);
	}

	template<typename T>
	void operator_uminus(T const * a, T * b) const
	{
		b[0] = -a[0];
		b[1] = G::adjoin(a[1]);
	}

	template<typename T>
	void operator_div(T const * a, T * b, int n) const
	{
		b[0] = a[0]/n;
	}
};

template<class S>
inline bool operator<(U1ChargeLPG<S> const & a, U1ChargeLPG<S> const & b)
{
	return tpl_ops_lpg_().operator_lt(a.begin(), b.begin());
}

template<class S>
inline bool operator>(U1ChargeLPG<S> const & a, U1ChargeLPG<S> const & b)
{
	return tpl_ops_lpg_().operator_gt(a.begin(), b.begin());
}

template<class S>
inline bool operator==(U1ChargeLPG<S> const & a, U1ChargeLPG<S> const & b)
{
	return tpl_ops_lpg_().operator_eq(a.begin(), b.begin());
}

template<class S>
inline bool operator!=(U1ChargeLPG<S> const & a, U1ChargeLPG<S> const & b)
{
	return !(a==b);
}

template<class S>
U1ChargeLPG<S> operator+(U1ChargeLPG<S> const & a, U1ChargeLPG<S> const & b)
{
	U1ChargeLPG<S> ret;
	tpl_arith_lpg_<U1LPG_tpl<S> >().operator_plus(a.begin(), b.begin(), ret.begin());
	return ret;
}

template<class S>
U1ChargeLPG<S> operator-(U1ChargeLPG<S> const & rhs)
{
	U1ChargeLPG<S> ret;
	tpl_arith_lpg_<U1LPG_tpl<S> >().operator_uminus(rhs.begin(), ret.begin());
	return ret;
}

template<class S>
U1ChargeLPG<S> operator/(U1ChargeLPG<S> const & a, int n)
{
	U1ChargeLPG<S> ret;
	tpl_arith_lpg_<U1LPG_tpl<S> >().operator_div(a.begin(), ret.begin(), n);
	return ret;
}


template<class S>
class U1LPG_tpl
{
public:
	typedef S subcharge;
	typedef U1ChargeLPG<S> charge;
	typedef std::vector<charge> charge_v;

	static const charge IdentityCharge;
    static const bool finite = false;
	static const alps::numeric::matrix<S> mult_table;
	static const std::vector<S> adjoin_table;

	static subcharge adjoin(subcharge I)
	{
		return adjoin_table[I];
	}
	
	static charge fuse(charge a, charge b) { return a + b; }
	
	template<int R> static charge fuse(boost::array<charge, R> const & v)
	{
		charge ret = v[0];
		for (int i = 1; i < R; ++i)
			ret = fuse(ret, v[i]);
		return ret;
	}
};

template<class S>
std::vector<S> generate_lpg_adjoin()
{
    int num_irreps = 128;
    std::vector<S> adjoin_table(num_irreps);

    // boson irreps
    adjoin_table[0] = 0;
    adjoin_table[63] = 63;

    for (int i = 1; i < num_irreps/2 - 1; ++i) {
        adjoin_table[i] = i + pow(-1,i+1);
    }

    // fermion irreps
    for (int i = num_irreps/2; i < num_irreps; ++i) {
        adjoin_table[i] = i + pow(-1,i);
    }

    return adjoin_table;
}

template<class S>
alps::numeric::matrix<S> generate_lpg_mult_table()
{
    // ************* TO DO ***************
    // check double group --> insert input param in the function?
    // go in the right case for the double group
    // initialize variables:
    // number of irreps, mult table, inverse and adjoints elements?,
    //
    // ***********************************

    // Cinfv double group mapped to C64
    // inverse and adjoint elements not implemented --> sebastian didn't

    int num_irreps = 128;
    if(num_irreps/2 % 2 == 1){
        throw std::logic_error("Number of boson and fermion irreps must be even\n");}
    int shift = num_irreps/2;
    int irrep = 1;
    int mj = 0;
    std::vector<S> mj2rep(num_irreps+1);
    alps::numeric::matrix<S> mult_table(num_irreps,num_irreps);
    mj2rep[shift+mj] = irrep;

    // populate mj2rep vector with boson and fermion irreps
    for(mj = 2; mj <= num_irreps/2-2; mj+=2){
        irrep++;
        mj2rep[shift+mj] = irrep;
        irrep++;
        mj2rep[shift-mj] = irrep;
    }

    mj = num_irreps/2;
    irrep++;
    mj2rep[shift+mj] = irrep;
    mj2rep[shift-mj] = irrep;

    for(mj = 1; mj <= num_irreps/2-1; mj+=2){
        irrep++;
        mj2rep[shift+mj] = irrep;
        irrep++;
        mj2rep[shift-mj] = irrep;
    }

    // build multiplication table
    int mij = 0;
    int jrrep = 0;
    int ijrrep = 0;
    for(int mi = -num_irreps/2; mi <= num_irreps/2; mi++){
        for(mj = -num_irreps/2; mj <= num_irreps/2; mj++){
            mij = mi + mj;
            if(mij <  -num_irreps/2){mij = mij + num_irreps;}
            if(mij >   num_irreps/2){mij = mij - num_irreps;}
            if(mij == -num_irreps/2){mij = num_irreps/2;}
            irrep  = mj2rep[shift+mi];
            jrrep  = mj2rep[shift+mj];
            ijrrep = mj2rep[shift+mij];
            mult_table(irrep-1,jrrep-1) = ijrrep-1;
        }
    }

    return mult_table;
}

template<class S> const typename U1LPG_tpl<S>::charge U1LPG_tpl<S>::IdentityCharge = typename U1LPG_tpl<S>::charge();
template<class S> const alps::numeric::matrix<S> U1LPG_tpl<S>::mult_table = generate_lpg_mult_table<S>();
template<class S> const std::vector<S> U1LPG_tpl<S>::adjoin_table = generate_lpg_adjoin<S>();

typedef U1LPG_tpl<int> U1LPG;

#endif

