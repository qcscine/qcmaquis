/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Sebastian Keller <sebkelleb@phys.ethz.ch>
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

#ifndef SYMMETRY_NU1LPG_H
#define SYMMETRY_NU1LPG_H

#include "utils/io.hpp"
#include <vector>
#include <list>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>

#include <alps/numeric/matrix.hpp>

#include <dmrg/block_matrix/symmetry/nu1pg.h>
#include <dmrg/block_matrix/symmetry/lpg_tables.h>

template<int N, class S>
class NU1LPG;

template<int N, class S>
class NU1ChargeLPG : public NU1ChargePG<N, S>
{
    typedef NU1ChargePG<N, S> base;

public:
    NU1ChargeLPG(S init = 0) : base(init) {}
    NU1ChargeLPG(boost::array<S, N> const & rhs) : base(rhs) {}

    S * begin() { return base::begin(); }
    S * end() { return base::end(); }

    S const * begin() const { return base::begin(); }
    S const * end() const { return base::end(); }

    S & operator[](std::size_t p) { return base::operator[](p); }
    S const & operator[](std::size_t p) const { return base::operator[](p); }

    template<class Archive>
    void save(Archive & ar) const
    {
        base::save(ar); 
    }

    template<class Archive>
    void load(Archive & ar)
    {
        base::load(ar);
    }

    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        base::serialize(ar, version);
    }

};

namespace boost {
    template <int N, class S>
    class hash<NU1ChargeLPG<N, S> >{
        public :
            size_t operator()(NU1ChargeLPG<N, S> const &Charge ) const {
                return hash<NU1ChargePG<N, S> >()(Charge);
            }
    };

    template <int N, class S>
    class hash<std::pair<NU1ChargeLPG<N, S>, NU1ChargeLPG<N, S> > >{
        public :
            size_t operator()(std::pair<NU1ChargeLPG<N, S>, NU1ChargeLPG<N, S> > const &Pair_of_charge ) const {
                return hash<std::pair<NU1ChargePG<N, S>, NU1ChargePG<N, S> > >()(Pair_of_charge);
            }
    };

}

template<int N, class S>
inline bool operator<(NU1ChargeLPG<N, S> const & a, NU1ChargeLPG<N, S> const & b)
{
    return tpl_ops_pg_<N, 0>().operator_lt(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator>(NU1ChargeLPG<N, S> const & a, NU1ChargeLPG<N, S> const & b)
{
    return tpl_ops_pg_<N, 0>().operator_gt(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator==(NU1ChargeLPG<N, S> const & a, NU1ChargeLPG<N, S> const & b)
{
    return tpl_ops_pg_<N, 0>().operator_eq(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator!=(NU1ChargeLPG<N, S> const & a, NU1ChargeLPG<N, S> const & b)
{
    return !(a==b);
}

template<int N, class S>
NU1ChargeLPG<N, S> operator+(NU1ChargeLPG<N, S> const & a,
                       NU1ChargeLPG<N, S> const & b)
{
    NU1ChargeLPG<N, S> ret;
    tpl_arith_<NU1LPG<N, S>, N, 0>().operator_plus(a.begin(), b.begin(), ret.begin());
    return ret;
}


template<int N, class S>
NU1ChargeLPG<N, S> operator-(NU1ChargeLPG<N, S> const & rhs)
{
    NU1ChargeLPG<N, S> ret;
    tpl_arith_<NU1LPG<N,S>, N, 0>().operator_uminus(rhs.begin(), ret.begin());
    return ret;
}

template<int N, class S>
NU1ChargeLPG<N, S> operator/(NU1ChargeLPG<N, S> const & a, int n)
{
    NU1ChargeLPG<N, S> ret;
    tpl_arith_<NU1LPG<N,S>, N, 0>().operator_div(a.begin(), ret.begin(), n);
    return ret;
}

template<int N, class S = int>
class NU1LPG
{
public:
    typedef S subcharge;
    typedef NU1ChargeLPG<N, S> charge;
    typedef std::vector<charge> charge_v;
	
	// TODO: Take it as an input
	static const std::size_t group = 1;

    static const charge IdentityCharge;
    static const bool finite = false;
    static const alps::numeric::matrix<S> mult_table;
    static const std::vector<S> adjoin_table;

    static subcharge adjoin(subcharge I)
    {
        return adjoin_table[I];
    }


    static charge fuse(charge a, charge b)
    {
        return a+b;
    }
    
    template<int R> static charge fuse(boost::array<charge, R> const & v)
    {
        charge ret = v[0];
        for (int i = 1; i < R; ++i)
            ret = fuse(ret, v[i]);
        return ret;
    }
};
  
template<class S>
std::vector<S> generate_lpg_adjoin_table(std::size_t group)
{
	std::vector<S> ret;
	// TODO: Add all the cases for all symmetry groups
	switch(group) {
		case(1):
			ret = generate_adjoin_table_Cinf<S>();
			break;
	}

	/*
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

    return adjoin_table;*/
	return ret;
}
  
template<class S>
alps::numeric::matrix<S> generate_lpg_mult_table(std::size_t group)
{
	alps::numeric::matrix<S> ret;
	// TODO: Add all the cases for all symmetry groups
	switch(group) {
		case(1):
			ret = generate_mult_table_Cinf<S>();
			break;
	}


    // ************* TO DO ***************
    // check double group --> insert input param in the function?
    // go in the right case for the double group
    // initialize variables:
    // number of irreps, mult table, inverse and adjoints elements?,
    //
    // ***********************************

    // Cinfv double group mapped to C64
    // inverse and adjoint elements not implemented --> sebastian didn't
   /*

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
*/
	return ret;
}

template<int N, class S> const typename NU1LPG<N,S>::charge NU1LPG<N,S>::IdentityCharge = typename NU1LPG<N,S>::charge();
template<int N, class S> const alps::numeric::matrix<S> NU1LPG<N,S>::mult_table = generate_lpg_mult_table<S>(NU1LPG<N,S>::group);
template<int N, class S> const std::vector<S> NU1LPG<N,S>::adjoin_table = generate_lpg_adjoin_table<S>(NU1LPG<N,S>::group);

typedef NU1LPG<1> U1LPG;

#endif
