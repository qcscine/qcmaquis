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

#ifndef SYMMETRY_NU1DG_H
#define SYMMETRY_NU1DG_H

#include "utils/io.hpp"
#include <vector>
#include <list>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>
#include <boost/algorithm/string.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>

#include <alps/numeric/matrix.hpp>

#include <dmrg/utils/BaseParameters.h>

#include <dmrg/block_matrix/symmetry/nu1pg.h>
#include <dmrg/block_matrix/symmetry/dg_tables.h>

template<int N, class S>
class NU1DG;

template<int N, class S>
class NU1ChargeDG : public NU1ChargePG<N, S>
{
    typedef NU1ChargePG<N, S> base;

public:
    NU1ChargeDG(S init = 0) : base(init) {}
    NU1ChargeDG(boost::array<S, N> const & rhs) : base(rhs) {}

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
    class hash<NU1ChargeDG<N, S> >{
        public :
            size_t operator()(NU1ChargeDG<N, S> const &Charge ) const {
                return hash<NU1ChargePG<N, S> >()(Charge);
            }
    };

    template <int N, class S>
    class hash<std::pair<NU1ChargeDG<N, S>, NU1ChargeDG<N, S> > >{
        public :
            size_t operator()(std::pair<NU1ChargeDG<N, S>, NU1ChargeDG<N, S> > const &Pair_of_charge ) const {
                return hash<std::pair<NU1ChargePG<N, S>, NU1ChargePG<N, S> > >()(Pair_of_charge);
            }
    };

}

template<int N, class S>
inline bool operator<(NU1ChargeDG<N, S> const & a, NU1ChargeDG<N, S> const & b)
{
    return tpl_ops_pg_<N, 0>().operator_lt(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator>(NU1ChargeDG<N, S> const & a, NU1ChargeDG<N, S> const & b)
{
    return tpl_ops_pg_<N, 0>().operator_gt(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator==(NU1ChargeDG<N, S> const & a, NU1ChargeDG<N, S> const & b)
{
    return tpl_ops_pg_<N, 0>().operator_eq(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator!=(NU1ChargeDG<N, S> const & a, NU1ChargeDG<N, S> const & b)
{
    return !(a==b);
}

template<int N, class S>
NU1ChargeDG<N, S> operator+(NU1ChargeDG<N, S> const & a,
                       NU1ChargeDG<N, S> const & b)
{
    NU1ChargeDG<N, S> ret;
    tpl_arith_<NU1DG<N, S>, N, 0>().operator_plus(a.begin(), b.begin(), ret.begin());
    return ret;
}


template<int N, class S>
NU1ChargeDG<N, S> operator-(NU1ChargeDG<N, S> const & rhs)
{
    NU1ChargeDG<N, S> ret;
    tpl_arith_<NU1DG<N,S>, N, 0>().operator_uminus(rhs.begin(), ret.begin());
    return ret;
}

template<int N, class S>
NU1ChargeDG<N, S> operator/(NU1ChargeDG<N, S> const & a, int n)
{
    NU1ChargeDG<N, S> ret;
    tpl_arith_<NU1DG<N,S>, N, 0>().operator_div(a.begin(), ret.begin(), n);
    return ret;
}

template<int N, class S = int>
class NU1DG
{
    template<class G, int A, int B> friend struct tpl_arith_;
public:
    typedef S subcharge;
    typedef NU1ChargeDG<N, S> charge;
    typedef std::vector<charge> charge_v;
private:
    static alps::numeric::matrix<S> mult_table;
    static std::vector<S> adjoin_table;
    static std::size_t group_id;
    static subcharge max_irrep;
public:	
    static const charge IdentityCharge;
    static const bool finite = false;

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

    static void initialize_dg_tables(BaseParameters & parms)
    {
        // open integral_file
        std::string integral_file = parms["integral_file"];
        std::ifstream integral_stream;
        integral_stream.open(integral_file.c_str());
        // get first line of integral_file
        std::string line;
        std::getline(integral_stream, line);
        // split it
        std::vector<std::string> split_line;
        boost::split(split_line, line, boost::is_any_of("="));
        std::string grp_str = *(--split_line.end());
        // isolate group number and read it
        boost::erase_all(grp_str,",");
        std::istringstream iss(grp_str);
        iss >> group_id;
        // close integral_file
        integral_stream.close();

        // set up multiplication table according to group_id (same numbering as in Dirac)
        switch(group_id) {
            // C2h, D2h
            case 4:
                max_irrep    = 8;
                mult_table   = generate_mult_table_C2h<S>(max_irrep);
                adjoin_table = generate_adjoin_table_C2h<S>(max_irrep);
                break;
            // C2, D2, C2v, Cs
            case 5:
                max_irrep    = 4;
                mult_table   = generate_mult_table_Cs_C2<S>(max_irrep);
                adjoin_table = generate_adjoin_table_Cs_C2<S>(max_irrep);
                break;
            // Cs --> group_ID = 5
            case 6:
                throw std::runtime_error("Double group Cs has ID = 5\n");
            // Ci
            case 7:
                max_irrep    = 4;
                mult_table   = generate_mult_table_Ci<S>(max_irrep);
                adjoin_table = generate_adjoin_table_Ci<S>(max_irrep);
                break;
            // C1
            case 8:
                max_irrep    = 4;
                mult_table   = generate_mult_table_C1<S>(max_irrep);
                adjoin_table = generate_adjoin_table_C1<S>(max_irrep);
                break;
            // D2h spinfree
            case 9:
                throw std::runtime_error("D2h spinfree not implemented yet!\n");
            // Cinfv (C2v + lin. symmetry) --> mapped to C32
            case 10:
                max_irrep    = 128;
                mult_table   = generate_mult_table_C32<S>(max_irrep);
                adjoin_table = generate_adjoin_table_C32<S>(max_irrep);
                break;
            // Dinfh (D2h + lin. symmetry) --> mapped to C16h
            case 11:
                max_irrep    = 128;
                mult_table   = generate_mult_table_C16h<S>(max_irrep);
                adjoin_table = generate_adjoin_table_C16h<S>(max_irrep);
                break;
            default:
                throw std::runtime_error("Double group not known!\nAvailable double groups are C1, Ci, C2h, D2h, C2, D2, C2v, Cs, Cinfv, Dinfh.\n");
        }
    }

    static subcharge const get_max_irrep()
    {
        return max_irrep;
    }

    template<S> friend std::vector<S> default_dg_adjoin_table();
    template<S> friend alps::numeric::matrix<S> default_dg_mult_table();
};
  
template<class S>
std::vector<S> default_dg_adjoin_table()
{
	return std::vector<S>();
}
  
template<class S>
alps::numeric::matrix<S> default_dg_mult_table()
{
    return alps::numeric::matrix<S>();
}

template<int N, class S> const typename NU1DG<N,S>::charge NU1DG<N,S>::IdentityCharge = typename NU1DG<N,S>::charge();
template<int N, class S> alps::numeric::matrix<S> NU1DG<N,S>::mult_table = default_dg_mult_table<S>();
template<int N, class S> std::vector<S> NU1DG<N,S>::adjoin_table = default_dg_adjoin_table<S>();
template<int N, class S> std::size_t NU1DG<N,S>::group_id = 0;
template<int N, class S> typename NU1DG<N,S>::subcharge NU1DG<N,S>::max_irrep = 0;

typedef NU1DG<1> U1DG;

#endif
