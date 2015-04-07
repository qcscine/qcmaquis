/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#ifndef SYMMETRY_NU1PG_H
#define SYMMETRY_NU1PG_H

#include "utils/io.hpp"
#include <vector>
#include <list>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>

#include <alps/numeric/matrix.hpp>

template<int N, class S>
class NU1PG;

template<int N, class S = int>
class NU1ChargePG
{   
public:
    NU1ChargePG(S init = 0)
    {
        for (S i = 0; i < N; ++i) (*this)[i] = init;
        (*this)[N] = 0;
    }
    
    NU1ChargePG(boost::array<S, N> const & rhs)
    {
        std::copy(rhs.begin(), rhs.end(), this->begin());
    }
    
    S * begin() { return &data_[0]; }
    S * end() { return &data_[N+1]; }

    S const * begin() const { return &data_[0]; }
    S const * end() const { return &data_[N+1]; }

    S & operator[](std::size_t p) { return data_[p]; }
    S const & operator[](std::size_t p) const { return data_[p]; }
   
    template<class Archive> 
    void save(Archive & ar) const
    {
        for (int i = 0; i < N+1; ++i)
            ar[boost::lexical_cast<std::string>(i)] << (*this)[i];
    }
    
    template<class Archive> 
    void load(Archive & ar)
    {
        for (int i = 0; i < N+1; ++i)
            ar[boost::lexical_cast<std::string>(i)] >> (*this)[i];
    }
    
    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & data_;
    }

private:
    S data_[N+1];
};

namespace boost {
    template <int N, class S>
    class hash<NU1ChargePG<N, S> >{
        public :
            size_t operator()(NU1ChargePG<N, S> const &Charge ) const {
                return boost::hash_range(&Charge[0],&Charge[N+1]);
            }
    };

    template <int N, class S>
    class hash<std::pair<NU1ChargePG<N, S>, NU1ChargePG<N, S> > >{
        public :
            size_t operator()(std::pair<NU1ChargePG<N, S>, NU1ChargePG<N, S> > const &Pair_of_charge ) const {
                std::size_t seed = 0;
                for(int i(0); i<N+1; i++){
                    boost::hash_combine(seed, (Pair_of_charge.first)[i]);
                    boost::hash_combine(seed, (Pair_of_charge.second)[i]);
                }
                return seed;
            }
    };

    /*
    template <>
    class hash<std::pair<int, int>  >{
        public :
            size_t operator()(std::pair<int, int> const &Pair_of_charge ) const {
                return boost::hash_value(Pair_of_charge);
            }
    };
    */

};

template<int N, class S>
std::ostream& operator<<(std::ostream& os, NU1ChargePG<N, S> const & c)
{
    os << "<";
    for (int i = 0; i < N+1; ++i) {
        os << c[i];
        if (i+1 < N+1)
            os << ",";
    }
    os << ">";
    return os;
}

template<int N, int I>
struct tpl_ops_pg_
{
    template<typename T>
    bool operator_lt(T const * a, T const * b) const
    {
        if (a[I] < b[I])
            return true;
        else if (a[I] > b[I])
            return false;
        else
            return tpl_ops_pg_<N, I+1>().operator_lt(a, b);
    }
    
    template<typename T>
    bool operator_gt(T const * a, T const * b) const
    {
        if (a[I] > b[I])
            return true;
        else if (a[I] < b[I])
            return false;
        else
            return tpl_ops_pg_<N, I+1>().operator_gt(a, b);
    }
    
    template<typename T>
    bool operator_eq(T const * a, T const * b) const
    {
        if (a[I] != b[I])
            return false;
        else
            return tpl_ops_pg_<N, I+1>().operator_eq(a, b);
    }
};

template<class G, int N, int I>
struct tpl_arith_
{
    template<typename T>
    void operator_plus(T const * a, T const * b, T * c) const
    {
        c[I] = a[I] + b[I];
        tpl_arith_<G, N, I+1>().operator_plus(a, b, c);
    }
    
    template<typename T>
    void operator_uminus(T const * a, T * b) const
    {
        b[I] = -a[I];
        tpl_arith_<G, N, I+1>().operator_uminus(a, b);
    }

    template<typename T>
    void operator_div(T const * a, T * b, int n) const
    {
        b[I] = a[I]/n;
        tpl_arith_<G, N, I+1>().operator_div(a, b, n);
    }

};

template<int N>
struct tpl_ops_pg_<N, N>
{
    template<typename T>
    bool operator_lt(T const * a , T const * b) const { return a[N] < b[N]; }
    
    template<typename T>
    bool operator_gt(T const * a, T const * b) const { return a[N] > b[N]; }
    
    template<typename T>
    bool operator_eq(T const * a, T const * b) const { return a[N] == b[N]; }
};

template<class G, int N>
struct tpl_arith_<G, N, N>
{
    template<typename T>
    void operator_plus(T const * a, T const * b, T * ret) const
    {
        ret[N] = G::mult_table(a[N], b[N]);
    }

    
    template<typename T>
    void operator_uminus(T const * a, T * b) const { b[N] = G::adjoin(a[N]); }

    template<typename T>
    void operator_div(T const *, T *, int) const { }
};

template<int N, class S>
inline bool operator<(NU1ChargePG<N, S> const & a, NU1ChargePG<N, S> const & b)
{
    return tpl_ops_pg_<N, 0>().operator_lt(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator>(NU1ChargePG<N, S> const & a, NU1ChargePG<N, S> const & b)
{
    return tpl_ops_pg_<N, 0>().operator_gt(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator==(NU1ChargePG<N, S> const & a, NU1ChargePG<N, S> const & b)
{
    return tpl_ops_pg_<N, 0>().operator_eq(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator!=(NU1ChargePG<N, S> const & a, NU1ChargePG<N, S> const & b)
{
    return !(a==b);
}

template<int N, class S>
NU1ChargePG<N, S> operator+(NU1ChargePG<N, S> const & a,
                       NU1ChargePG<N, S> const & b)
{
    NU1ChargePG<N, S> ret;
    tpl_arith_<NU1PG<N, S>, N, 0>().operator_plus(a.begin(), b.begin(), ret.begin());
    return ret;
}

template<int N, class S>
NU1ChargePG<N, S> operator-(NU1ChargePG<N, S> const & rhs)
{
    NU1ChargePG<N, S> ret;
    tpl_arith_<NU1PG<N,S>, N, 0>().operator_uminus(rhs.begin(), ret.begin());
    return ret;
}

template<int N, class S>
NU1ChargePG<N, S> operator/(NU1ChargePG<N, S> const & a, int n)
{
    NU1ChargePG<N, S> ret;
    tpl_arith_<NU1PG<N,S>, N, 0>().operator_div(a.begin(), ret.begin(), n);
    return ret;
}


template<int N, class S = int>
class NU1PG
{
public:
    typedef S subcharge;
    typedef NU1ChargePG<N, S> charge;
    typedef std::vector<charge> charge_v;

    static const charge IdentityCharge;
    static const bool finite = false;
    static const alps::numeric::matrix<S> mult_table;

    static subcharge particleNumber(charge rhs) { return std::accumulate(rhs.begin(), &rhs[N], 0); }

    static subcharge & irrep(charge & rhs) { return rhs[N]; }
    static subcharge const & irrep(charge const & rhs) { return rhs[N]; }

    static charge fuse(charge a, charge b)
    {
        return a+b;
    }
    
    static subcharge adjoin(subcharge I)
    {
        return I;
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
alps::numeric::matrix<S> generate_mult_table()
{
    alps::numeric::matrix<S> r(8,8);
    r(0,0) = 0; r(0,1) = 1; r(0,2) = 2; r(0,3) = 3;   r(0,4) = 4; r(0,5) = 5; r(0,6) = 6; r(0,7) = 7;
    r(1,0) = 1; r(1,1) = 0; r(1,2) = 3; r(1,3) = 2;   r(1,4) = 5; r(1,5) = 4; r(1,6) = 7; r(1,7) = 6;
    r(2,0) = 2; r(2,1) = 3; r(2,2) = 0; r(2,3) = 1;   r(2,4) = 6; r(2,5) = 7; r(2,6) = 4; r(2,7) = 5;
    r(3,0) = 3; r(3,1) = 2; r(3,2) = 1; r(3,3) = 0;   r(3,4) = 7; r(3,5) = 6; r(3,6) = 5; r(3,7) = 4;

    r(4,0) = 4; r(4,1) = 5; r(4,2) = 6; r(4,3) = 7;   r(4,4) = 0; r(4,5) = 1; r(4,6) = 2; r(4,7) = 3;
    r(5,0) = 5; r(5,1) = 4; r(5,2) = 7; r(5,3) = 6;   r(5,4) = 1; r(5,5) = 0; r(5,6) = 3; r(5,7) = 2;
    r(6,0) = 6; r(6,1) = 7; r(6,2) = 4; r(6,3) = 5;   r(6,4) = 2; r(6,5) = 3; r(6,6) = 0; r(6,7) = 1;
    r(7,0) = 7; r(7,1) = 6; r(7,2) = 5; r(7,3) = 4;   r(7,4) = 3; r(7,5) = 2; r(7,6) = 1; r(7,7) = 0;
    return r;
}

template<int N, class S> const typename NU1PG<N,S>::charge NU1PG<N,S>::IdentityCharge = typename NU1PG<N,S>::charge();
template<int N, class S> const alps::numeric::matrix<S> NU1PG<N,S>::mult_table = generate_mult_table<S>();

typedef NU1PG<2> TwoU1PG;

#endif
