/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelleb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SYMMETRY_NU1PG_H
#define SYMMETRY_NU1PG_H

#include <iostream>
#include <vector>
#include <list>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>


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

template<int N, int I>
struct tpl_arith_
{
    template<typename T>
    void operator_plus(T const * a, T const * b, T * c) const
    {
        c[I] = a[I] + b[I];
        tpl_arith_<N, I+1>().operator_plus(a, b, c);
    }
    
    template<typename T>
    void operator_uminus(T const * a, T * b) const
    {
        b[I] = -a[I];
        tpl_arith_<N, I+1>().operator_uminus(a, b);
    }

    template<typename T>
    void operator_div(T const * a, T * b, int n) const
    {
        b[I] = a[I]/n;
        tpl_arith_<N, I+1>().operator_div(a, b, n);
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

template<int N>
struct tpl_arith_<N, N>
{
    template<typename T>
    void operator_plus(T const *, T const *, T *) const { }
    
    template<typename T>
    void operator_uminus(T const *, T *) const { }

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
    tpl_arith_<N, 0>().operator_plus(a.begin(), b.begin(), ret.begin());
    return ret;
}

template<int N, class S>
NU1ChargePG<N, S> operator-(NU1ChargePG<N, S> const & rhs)
{
    NU1ChargePG<N, S> ret;
    tpl_arith_<N, 0>().operator_uminus(rhs.begin(), ret.begin());
    return ret;
}

template<int N, class S>
NU1ChargePG<N, S> operator/(NU1ChargePG<N, S> const & a, int n)
{
    NU1ChargePG<N, S> ret;
    tpl_arith_<N, 0>().operator_div(a.begin(), ret.begin(), n);
    return ret;
}
template<int N, class S>
NU1ChargePG<N, S> operator/(int n, NU1ChargePG<N, S> const & a) { return a/n; }


template<int N, class S = int>
class NU1PG
{
public:
    typedef NU1ChargePG<N, S> charge;
    typedef std::vector<charge> charge_v;
    
    static const charge IdentityCharge;
    static const bool finite = false;

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

template<int N, class S> const typename NU1PG<N,S>::charge NU1PG<N,S>::IdentityCharge = typename NU1PG<N,S>::charge();

typedef NU1PG<2> TwoU1PG;

#endif
