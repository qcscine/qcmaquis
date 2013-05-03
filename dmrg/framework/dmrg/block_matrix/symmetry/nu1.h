/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SYMMETRY_NU1_H
#define SYMMETRY_NU1_H

#include <iostream>
#include <vector>
#include <list>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>


template<int N>
class NU1Charge
{   
public:
    NU1Charge(int init = 0)
    {
        for (int i = 0; i < N; ++i) (*this)[i] = init;
    }
    
    NU1Charge(boost::array<int, N> const & rhs)
    {
        std::copy(rhs.begin(), rhs.end(), this->begin());
    }
    
    int * begin() { return &data_[0]; }
    int * end() { return &data_[N]; }
    
    int const * begin() const { return &data_[0]; }
    int const * end() const { return &data_[N]; }
    
    int & operator[](std::size_t p) { return data_[p]; }
    int const & operator[](std::size_t p) const { return data_[p]; }
    
#ifdef HAVE_ALPS_HDF5
    void save(alps::hdf5::archive & ar) const
    {
        for (int i = 0; i < N; ++i)
            ar << alps::make_pvp(boost::lexical_cast<std::string>(i), (*this)[i]);
    }
    
    void load(alps::hdf5::archive & ar)
    {
        for (int i = 0; i < N; ++i)
            ar >> alps::make_pvp(boost::lexical_cast<std::string>(i), (*this)[i]);
    }
#endif
    
    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & data_;
    }

private:
    int data_[N];
};

namespace boost {
    template <int N>
    class hash<NU1Charge<N> >{
        public :
            size_t operator()(NU1Charge<N> const &Charge ) const {
                return boost::hash_range(&Charge[0],&Charge[N]);
            }
    };

    template <int N>
    class hash<std::pair<NU1Charge<N>, NU1Charge<N> > >{
        public :
            size_t operator()(std::pair<NU1Charge<N>,NU1Charge<N> > const &Pair_of_charge ) const {
                std::size_t seed = 0;
                for(int i(0); i<N; i++){
                    boost::hash_combine(seed, (Pair_of_charge.first)[i]);
                    boost::hash_combine(seed, (Pair_of_charge.second)[i]);
                }
                return seed;
            }
    };

    template <>
    class hash<std::pair<int,int>  >{
        public :
            size_t operator()(std::pair<int, int> const &Pair_of_charge ) const {
                return boost::hash_value(Pair_of_charge);
            }
    };

};

template<int N>
std::ostream& operator<<(std::ostream& os, NU1Charge<N> const & c)
{
    os << "<";
    for (int i = 0; i < N; ++i) {
        os << c[i];
        if (i+1 < N)
            os << ",";
    }
    os << ">";
    return os;
}

template<int N, int I>
struct tpl_ops_
{
    template<typename T>
    bool operator_lt(T const * a, T const * b) const
    {
        if (a[I] < b[I])
            return true;
        else if (a[I] > b[I])
            return false;
        else
            return tpl_ops_<N, I+1>().operator_lt(a, b);
    }
    
    template<typename T>
    bool operator_gt(T const * a, T const * b) const
    {
        if (a[I] > b[I])
            return true;
        else if (a[I] < b[I])
            return false;
        else
            return tpl_ops_<N, I+1>().operator_gt(a, b);
    }
    
    template<typename T>
    bool operator_eq(T const * a, T const * b) const
    {
        if (a[I] != b[I])
            return false;
        else
            return tpl_ops_<N, I+1>().operator_eq(a, b);
    }
    
    template<typename T>
    void operator_plus(T const * a, T const * b, T * c) const
    {
        c[I] = a[I] + b[I];
        tpl_ops_<N, I+1>().operator_plus(a, b, c);
    }
    
    template<typename T>
    void operator_uminus(T const * a, T * b) const
    {
        b[I] = -a[I];
        tpl_ops_<N, I+1>().operator_uminus(a, b);
    }

    template<typename T>
    void operator_div(T const * a, T * b, int n) const
    {
        b[I] = a[I]/n;
        tpl_ops_<N, I+1>().operator_div(a, b, n);
    }
};

template<int N>
struct tpl_ops_<N, N>
{
    template<typename T>
    bool operator_lt(T const *, T const *) const { return false; }
    
    template<typename T>
    bool operator_gt(T const *, T const *) const { return false; }
    
    template<typename T>
    bool operator_eq(T const *, T const *) const { return true; }
    
    template<typename T>
    void operator_plus(T const *, T const *, T *) const { }
    
    template<typename T>
    void operator_uminus(T const *, T *) const { }

    template<typename T>
    void operator_div(T const *, T *, int) const { }
};

template<int N>
inline bool operator<(NU1Charge<N> const & a, NU1Charge<N> const & b)
{
    return tpl_ops_<N, 0>().operator_lt(a.begin(), b.begin());
}

template<int N>
inline bool operator>(NU1Charge<N> const & a, NU1Charge<N> const & b)
{
    return tpl_ops_<N, 0>().operator_gt(a.begin(), b.begin());
}

template<int N>
inline bool operator==(NU1Charge<N> const & a, NU1Charge<N> const & b)
{
    return tpl_ops_<N, 0>().operator_eq(a.begin(), b.begin());
}

template<int N>
inline bool operator!=(NU1Charge<N> const & a, NU1Charge<N> const & b)
{
    return !(a==b);
}

template<int N>
NU1Charge<N> operator+(NU1Charge<N> const & a,
                       NU1Charge<N> const & b)
{
    NU1Charge<N> ret;
    tpl_ops_<N, 0>().operator_plus(a.begin(), b.begin(), ret.begin());
    return ret;
}

template<int N>
NU1Charge<N> operator-(NU1Charge<N> const & rhs)
{
    NU1Charge<N> ret;
    tpl_ops_<N, 0>().operator_uminus(rhs.begin(), ret.begin());
    return ret;
}

template<int N>
NU1Charge<N> operator/(NU1Charge<N> const & a, int n)
{
    NU1Charge<N> ret;
    tpl_ops_<N, 0>().operator_div(a.begin(), ret.begin(), n);
    return ret;
}
template<int N>
NU1Charge<N> operator/(int n, NU1Charge<N> const & a) { return a/n; }


template<int N>
class NU1
{
public:
    typedef NU1Charge<N> charge;
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

template<int N> const typename NU1<N>::charge NU1<N>::IdentityCharge = typename NU1<N>::charge();

typedef NU1<2> TwoU1;
typedef NU1<3> ThreeU1;

#endif
