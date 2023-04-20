/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SYMMETRY_NU1_TEMPLATE_H
#define SYMMETRY_NU1_TEMPLATE_H

#include <iostream>
#include <vector>
#include <list>
#include <numeric>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>

// =================
//  NU1Charge class
// =================
// Identifies a quantum number for the NU1 symmetry point group,
// which is defined in the following.
// The NU1Charge is basically a vector of type S (default is integer)
// of N elements, where N is set in the CMake file. The i-th element
// of the vector corresponds to the number of particles in the sites
// of type i

template<int N, class S = int>
class NU1Charge
{
public:
    // N is the dimension of the vector (i.e. of the group). It is made
    // both static and const
    static const int static_size = N;
    // -- Constructors --
    NU1Charge(S init = 0)
    {
        for (S i = 0; i < N; ++i)
            (*this)[i] = init;
    }
    NU1Charge(std::vector<S> const & rhs) : NU1Charge(0)
    {
        std::copy(rhs.begin(), rhs.end(), this->begin());
    }
    NU1Charge(std::vector<unsigned int> const & temp) : NU1Charge(0)
    {
        std::vector<int> rhs(temp.begin(), temp.end());
        std::copy(rhs.begin(), rhs.end(), this->begin());
    }

    NU1Charge(boost::array<S, N> const & rhs)
    {
        std::copy(rhs.begin(), rhs.end(), this->begin());
    }

    NU1Charge(std::array<S, N> const & rhs)
    {
        std::copy(rhs.begin(), rhs.end(), this->begin());
    }

    // -- Access methods --
    S * begin() { return &data_[0]; }
    S * end() { return &data_[N]; }
    S const * begin() const { return &data_[0]; }
    S const * end() const { return &data_[N]; }
    S & operator[](std::size_t p) { return data_[p]; }
    S const & operator[](std::size_t p) const { return data_[p]; }
    // -- Extracts info from archive
    template<class Archive>
    void save(Archive & ar) const
    {
        for (int i = 0; i < N; ++i)
            ar[boost::lexical_cast<std::string>(i)] << (*this)[i];
    }
    //
    template<class Archive>
    void load(Archive & ar)
    {
        for (int i = 0; i < N; ++i)
            ar[boost::lexical_cast<std::string>(i)] >> (*this)[i];
    }
    //
    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & data_;
    }

private:
    // -- MEMBERS --
    // The only member is the vector with the QN
    S data_[N];
};

namespace boost {
    template <int N, class S>
    class hash<NU1Charge<N, S> >{
        public :
            size_t operator()(NU1Charge<N, S> const &Charge ) const {
                return boost::hash_range(&Charge[0],&Charge[N]);
            }
    };

    template <int N, class S>
    class hash<std::pair<NU1Charge<N, S>, NU1Charge<N, S> > >{
        public :
            size_t operator()(std::pair<NU1Charge<N, S>, NU1Charge<N, S> > const &Pair_of_charge ) const {
                std::size_t seed = 0;
                for(int i(0); i<N; i++){
                    boost::hash_combine(seed, (Pair_of_charge.first)[i]);
                    boost::hash_combine(seed, (Pair_of_charge.second)[i]);
                }
                return seed;
            }
    };

    // Leon: Commented out as it breaks on later Boost versions!
    // Anyway to whoever did that: overloading functions in standard namespaces is bad practice!
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
std::ostream& operator<<(std::ostream& os, NU1Charge<N, S> const & c)
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

template<int N, class S>
inline bool operator<(NU1Charge<N, S> const & a, NU1Charge<N, S> const & b)
{
    return tpl_ops_<N, 0>().operator_lt(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator>(NU1Charge<N, S> const & a, NU1Charge<N, S> const & b)
{
    return tpl_ops_<N, 0>().operator_gt(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator==(NU1Charge<N, S> const & a, NU1Charge<N, S> const & b)
{
    return tpl_ops_<N, 0>().operator_eq(a.begin(), b.begin());
}

template<int N, class S>
inline bool operator!=(NU1Charge<N, S> const & a, NU1Charge<N, S> const & b)
{
    return !(a==b);
}

template<int N, class S>
NU1Charge<N, S> operator+(NU1Charge<N, S> const & a,
                          NU1Charge<N, S> const & b)
{
    NU1Charge<N, S> ret;
    tpl_ops_<N, 0>().operator_plus(a.begin(), b.begin(), ret.begin());
    return ret;
}

template<int N, class S>
NU1Charge<N, S> operator-(NU1Charge<N, S> const & rhs)
{
    NU1Charge<N, S> ret;
    tpl_ops_<N, 0>().operator_uminus(rhs.begin(), ret.begin());
    return ret;
}

template<int N, class S>
NU1Charge<N, S> operator/(NU1Charge<N, S> const & a, int n)
{
    NU1Charge<N, S> ret;
    tpl_ops_<N, 0>().operator_div(a.begin(), ret.begin(), n);
    return ret;
}
template<int N, class S>
NU1Charge<N, S> operator/(int n, NU1Charge<N, S> const & a) { return a/n; }

// ====================
//  NU1_TEMPLATE CLASS
// ====================
// Class associated to the NU1 symmetry point group.

template<int N, class S = int>
class NU1_template
{
public:
    typedef S subcharge;
    typedef NU1Charge<N, S> charge;
    typedef std::vector<charge> charge_v;

    static const charge IdentityCharge;
    static const bool finite = false;

    static subcharge particleNumber(charge rhs) { return std::accumulate(rhs.begin(), rhs.end(), 0); }

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

    static charge revert(charge& input)
    {
        charge output = charge(0) ;
        for (std::size_t idx = 0; idx < N; idx++)
            output[idx] = input[N-idx-1] ;
        return output ;
    }

    int get_dimension() { return N ; }
};


template<int N, class S> const typename NU1_template<N,S>::charge NU1_template<N,S>::IdentityCharge = typename NU1_template<N,S>::charge();

#endif
