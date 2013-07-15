/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SYMMETRY_ZQ_H
#define SYMMETRY_ZQ_H

#include <iostream>
#include <vector>
#include <list>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>


template<int Q>
class ZqCharge
{
    friend ZqCharge<Q> operator+(ZqCharge<Q> a, ZqCharge<Q> b)
    {
        return ZqCharge<Q>((a.c_+b.c_)%Q);
    }

    friend std::ostream& operator<<(std::ostream& os, ZqCharge<Q> c)
    {
        os << c.get();
        return os;
    }

    friend bool operator<(ZqCharge<Q> a, ZqCharge<Q> b)
    {
        return a.c_ < b.c_;
    }

    friend bool operator>(ZqCharge<Q> a, ZqCharge<Q> b)
    {
        return a.c_ > b.c_;
    }
    
public:
    ZqCharge(unsigned int c = 0) : c_(c) { }
    
    bool operator==(ZqCharge<Q> b) const { return c_ == b.c_; }
    bool operator!=(ZqCharge<Q> b) const { return c_ != b.c_; }
    ZqCharge operator-() const
    {
        if (c_ == 0)
            return 0;
        else
            return Q-c_;
    }
    
    unsigned int get() { return c_; }
    
    template<class Archive> void save(Archive & ar) const { ar["c"] << c_; }
    template<class Archive> void load(Archive & ar) { ar["c"] >> c_; }
    
    template <class Archive>
    void serialize(Archive & ar, const unsigned int version)
    {
        ar & c_;
    }
    
protected:
    unsigned int c_;
};

template<int Q>
class Zq
{
public:
    typedef ZqCharge<Q> charge;
    
    static const charge IdentityCharge;
    static const bool finite = true;

    static const int q = Q;
    
    static charge fuse(charge a, charge b)
    {
        return a+b;
    }
    
    template<int R> static charge fuse(const boost::array<charge, R> &v)
    {
        charge ret = v[0];
        for (int i = 1; i < R; i++)
            ret = fuse(ret, v[i]);
        return ret;
    }
};

template<int Q> const typename Zq<Q>::charge Zq<Q>::IdentityCharge = ZqCharge<Q>(0);

#endif
