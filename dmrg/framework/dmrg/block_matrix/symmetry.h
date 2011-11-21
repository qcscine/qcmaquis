/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <iostream>
#include <vector>
#include <list>

#include <boost/lexical_cast.hpp>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

#include <boost/lexical_cast.hpp>

class TrivialGroup
{
public:
	typedef enum { Plus } charge;
	static const charge IdentityCharge = Plus;
    
	static inline charge fuse(charge a, charge b) { return Plus; }
	template<int R> static charge fuse(boost::array<charge, R>) { return Plus; }
};

#ifdef HAVE_ALPS_HDF5
namespace alps {
    namespace hdf5 {
        
        inline void save(
                         alps::hdf5::archive & ar
                         , std::string const & path
                         , TrivialGroup::charge const & value
                         , std::vector<std::size_t> size = std::vector<std::size_t>()
                         , std::vector<std::size_t> chunk = std::vector<std::size_t>()
                         , std::vector<std::size_t> offset = std::vector<std::size_t>())
        {
            int k = 1;
            ar << alps::make_pvp(path, k);            
        }
        inline void load(
                         alps::hdf5::archive & ar
                         , std::string const & path
                         , TrivialGroup::charge & value
                         , std::vector<std::size_t> chunk = std::vector<std::size_t>()
                         , std::vector<std::size_t> offset = std::vector<std::size_t>())
        {
            value = TrivialGroup::Plus;
        }
        
    }
}
#endif

inline TrivialGroup::charge operator-(TrivialGroup::charge a) { return a; }

class Ztwo
	{
	public:
		typedef enum { Plus, Minus } charge;
		
		static const charge IdentityCharge = Plus;
		
		static inline charge fuse(charge a, charge b)
		{
			if (a == b)
				return Plus;
			else
				return Minus;
		}
		
		template<int R>
		static charge fuse(boost::array<charge, R> v)
		{
			// this operation actually could be rearranged into a tree
			for (int i = 1; i < R; i++)
				v[0] = fuse(v[0], v[i]);
			return v[0];
		}
	};	

#ifdef HAVE_ALPS_HDF5
void save(alps::hdf5::archive & ar,
          std::string const & p,
          Ztwo::charge const & v,
          std::vector<std::size_t> size = std::vector<std::size_t>(),
          std::vector<std::size_t> chunk = std::vector<std::size_t>(),
          std::vector<std::size_t> offset = std::vector<std::size_t>());
void load(alps::hdf5::archive & ar,
          std::string const & p,
          Ztwo::charge & v,
          std::vector<std::size_t> size = std::vector<std::size_t>(),
          std::vector<std::size_t> chunk = std::vector<std::size_t>(),
          std::vector<std::size_t> offset = std::vector<std::size_t>());
#endif

inline Ztwo::charge operator-(Ztwo::charge a) { return a; }

inline std::ostream& operator<<(std::ostream& ost, Ztwo::charge c)
{
	if (c == Ztwo::Plus)
		ost << "Plus";
	else if (c == Ztwo::Minus)
		ost << "Minus";
	else
		ost << "???";
	return ost;
}
inline std::ostream& operator<<(std::ostream& ost, const std::vector<Ztwo::charge> &c)
{
	ost << "[ ";
	for (std::vector<Ztwo::charge>::const_iterator it = c.begin();
		 it != c.end();
		 it++)
		ost << ", " << *it;
	ost << " ]";
	return ost;
}

class U1
{
public:
	typedef int charge;

	static const charge IdentityCharge = 0;
	
	static charge fuse(charge a, charge b) { return a + b; }
	
	template<int R> static charge fuse(const boost::array<charge, R> &v)
	{
		charge ret = 0;
		for (int i = 0; i < R; i++)
			ret += v[i];
		return ret;
	}
};

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
    
#ifdef HAVE_ALPS_HDF5
    void save(alps::hdf5::archive & ar) const { ar << alps::make_pvp("c", c_); }
    void load(alps::hdf5::archive & ar) { ar >> alps::make_pvp("c", c_); }
#endif
    
protected:
    unsigned int c_;
};

template<int Q>
class Zq
{
public:
    typedef ZqCharge<Q> charge;
    
    static const charge IdentityCharge;
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
    
private:
    int data_[N];
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
class NU1
{
public:
    typedef NU1Charge<N> charge;
    typedef std::vector<charge> charge_v;
    
    static const charge IdentityCharge;

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
