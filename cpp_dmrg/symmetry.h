#ifndef SYMMETRY_H
#define SYMMETRY_H

#include <iostream>
#include <vector>
#include <list>

template<class C> class NullMap
{
public:
	int operator[](C) const { return 0; }
    
#ifdef HAVE_ALPS_HDF5
    void serialize(alps::hdf5::oarchive & ar) const {
        ar << alps::make_pvp("dummy", std::string("This is a dummy."));
    }
    void serialize(alps::hdf5::iarchive & ar) { }
#endif
};

class NullGroup
{
public:
	typedef enum { Plus } charge;
	typedef std::vector<charge> charge_v;
	typedef NullMap<charge> map;
	static map get_map(const charge_v&) { return NullMap<charge>(); }
	static const charge SingletCharge = Plus;
	static inline charge fuse(charge a, charge b) { return Plus; }
	template<int R> static charge fuse(boost::array<charge, R>) { return Plus; }
	
	static charge_v default_charges()
	{
		charge_v ret(1); ret[0] = Plus; return ret;
	}
	
	static bool overflow(charge) { return false; }
};

#ifdef HAVE_ALPS_HDF5
inline alps::hdf5::oarchive & serialize(alps::hdf5::oarchive & ar,
                                        std::string const & p,
                                        NullGroup::charge const & v)
{
    int k = 1;
    ar << alps::make_pvp(p, k);
    return ar;
}
        
inline alps::hdf5:: iarchive & serialize(alps::hdf5::iarchive & ar,
                                         std::string const & p,
                                         NullGroup::charge & v) {
    v = NullGroup::Plus;
    return ar;
}
#endif

inline NullGroup::charge operator-(NullGroup::charge a) { return a; }

template<typename key_type>
class ZtwoMap
{
public:
	ZtwoMap(bool _zero = false) : zero(_zero)
	{
	}
	
	int operator[](const key_type& k) const
	{
		if (zero)
			return 0;
		else
			return k;
	}
    
#ifdef HAVE_ALPS_HDF5
    void serialize(alps::hdf5::oarchive & ar) const
    {
        using namespace alps;
        ar << make_pvp("zero", static_cast<int>(zero));
    }
    
    void serialize(alps::hdf5::iarchive & ar)
    {
        using namespace alps;
        int tmp;
        ar >> make_pvp("zero", tmp);
        zero = tmp;
    }
#endif
	
private:
	bool zero;
};

class Ztwo
	{
	public:
		typedef enum { Plus, Minus } charge;
		typedef std::vector<charge> charge_v;
		typedef ZtwoMap<charge> map;
		static map get_map(const charge_v &v)
		{
			return map(v.size() == 1);
		}
		
		static const charge SingletCharge = Plus;
		
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
		
		static charge_v default_charges()
		{
			charge_v ret(2);
			ret[0] = Plus;
			ret[1] = Minus;
			return ret;
		}
	};	

#ifdef HAVE_ALPS_HDF5
alps::hdf5::oarchive & serialize(alps::hdf5::oarchive & ar,
                                 std::string const & p,
                                 Ztwo::charge const & v);
alps::hdf5::iarchive & serialize(alps::hdf5::iarchive & ar,
                                 std::string const & p,
                                 Ztwo::charge & v);
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

// this is a map with HUGE memory, but only little performance overhead
template<class T>
class U1map
{
public:
	U1map() { }
	
	U1map(const std::vector<T> cv)
	{
		assert(cv.size() > 0);
		
		T _max = cv[0], _min = cv[0];
		for (typename std::vector<T>::const_iterator it = cv.begin()+1;
             it != cv.end();
             it++)
		{
			_max = std::max(*it, _max);
			_min = std::min(*it, _min);
		}
		
		_map.resize(_max-_min+1);
		offset = -_min;
		
		for (unsigned int i = 0; i < cv.size(); i++)
			_map[cv[i]+offset] = i;
	}
	
	inline int operator[](int q) const
    {
        return _map[q+offset];
    }
    
#ifdef HAVE_ALPS_HDF5
    void serialize(alps::hdf5::oarchive & ar) const
    {
        using namespace alps;
        ar << make_pvp(std::string("_map"), _map);
        ar << make_pvp(std::string("offset"), offset);
    }
    
    void serialize(alps::hdf5::iarchive & ar)
    {
        using namespace alps;
        ar >> make_pvp(std::string("_map"), _map);
        ar >> make_pvp(std::string("offset"), offset);
    }
#endif
    
private:
	std::vector<T> _map;
	int offset;
};

class U1
{
public:
	typedef int charge;
	typedef std::vector<charge> charge_v;
	typedef U1map<int> map;
	
	static map get_map(const charge_v &v) { return map(v); }

	static const charge SingletCharge = 0;
	
	static charge fuse(charge a, charge b) { return a + b; }
	
	template<int R> static charge fuse(const boost::array<charge, R> &v)
	{
		charge ret = 0;
		for (int i = 0; i < R; i++)
			ret += v[i];
		return ret;
	}
	
	static charge_v default_charges() { return charge_v(1,SingletCharge); }
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
    void serialize(alps::hdf5::oarchive & ar) const { ar << alps::make_pvp("c", c_); }
    void serialize(alps::hdf5::iarchive & ar) { ar >> alps::make_pvp("c", c_); }
#endif
    
protected:
    unsigned int c_;
};

template<int Q>
class ZqMap
{
public:
    template<class Iterable>
    ZqMap(Iterable v) : map_(Q, 0)
    {
        int c = 0;
        for (typename Iterable::iterator it = v.begin();
             it != v.end();
             ++it)
            map_[it->get()] = c++;
    }
    
    ZqMap() { }
    
    inline int operator[](ZqCharge<Q> q) const
    {
        return map_[q.get()];
    }
    
#ifdef HAVE_ALPS_HDF5
    void serialize(alps::hdf5::oarchive & ar) const { ar << alps::make_pvp("map_", map_); }
    void serialize(alps::hdf5::iarchive & ar) { ar >> alps::make_pvp("map_", map_); }
#endif
    
protected:
    std::vector<int> map_;
};

template<int Q>
class Zq
{
public:
    typedef ZqCharge<Q> charge;
    typedef std::vector<charge> charge_v;
    typedef ZqMap<Q> map;
    
    static map get_map(const charge_v &v) { return map(v); }
    
    static const charge SingletCharge;
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
    
    static charge_v default_charges() {
        std::vector<charge> ret(Q);
        for (int i = 0; i < Q; ++i)
            ret[i] = i;
        return ret;
    }
};

template<int Q> const typename Zq<Q>::charge Zq<Q>::SingletCharge = ZqCharge<Q>(0);

#endif
