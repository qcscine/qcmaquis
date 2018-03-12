#ifndef INDEXING_WRAPPERS_H
#define INDEXING_WRAPPERS_H

#include <iostream>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"

#include <boost/python.hpp>

using namespace boost::python;

template<class Enum>
struct EnumToClass
{
    EnumToClass() { }
    EnumToClass(int i) : value(static_cast<Enum>(i)) { }
    
    std::string str()
    {
        std::ostringstream r;
        r << value;
        return r.str();
    }
    
    void set(int i) { value = static_cast<Enum>(i); }
    
    friend EnumToClass operator-(EnumToClass cp)
    {
        cp.value = -cp.value;
        return cp;
    }
    
    Enum value;
    
    operator Enum() const { return value; }
};

template<class sgrp> struct charge_wrapped_as { };
template<> struct charge_wrapped_as<TrivialGroup> { typedef EnumToClass<TrivialGroup::charge> type; };
template<> struct charge_wrapped_as<Ztwo> { typedef EnumToClass<Ztwo::charge> type; };

template<class sgrp>
typename charge_wrapped_as<sgrp>::type wrapped_fuse(object o1, object o2)
{
    typedef typename charge_wrapped_as<sgrp>::type wt;
    typedef typename sgrp::charge charge;
    
    extract<wt> w1(o1), w2(o2);
    if (!(w1.check() && w2.check()))
        maquis::cerr << "Conversion error!" << std::endl;
    
    return sgrp::fuse(w1(), w2());
}

template<class sgrp>
struct wrapped_pair
{
    typedef std::pair<typename sgrp::charge, std::size_t> raw_pair;
    typedef typename charge_wrapped_as<sgrp>::type wrapped_charge;
    
    wrapped_pair() { }
    
    wrapped_pair(raw_pair const & p)
    : data_(p) { }
    
    /* this is not necessary, surprisingly
     wrapped_pair(object o1, object o2)
     {
     extract<std::size_t> w1(o1);
     extract<wrapped_charge> w2(o2);
     if (!(w1.check() && w2.check()))
     maquis::cerr << "Conversion error!" << std::endl;
     
     data_.first = w1();
     data_.second = w2();
     } */
    
    wrapped_pair(wrapped_charge c, std::size_t s)
    {
        data_.first = static_cast<typename sgrp::charge>(c);
        data_.second = s;
    }
    
    wrapped_charge get_first() const
    {
        return wrapped_charge(data_.first);
    }
    
    std::size_t get_second() const
    {
        return data_.second;
    }
    
    
    void set_first(object o)
    {
        extract<wrapped_charge> e(o);
        data_.first = e();
    }
    
    void set_second(object o)
    {
        extract<std::size_t> e(o);
        data_.second = e();
    }
    
    raw_pair data_;
    
    std::string str() const
    {
        std::ostringstream oss;
        oss << "(" << data_.first << "," << data_.second << ")";
        return oss.str();
    }
};

#endif
