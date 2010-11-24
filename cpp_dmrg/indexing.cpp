#include <iostream>

#include "indexing.h"
#include "symmetry.h"

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
        std::cout << "Called!" << std::endl;
        cp.value = -cp.value;
        return cp;
    }
    
    Enum value;
    
    operator Enum() const { return value; }
};

template<class sgrp> struct charge_wrapped_as { };
template<> struct charge_wrapped_as<NullGroup> { typedef EnumToClass<NullGroup::charge> type; };
template<> struct charge_wrapped_as<Ztwo> { typedef EnumToClass<Ztwo::charge> type; };

template<class sgrp>
typename charge_wrapped_as<sgrp>::type wrapped_fuse(object o1, object o2)
{
    typedef typename charge_wrapped_as<sgrp>::type wt;
    typedef typename sgrp::charge charge;
    
    extract<wt> w1(o1), w2(o2);
    if (!(w1.check() && w2.check()))
        std::cerr << "Conversion error!" << std::endl;
    
    return sgrp::fuse(w1(), w2());
}

BOOST_PYTHON_MODULE(indexing) {
#define EXPORT_CHARGE_ENUM(charge_enum, name) \
class_<EnumToClass<charge_enum> >(name) \
.def(init<int>()) \
.def("__str__", &EnumToClass<charge_enum>::str) \
.def("__repr__", &EnumToClass<charge_enum>::str) \
.def("set", &EnumToClass<charge_enum>::set) \
.def(-self);
    EXPORT_CHARGE_ENUM(NullGroup::charge, "NG_charge");
    EXPORT_CHARGE_ENUM(Ztwo::charge, "Z2_charge");
#undef EXPORT_CHARGE_ENUM
    
#define EXPORT_SYMM_GROUP(sgrp, name) \
class_<sgrp>(name) \
.def("fuse", &wrapped_fuse<sgrp>) \
.staticmethod("fuse");
    EXPORT_SYMM_GROUP(NullGroup, "NG");
    EXPORT_SYMM_GROUP(Ztwo, "Z2");
#undef EXPORT_SYMM_GROUP
    
    /* what needs to be done:
     ii) Export pair to tuple
     iii) Export index
     
     Questions:
     i) Is there something for pair in Boost already?
     */
}
