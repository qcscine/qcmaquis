#include <iostream>

#define PYTHON_EXPORTS

#include "indexing.h"
#include "symmetry.h"

#include <boost/python.hpp>

using namespace boost::python;

#include "indexing_wrappers.h"

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
    
#define EXPORT_BLOCK_PAIR(sgrp, name) \
class_<wrapped_pair<sgrp> >(name) \
.def(init<wrapped_pair<sgrp>::wrapped_charge, std::size_t>()) \
.add_property("first", &wrapped_pair<sgrp>::get_first, &wrapped_pair<sgrp>::set_first) \
.add_property("second", &wrapped_pair<sgrp>::get_second, &wrapped_pair<sgrp>::set_second) \
.def("__str__", &wrapped_pair<sgrp>::str) \
.def("__repr__", &wrapped_pair<sgrp>::str)
    EXPORT_BLOCK_PAIR(NullGroup, "NG_pair");
    EXPORT_BLOCK_PAIR(Ztwo, "Z2_pair");
#undef EXPORT_BLOCK_PAIR
    
#define EXPORT_INDEX(sgrp, name) \
class_<Index<sgrp> >(name) \
.def("insert", &Index<sgrp>::py_insert)
    EXPORT_INDEX(NullGroup, "NG_index");
#undef EXPORT_INDEX
    
    /* what needs to be done:
     iii) Export index
      * Element access
      * sizes(), charges()
     */
}
