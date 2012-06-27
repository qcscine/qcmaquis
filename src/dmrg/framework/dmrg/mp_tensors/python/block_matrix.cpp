#include <iostream>

#define PYTHON_EXPORTS

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/detail/alps_matrix.hpp"

#include <boost/python.hpp>

using namespace boost::python;

#include "dmrg/mp_tensors/python/wrappers.h"

BOOST_PYTHON_MODULE(block_matrix) {
#define EXPORT_CHARGE_ENUM(charge_enum, name) \
class_<EnumToClass<charge_enum> >(name) \
.def(init<int>()) \
.def("__str__", &EnumToClass<charge_enum>::str) \
.def("__repr__", &EnumToClass<charge_enum>::str) \
.def("set", &EnumToClass<charge_enum>::set) \
.def(-self);
    EXPORT_CHARGE_ENUM(TrivialGroup::charge, "NG_charge");
    EXPORT_CHARGE_ENUM(Ztwo::charge, "Z2_charge");
#undef EXPORT_CHARGE_ENUM
    
#define EXPORT_SYMM_GROUP(sgrp, name) \
class_<sgrp>(name) \
.def("fuse", &wrapped_fuse<sgrp>) \
.staticmethod("fuse");
    EXPORT_SYMM_GROUP(TrivialGroup, "NG");
    EXPORT_SYMM_GROUP(Ztwo, "Z2");
#undef EXPORT_SYMM_GROUP
    
#define EXPORT_BLOCK_PAIR(sgrp, name) \
class_<wrapped_pair<sgrp> >(name) \
.def(init<wrapped_pair<sgrp>::wrapped_charge, std::size_t>()) \
.add_property("first", &wrapped_pair<sgrp>::get_first, &wrapped_pair<sgrp>::set_first) \
.add_property("second", &wrapped_pair<sgrp>::get_second, &wrapped_pair<sgrp>::set_second) \
.def("__str__", &wrapped_pair<sgrp>::str) \
.def("__repr__", &wrapped_pair<sgrp>::str)
    EXPORT_BLOCK_PAIR(TrivialGroup, "NG_pair");
    EXPORT_BLOCK_PAIR(Ztwo, "Z2_pair");
#undef EXPORT_BLOCK_PAIR
    
#define EXPORT_INDEX(sgrp, name) \
class_<Index<sgrp> >(name) \
.def("insert", &Index<sgrp>::py_insert) \
.def("sizes", &Index<sgrp>::py_sizes) \
.def("charges", &Index<sgrp>::py_charges)
    EXPORT_INDEX(TrivialGroup, "NG_index");
    EXPORT_INDEX(Ztwo, "Z2_index");
#undef EXPORT_INDEX
    
#define EXPORT_BLOCK_MATRIX(matrix, sgrp, name) \
class_<block_matrix<matrix, sgrp> >(name) \
.def("left_basis", &block_matrix<matrix, sgrp>::left_basis) \
.def("right_basis", &block_matrix<matrix, sgrp>::right_basis)
    EXPORT_BLOCK_MATRIX(alps::numeric::matrix<double>, TrivialGroup, "NG_d_matrix");
    EXPORT_BLOCK_MATRIX(alps::numeric::matrix<double>, Ztwo, "Z2_d_matrix");
#undef EXPORT_BLOCK_MATRIX
}
