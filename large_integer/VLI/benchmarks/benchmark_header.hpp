#ifndef VLI_BENCH_HEADER_HPP
#define VLI_BENCH_HEADER_HPP

#define BOOST_TEST_MODULE vli_cpu

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
 
#ifdef VLI_USE_GPU
#include "vli/detail/gpu/inner_product_gpu_booster.hpp"
#include "vli/utils/gpu_manager.h"
#include "vli/utils/gpu_manager.hpp"
#endif //VLI_USE_GPU

#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial.hpp"
#include "vli/vli_cpu.h"
#include "vli/vli_traits.hpp"
#include "utils/timings.h"

#include "utils/tools.h"

using vli::vli_cpu;
using vli::polynomial;
using vli::vector_polynomial;
#define SIZE_VEC 128


#define VLI_VARIABLE0 'x'
#define VLI_VARIABLE1 'y'
#define VLI_VARIABLE2 'z'
#define VLI_VARIABLE3 'w'

#define VLI_GET_VAR(z,n,unused) BOOST_PP_COMMA_IF(n) vli::var<BOOST_PP_CAT(VLI_VARIABLE,n)> 
#define VLI_GET_NULVAR(z,n,unused) BOOST_PP_COMMA_IF(BOOST_PP_SUB(4,n)) vli::no_variable
#define VLI_EXPEND_VAR(ARG) BOOST_PP_REPEAT(ARG,VLI_GET_VAR,~)  BOOST_PP_REPEAT(BOOST_PP_SUB(4,ARG),VLI_GET_NULVAR,~)


//typedef vli
typedef vli_cpu< unsigned long int, VLI_SIZE> vli_type_cpu;
//typedef poly max order each
typedef vli::polynomial< vli_type_cpu, vli::max_order_each<ORDER_POLY>, VLI_EXPEND_VAR(VARIABLE_POLY) > Polynomial_type_each;
//typedef poly max order combined
typedef vli::polynomial< vli_type_cpu, vli::max_order_combined<ORDER_POLY>, VLI_EXPEND_VAR(VARIABLE_POLY) > Polynomial_type_combined;


#endif // VLI_BENCH_HEADER_HPP
