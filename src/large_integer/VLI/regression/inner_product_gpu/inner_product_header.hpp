#ifndef VLI_INNERPRODUCT_HEADER_HPP
#define VLI_INNERPRODUCT_HEADER_HPP

#define BOOST_TEST_MODULE innerproduct 

#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
 
#ifdef VLI_USE_GPU
#include "vli/detail/gpu/inner_product_gpu_accelerator.hpp"
#endif //VLI_USE_GPU

#include "vli/polynomial/vector_polynomial_cpu.hpp"
#include "vli/polynomial/polynomial.hpp"
#include "vli/vli.hpp"
#include "utils/timings.h"

#include "utils/tools.h"


#define Size_vec 8
#define Order 10

using vli::polynomial;
using vli::vector_polynomial;
//typedef vli
typedef vli::vli<64*VLI_SIZE> integer_type_cpu;
//typedef poly max order each
typedef vli::polynomial< integer_type_cpu, vli::max_order_each<Order>, vli::var<'x'>  >polynomial_type_each_x;
typedef vli::polynomial< integer_type_cpu, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>  >polynomial_type_each_xy;
typedef vli::polynomial< integer_type_cpu, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>  >polynomial_type_each_xyz;
typedef vli::polynomial< integer_type_cpu, vli::max_order_each<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'>  >polynomial_type_each_xyzw;
//typedef poly max order combined
typedef vli::polynomial< integer_type_cpu, vli::max_order_combined<Order>, vli::var<'x'> > polynomial_type_combined_x;
typedef vli::polynomial< integer_type_cpu, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'> > polynomial_type_combined_xy;
typedef vli::polynomial< integer_type_cpu, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'> > polynomial_type_combined_xyz;
typedef vli::polynomial< integer_type_cpu, vli::max_order_combined<Order>, vli::var<'x'>, vli::var<'y'>, vli::var<'z'>, vli::var<'w'> > polynomial_type_combined_xyzw;

//polynomial_type_each_xyzw, // buffer too large
typedef boost::mpl::list<  polynomial_type_each_x,
                           polynomial_type_each_xy,
                           polynomial_type_each_xyz
                           > polynomial_list_max_order_each;

typedef boost::mpl::list<  polynomial_type_combined_x,
                           polynomial_type_combined_xy,
                           polynomial_type_combined_xyz,
                           polynomial_type_combined_xyzw
                          > polynomial_list_max_order_combined;
#endif 
