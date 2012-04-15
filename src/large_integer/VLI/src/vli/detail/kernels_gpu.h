#ifndef VLI_KERNELS_GPU_H
#define VLI_KERNELS_GPU_H

#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include "vli/vli_config.h"

namespace vli {
namespace detail {


template <std::size_t Size, unsigned int Order>
struct vli_size_tag
{
    enum {size = Size};
    enum {order = Order};
};


#define VLI_DECLARE_GPU_FUNCTIONS(TYPE, VLI_SIZE, POLY_ORDER) \
    void inner_product_vector(vli_size_tag<VLI_SIZE, POLY_ORDER>, std::size_t vector_size, TYPE const* A, TYPE const* B, TYPE* C ); \
    void addition_gpu(vli_size_tag<2*VLI_SIZE,POLY_ORDER>, TYPE* A, TYPE const* B); \
    void multiplication_gpu(vli_size_tag<VLI_SIZE,POLY_ORDER>, TYPE* A, TYPE const* B, TYPE const* C);

#define VLI_DECLARE_GPU_FUNCTIONS_FOR(r, data, BASEINT_SIZE_ORDER_TUPLE) \
    VLI_DECLARE_GPU_FUNCTIONS( BOOST_PP_TUPLE_ELEM(3,0,BASEINT_SIZE_ORDER_TUPLE), BOOST_PP_TUPLE_ELEM(3,1,BASEINT_SIZE_ORDER_TUPLE), BOOST_PP_TUPLE_ELEM(3,2,BASEINT_SIZE_ORDER_TUPLE) )

BOOST_PP_SEQ_FOR_EACH(VLI_DECLARE_GPU_FUNCTIONS_FOR, _, VLI_COMPILE_BASEINT_SIZE_ORDER_TUPLE_SEQ)

#undef VLI_DECLARE_GPU_FUNCTIONS_FOR
#undef VLI_DECLARE_GPU_FUNCTIONS 


} //namespace detail

} //namespace vli

#endif
