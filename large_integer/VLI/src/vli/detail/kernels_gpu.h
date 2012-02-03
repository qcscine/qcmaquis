#ifndef VLI_KERNELS_GPU_H
#define VLI_KERNELS_GPU_H

#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include "vli/vli_config.h"

namespace vli {
namespace detail {


template <std::size_t N>
struct vli_size_tag
{
    enum {size = N};
};


#define VLI_DECLARE_GPU_FUNCTIONS(TYPE, VLI_SIZE) \
    void inner_product_vector(vli_size_tag<VLI_SIZE>, unsigned int max_order, std::size_t vector_size, TYPE const* A, TYPE const* B, TYPE* C ); \

#define VLI_DECLARE_GPU_FUNCTIONS_FOR(r, data, BASEINT_SIZE_PAIR) \
    VLI_DECLARE_GPU_FUNCTIONS( BOOST_PP_TUPLE_ELEM(2,0,BASEINT_SIZE_PAIR), BOOST_PP_TUPLE_ELEM(2,1,BASEINT_SIZE_PAIR) )

BOOST_PP_SEQ_FOR_EACH(VLI_DECLARE_GPU_FUNCTIONS_FOR, _, VLI_COMPILE_BASEINT_SIZE_PAIRS_SEQ)

#undef VLI_DECLARE_GPU_FUNCTIONS_FOR
#undef VLI_DECLARE_GPU_FUNCTIONS 


} //namespace detail

} //namespace vli

#endif
