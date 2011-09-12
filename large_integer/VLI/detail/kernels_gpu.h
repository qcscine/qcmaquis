#ifndef VLI_KERNELS_GPU_H
#define VLI_KERNELS_GPU_H

#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include "vli_utils/vli_config.h"

namespace vli {
namespace detail {


template <std::size_t N>
struct vli_size_tag
{
    enum {size = N};
};


#define VLI_DECLARE_GPU_FUNCTIONS(TYPE, VLI_SIZE) \
    void negate(vli_size_tag<VLI_SIZE>, TYPE* A); \
    void plus_assign(vli_size_tag<VLI_SIZE>, TYPE* A, TYPE const* B); \
    void plus_assign(vli_size_tag<VLI_SIZE>, TYPE* A, TYPE B); \
    void minus_assign(vli_size_tag<VLI_SIZE>, TYPE* A, TYPE const* B); \
    void multiplies_assign(vli_size_tag<VLI_SIZE>, TYPE* A, int a); \
    void entrywise_multiplies(vli_size_tag<VLI_SIZE>, TYPE const* A, TYPE const*  B, TYPE* C); \
    void poly_mono_multiply(vli_size_tag<VLI_SIZE>, TYPE const* A, TYPE const* B, TYPE* C, std::size_t j_exp, std::size_t h_exp); \
    void plus_assign_poly_int(vli_size_tag<VLI_SIZE>, TYPE* A, int a); \
    void plus_assign_poly(vli_size_tag<VLI_SIZE>, TYPE* A, TYPE const* B); \
    void minus_assign_poly(vli_size_tag<VLI_SIZE>, TYPE* A, TYPE const* B); \
    void poly_poly_multiply(vli_size_tag<VLI_SIZE>, TYPE const* A, TYPE const* B, TYPE* C); \
    void inner_product_vector(vli_size_tag<VLI_SIZE>, TYPE const* A, TYPE const* B, TYPE* C, std::size_t vector_size, std::size_t threads_per_block); \
    void vector_reduction(vli_size_tag<VLI_SIZE>, TYPE const* A, TYPE* B, std::size_t vector_size); 

#define VLI_DECLARE_GPU_FUNCTIONS_FOR(r, data, BASEINT_SIZE_PAIR) \
    VLI_DECLARE_GPU_FUNCTIONS( BOOST_PP_TUPLE_ELEM(2,0,BASEINT_SIZE_PAIR), BOOST_PP_TUPLE_ELEM(2,1,BASEINT_SIZE_PAIR) )

//#define VLI_COMPILE_BASEINT_SIZE_PAIRS_SEQ ((unsigned int,3)) ((unsigned long int,4)) ((unsigned long int,8))

BOOST_PP_SEQ_FOR_EACH(VLI_DECLARE_GPU_FUNCTIONS_FOR, _, VLI_COMPILE_BASEINT_SIZE_PAIRS_SEQ)

#undef VLI_DECLARE_GPU_FUNCTIONS_FOR
#undef VLI_DECLARE_GPU_FUNCTIONS 


} //namespace detail

} //namespace vli

#endif
