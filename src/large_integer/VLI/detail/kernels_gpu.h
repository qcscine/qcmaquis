#ifndef VLI_KERNELS_GPU_H
#define VLI_KERNELS_GPU_H

#include <boost/preprocessor/seq/for_each.hpp>

namespace vli {
namespace detail {


#define VLI_GPU_BASE_INT_TYPES_SEQ (unsigned int) (unsigned long int)

#define VLI_DECLARE_GPU_KERNELS_FOR(r, data, TYPE) \
    void negate_gpu(TYPE* A, int vli_size); \
    void plus_assign_gpu(TYPE* A, TYPE const* B, int num_integers, int vli_size); \
    void minus_assign_gpu(TYPE* A, TYPE* B, int vli_size);  \
    void entrywise_multiplies_gpu(TYPE const* A, TYPE const*  B, TYPE* C, int vli_size);  \
    void inner_product_gpu(TYPE const* A, TYPE const* B, TYPE* C, int num_integers, int vli_size);  \
    void poly_multiply_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size, int max_order);  \
    void poly_addition_gpu(TYPE* A, TYPE const* B, int vli_size, int max_order);  \
    void poly_substraction_gpu(TYPE* A, TYPE const* B, int vli_size, int max_order);  \
    void poly_mono_multiply_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size, int max_order);  \
    void inner_product_vector_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size, int max_order, int vector_size); \
    void vector_reduction_gpu(TYPE const* A, TYPE * B,  int vli_size, int max_order, int vector_size);
    
BOOST_PP_SEQ_FOR_EACH(VLI_DECLARE_GPU_KERNELS_FOR, _, VLI_GPU_BASE_INT_TYPES_SEQ)

#undef VLI_DECLARE_GPU_KERNELS_FOR

} //namespace detail

} //namespace vli

#endif
