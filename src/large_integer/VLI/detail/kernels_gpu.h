#ifndef VLI_KERNELS_GPU_H
#define VLI_KERNELS_GPU_H

#include <boost/preprocessor/seq/for_each.hpp>
#include "detail/vli_size_param.hpp"

namespace vli {
namespace detail {

    /**
    * VLI declarations
    */

    /** vli_gpu -> -vli_gpu */
    void negate(type_vli::BaseInt* A); 
    
    /** vli_gpu += vli_gpu */
    void plus_assign(type_vli::BaseInt* A, type_vli::BaseInt const* B); 
 
    /** vli_gpu += lambda */
    void plus_assign(type_vli::BaseInt* A, type_vli::BaseInt a); 

    /** vli_gpu -= vli_gpu */
    void minus_assign(type_vli::BaseInt* A, type_vli::BaseInt const* B);  

    /** vli_gpu = vli_gpu * vli_gpu*/
    void entrywise_multiplies(type_vli::BaseInt const* A, type_vli::BaseInt const*  B, type_vli::BaseInt* C);  

    /** vli_gpu *= lambda */
    void multiplies_assign(type_vli::BaseInt* A, int a); 

    /**
    *  MONO VLI declaration
    */

    /** vli_poly_gpu *= vli_mono_gpu */
    void poly_mono_multiply(type_vli::BaseInt const* A, type_vli::BaseInt const* B, type_vli::BaseInt* C);

    /** vli_poly_gpu += vli_mono_gpu */
    // void plus_assign_poly_mono(type_vli::BaseInt* A, type_vli::BaseInt const* B);

    /**
    *  POLY VLI declaration
    */

    /** vli_poly_gpu += vli_poly_gpu */
    void plus_assign_poly(type_vli::BaseInt* A, type_vli::BaseInt const* B);  

    /** vli_poly_gpu -= vli_poly_gpu */
    void minus_assign_poly(type_vli::BaseInt* A, type_vli::BaseInt const* B);  

    /** vli_poly_gpu *= vli_poly_gpu */
    void poly_poly_multiply(type_vli::BaseInt const* A, type_vli::BaseInt const* B, type_vli::BaseInt* C);

    /**
    *  VECTOR VLI declaration
    *  by definition vectors are dynamics so no templated
    */
    
    /** vli_vector_gpu[i] * vli_vector_gpu[i] */
    void inner_product_vector(type_vli::BaseInt const* A, type_vli::BaseInt const* B, type_vli::BaseInt* C, std::size_t SizeVector, std::size_t NumThreads); 
    
    /** vli_poly_gpu += vli_vector_gpu[i] */
    void vector_reduction(type_vli::BaseInt const* A, type_vli::BaseInt* B, std::size_t SizeVector);



/*

#define VLI_GPU_BASE_INT_TYPES_SEQ (unsigned int) (unsigned long int)
#define VLI_DECLARE_GPU_KERNELS_FOR(r, data, TYPE) \
    void inner_product_gpu(TYPE const* A, TYPE const* B, TYPE* C, int num_integers, int vli_size);  \
    void inner_product_vector_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size, int max_order, TYPE vector_size); \
    void vector_reduction_gpu(TYPE const* A, TYPE * B,  int vli_size, int max_order, TYPE vector_size);
BOOST_PP_SEQ_FOR_EACH(VLI_DECLARE_GPU_KERNELS_FOR, _, VLI_GPU_BASE_INT_TYPES_SEQ)

#undef VLI_DECLARE_GPU_KERNELS_FOR
*/
} //namespace detail

} //namespace vli

#endif
