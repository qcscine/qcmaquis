#ifndef VLI_KERNELS_GPU_H
#define VLI_KERNELS_GPU_H

/** TO DO BOOST PREPROC FOR GENERATING ALL THIS LISTING FOR INT AND LONG INT**/
namespace vli {
namespace detail {

typedef unsigned int TYPE;

void plus_assign_gpu(TYPE* A, TYPE const* B, int num_integers, int vli_size);
    
void minus_assign_gpu(TYPE* A, TYPE* B, int vli_size);
    
void entrywise_multiplies_assign_gpu(TYPE* A, TYPE const* B, int num_integers, int vli_size);

void entrywise_multiplies_gpu(TYPE const* A, TYPE const*  B, TYPE* C, int num_integers, int vli_size);

//void multiply_gpu(const int * A, const int*  B, int* C, int num_integers, int vli_size);  

void inner_product_gpu(TYPE const* A, TYPE const* B, TYPE* C, int num_integers, int vli_size);
/** multiplication polynome **/
void poly_multiply_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size, int max_order);
/** addition polynome **/
void poly_addition_gpu(TYPE* A, TYPE const* B, int vli_size, int max_order);
/** substraction polynome **/
void poly_substraction_gpu(TYPE* A, TYPE const* B, int vli_size, int max_order);    
/** multiplication polynome-monome **/
void poly_mono_multiply_gpu(TYPE const* A, TYPE const* B, TYPE* C, int vli_size, int max_order);
/** inner product of vector of polynome **/
void inner_product_vector_gpu(TYPE const* A, TYPE const* B, TYPE* C, TYPE* D, int vli_size, int max_order, int vector_size);
    
} //namespace detail

} //namespace vli



#endif
