#ifndef VLI_KERNELS_GPU_H
#define VLI_KERNELS_GPU_H

namespace vli {
namespace detail {

void addition_gpu(int* A, const int*  B, int num_integers, int vli_size);  

void multiply_gpu(const int * A, const int*  B, int* C, int num_integers, int vli_size);  

void inner_prod_gpu(const int* A, const int* B, int* C, int num_integers, int vli_size);


} //namespace detail

} //namespace vli



#endif
