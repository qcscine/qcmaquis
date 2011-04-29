#ifndef VLI_KERNELS_GPU_H
#define VLI_KERNELS_GPU_H

namespace vli {
namespace detail {

void addition_gpu(int* A, const int*  B, int num_integer, int ld);  

void multiply_gpu(const int * A, const int*  B, int* C, int num_integer, int ld);  



} //namespace detail

} //namespace vli



#endif
