#ifndef VLI_KERNELS_GPU_H
#define VLI_KERNELS_GPU_H

/** TO DO BOOST PREPROC FOR GENERATING ALL THIS LISTING FOR INT AND LONG INT**/
namespace vli {
namespace detail {

void plus_assign_gpu(int* A, int const* B, int num_integers, int vli_size);
//void addition_gpu(int* A, const int*  B, int num_integers, int vli_size);  

void entrywise_multiplies_assign_gpu(int* A, int const* B, int num_integers, int vli_size);

void entrywise_multiplies_gpu(int const* A, int const*  B, int* C, int num_integers, int vli_size);

//void multiply_gpu(const int * A, const int*  B, int* C, int num_integers, int vli_size);  

void inner_product_gpu(int const* A, int const* B, int* C, int num_integers, int vli_size);
/** multiplication polynome **/
void poly_multiply_gpu(int const* A, int const* B, int* C, int vli_size, int max_order);
/** addition polynome **/
void poly_addition_gpu(int* A, int const* B, int vli_size, int max_order);
/** multiplication polynome-monome **/
void poly_mono_multiply_gpu(int const* A, int const* B, int* C, int vli_size, int max_order);
/** equality between to gpu buffer e.g. VLI, mono, poly, dut to latency is it usefull **/
void equal_gpu(int const* A, int const* B, int vli_size, int* T);    
    
    
} //namespace detail

} //namespace vli



#endif
