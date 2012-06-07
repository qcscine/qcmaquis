namespace vli {
    namespace detail {

                     // declaration ASM functions  
                     __device__ void add384_384_gpu(unsigned int* x /* shared */, unsigned int const* y /* global */);
                     __device__ void mul384_384_gpu(unsigned int* x /* local */, unsigned int const* y /* shared */, unsigned int const* z /* shared */);
                    
                     // Algo order*order threads, based on full diagonals
                     template <typename BaseInt, std::size_t Size, unsigned int Order>
                     __device__ void diag_algo(unsigned int ThreadId, BaseInt const* a,  BaseInt const* b, BaseInt* c);

                } // end namespace detail
           } // end namespace vli
