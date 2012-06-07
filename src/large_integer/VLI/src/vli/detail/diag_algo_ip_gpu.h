
//#include "bit_masks.hpp"


namespace vli {
    namespace detail {

                 template <typename BaseInt, std::size_t Size, unsigned int Order>
                 void diag_algo(unsigned int threadid, BaseInt const* a, BaseInt const* b, BaseInt *c)
                 {
                     int qa,ra,qb,rb,pos; // find all indexes
                     int j = Order*Order-1;

                     BaseInt sc[2*Size];
                     __shared__ BaseInt sa[Size*Order*Order];
                     __shared__ BaseInt sb[Size*Order*Order];
                     __shared__ BaseInt scc[2*Size*2*Order*2*Order];
                
                     #pragma unroll
                     for(std::size_t k=0 ; k < Size ;++k){
                         sa[threadid+Order*Order*k] = a[threadid+Order*Order*k];
                         sb[threadid+Order*Order*k] = b[threadid+Order*Order*k];
                     }

                     #pragma unroll
                     for(std::size_t k=0 ; k < 8*Size ;++k)
                         scc[threadid+Order*Order*k] = 0;
                
                     __syncthreads(); // we sync to be sure sa, sb and sc are loaded fully
                     
                     
                     for(int i(0); i < Order*Order; ++i)
                     {
                         #pragma unroll
                         for(int k=0;k<2*Size;++k)
                             sc[k] = 0;                      

                         qa = i/Order;
                         ra = i%Order;
                         qb = ((i <= threadid) ? (threadid - i) : (Order*Order-1) - i) / Order;
                         rb = ((i <= threadid) ? (threadid - i) : (Order*Order-1) - i) % Order;
                         int offset = ((i <= threadid) ? (thradid - i) : (Order*Order-1) - i ) * Size; 

                         pos = 2*(qa+qb)*Order + (ra+rb);
                         mul384_384_gpu(&sc[0],&sa[offset],&sb[i*Size]);
                         add384_384_gpu(&scc[2*Size*pos],&sc[0]);
                     }
                     for(int i
                     __syncthreads(); // we sync to be sure sa, sb and sc are loaded fully

                     #pragma unroll
                     for(std::size_t k=0 ; k < 8*Size ;++k)
                         c[threadid+Order*Order*k] = scc[threadid+Order*Order*k];
             
                    }
         } // end namespace detail
    } // end namespace vli
