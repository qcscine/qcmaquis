#ifndef VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP
#define VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP

#include "vli/polynomial/algorithms_polynomial_cpu.hpp"
#include "vli/vli_cpu.h"

namespace vli{

template<class Vli, unsigned int Order>
class polynomial_cpu;
   
template<class Vli>
struct monomial;

template<class BaseInt, std::size_t Size>
class vli_cpu;
    
/** Algo based on triangle + diagnoal decomposition, nthreads maximum **/
template <class BaseInt, std::size_t Size, unsigned int Order>
void poly_multiply_block_algo(polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> & result, 
                        polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                        polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2) {
    // first PASS, half top right corner, 
    unsigned int n(0);
    for(unsigned int i=0;i< Order;++i){ // i will be a thread here, independence loop
        for(unsigned int j=0; j<=n; ++j){
            block_algo(j,n-j,result,p1,p2);
        }
        n++; // thread num
    }
    
    // second PASS, half bottom left corner, 
    n=1;
    for(unsigned int i=1; i<Order;++i){  // i will be a thread here, independence loop
        for(unsigned int j=n; j<Order; ++j){
            block_algo(j,Order-j+n-1,result,p1,p2);
        }
        n++; // thread num
    }
}
    
/** Algo based on diagonal decomposition, nthreads*nthreads maximum **/    
template <class BaseInt, std::size_t Size, unsigned int Order>
void poly_multiply_diag_algo(polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> & result, 
                             polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                             polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2) {   
    //first pass
    for(unsigned int i(0); i < Order*Order ; ++i)
        diagonal_up(i,result,p1,p2);
    //second pass    
    for(unsigned int i(0); i < Order*Order ; ++i)
        diagonal_down(Order*Order - i,result,p1,p2); 
        
}
    
    
    
        
}// end namespace 

#endif
