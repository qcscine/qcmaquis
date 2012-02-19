#ifndef VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP
#define VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP

#include "vli/polynomial/algorithms_polynomial_cpu.hpp"
#include "vli/vli_cpu.hpp"

namespace vli{

template<class Vli, unsigned int Order>
class polynomial_cpu;
   
template<class Vli>
struct monomial;

template<class BaseInt, std::size_t Size>
class vli_cpu;
    
    
template <class BaseInt, std::size_t Size, unsigned int Order>
void poly_multiply(polynomial_cpu<vli_cpu<BaseInt, 2*Size>, 2*Order> & result, 
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{
    typedef typename polynomial_cpu<vli_cpu<BaseInt,Size>,Order>::exponent_type exponent_type;
 
    // C - we can not still use the operator *, += still ok 
//    vli_cpu<BaseInt, 2*Size> inter;
    
    for(exponent_type je1 = 0; je1 < Order; ++je1)
    {
        for(exponent_type je2 = 0; je2 < Order; ++je2)
        {
            for(exponent_type he1 = 0; he1 < Order; ++he1)
            {
                for(exponent_type he2 = 0; he2 < Order; ++he2)
                {
                    multi_nt(result.coeffs_[(je1+je2)*2*Order + he1+he2 ],p1.coeffs_[je1*Order+he1],p2.coeffs_[je2*Order+he2]);
//                    result.coeffs_[(je1+je2)*2*Order + he1+he2 ] += inter;
//                    memset((void*)&inter[0],0,2*Size*sizeof(BaseInt));
                    // original version
                    // result.coeffs_[(je1+je2)*2*Order + he1+he2 ] += p1.coeffs_[je1*Order+he1] * p2.coeffs_[je2*Order+he2];
                }
            }
        }    
    }
}

/** Algo based on triangle + diagnoal decomposition, nthreads maximum **/
template <class BaseInt, std::size_t Size, unsigned int Order>
void poly_multiply_block_algo(polynomial_cpu<vli_cpu<BaseInt, Size>, 2*Order> & result, 
                        polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                        polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{
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
                             polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{   
    //first pass
    for(unsigned int i(0); i < Order*Order ; ++i){
        diagonal_up(i,result,p1,p2);
    }
    //second pass    
    for(unsigned int i(0); i < Order*Order ; ++i){
        diagonal_down(Order*Order - i,result,p1,p2); 
    }    
}
    
    
    
        
}// end namespace 

#endif
