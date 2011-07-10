#ifndef VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP
#define VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP

#include "detail/kernels_cpu.hpp"

namespace vli{

template<class Vli, int Order>
class polynomial_cpu;
   
template<class Vli>
class monomial;
    
    
template <class BaseInt, int Size, int Order>
void poly_multiply(polynomial_cpu<vli_cpu<BaseInt, Size>, Order> & result, 
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{
    typedef typename polynomial_cpu<vli_cpu<BaseInt,Size>,Order>::size_type size_type;
 
    for(size_type je1 = 0; je1 < Order; ++je1)
        for(size_type he1 = 0; he1 < Order; ++he1)
        {
            for(size_type je2 = 0; je2 < Order - je1; ++je2)
                for(size_type he2 = 0; he2 < Order - he1; ++he2)
                    result.coeffs[ (je1+je2)*Order + he1+he2 ] += p1.coeffs[je1*Order+he1] * p2.coeffs[je2*Order+he2];
        }    
}
        
template <class BaseInt, int Size, int Order>
polynomial_cpu<vli_cpu<BaseInt, Size>, Order> operator * (polynomial_cpu<vli_cpu<BaseInt, Size>, Order>  const& p, monomial<vli_cpu<BaseInt, Size> > const& m)
{
    typedef typename polynomial_cpu<vli_cpu<BaseInt,Size>,Order>::size_type size_type;
    
    polynomial_cpu<vli_cpu<BaseInt, Size>, Order> r;
    // TODO perhaps there is a better way to write these loops,
    //      so that they can be unrolled.
    for(size_type je = 0; je < Order-m.j_exp; ++je)
        for(size_type he = 0; he < Order-m.h_exp; ++he)
            r(je+m.j_exp,he+m.h_exp) = p(je,he) * m.coeff;
    return r;
}
    

}// end namespace 

#endif
