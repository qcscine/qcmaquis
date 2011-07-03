#ifndef VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP
#define VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP

#include "kernels_cpu.hpp"

namespace vli{

template<class Vli, int Order>
class polynomial;
   
template<class Vli>
class monomial;
    
    
template <class BaseInt, int Size, int Order>
void poly_multiply(polynomial<vli_cpu<BaseInt, Size>, Order> & result, 
                   polynomial<vli_cpu<BaseInt, Size>, Order> const & p1, 
                   polynomial<vli_cpu<BaseInt, Size>, Order> const & p2)
{
    std::size_t max_order = p1.max_order;
    for(std::size_t je1 = 0; je1 < max_order; ++je1)
        for(std::size_t he1 = 0; he1 < max_order; ++he1)
        {
            for(std::size_t je2 = 0; je2 < max_order - je1; ++je2)
                for(std::size_t he2 = 0; he2 < max_order - he1; ++he2)
                    result.coeffs[ (je1+je2)*max_order + he1+he2 ] += p1.coeffs[je1*max_order+he1] * p2.coeffs[je2*max_order+he2];
        }    
}
        
template <class BaseInt, int Size, int Order>
polynomial<vli_cpu<BaseInt, Size>, Order> operator * (polynomial<vli_cpu<BaseInt, Size>, Order>  const& p, monomial<vli_cpu<BaseInt, Size> > const& m)
{
    polynomial<vli_cpu<BaseInt, Size>, Order> r;
    std::size_t max_order = r.max_order;
    for(std::size_t je = 0; je < max_order-m.j_exp; ++je)
        for(std::size_t he = 0; he < max_order-m.h_exp; ++he)
            r(je+m.j_exp,he+m.h_exp) = p(je,he) * m.coeff;
    return r;
}
    

}// end namespace 

#endif