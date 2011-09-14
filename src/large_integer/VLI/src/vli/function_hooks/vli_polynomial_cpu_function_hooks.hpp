#ifndef VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP
#define VLI_POLYNOMIAL_CPU_FUNCTION_HOOKS_HPP

namespace vli{

template<class Vli, unsigned int Order>
class polynomial_cpu;
   
template<class Vli>
struct monomial;

template<class BaseInt, std::size_t Size>
class vli_cpu;
    
    
template <class BaseInt, std::size_t Size, unsigned int Order>
void poly_multiply(polynomial_cpu<vli_cpu<BaseInt, Size>, Order> & result, 
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p1, 
                   polynomial_cpu<vli_cpu<BaseInt, Size>, Order> const & p2)
{
    typedef typename polynomial_cpu<vli_cpu<BaseInt,Size>,Order>::exponent_type exponent_type;
 
    for(exponent_type je1 = 0; je1 < Order; ++je1)
    {
        for(exponent_type he1 = 0; he1 < Order; ++he1)
        {
            for(exponent_type je2 = 0; je2 < Order - je1; ++je2)
            {
                for(exponent_type he2 = 0; he2 < Order - he1; ++he2)
                {
                     result.coeffs_[(je1+je2)*Order + he1+he2 ] += p1.coeffs_[je1*Order+he1] * p2.coeffs_[je2*Order+he2];
                }
            }
        }    
    }
}
        
}// end namespace 

#endif
