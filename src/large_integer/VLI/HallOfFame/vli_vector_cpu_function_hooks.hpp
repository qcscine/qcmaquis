#ifndef  VLI_VECTOR_CPU_FUNCTION_HOOKS_HPP 
#define  VLI_VECTOR_CPU_FUNCTION_HOOKS_HPP

#include "kernels_cpu.hpp"

namespace vli {

    /**
    template forward declaration 
    */
    template <class Vli, int Size>
    class vli_vector_cpu;

namespace detail {

template <class Vli, int Size>
void plus_assign(vli_vector_cpu<Vli, Size> & v_a, vli_vector_cpu<Vli, Size> const& v_b )
{
    assert( v_a.size() == v_b.size() );

    std::size_t size = v_a.size();    
    #pragma omp parallel for private(i)
    for(size_type i = 0; i < size; ++i)
        addition_classic_cpu(v_a[i], v_b[i]);
}

template <class Vli, int Size>
void multiplies_assign(vli_vector_cpu<Vli, Size> & v_a, Vli const& vli )
{
    typedef typename VliType::size_type size_type;

    #pragma omp parallel for private(i)
    for(size_type i = 0; i < v_a.size(); ++i)
        multiplication_classic_cpu( v_a[i], vli);
}

template <class VliType>
void entrywise_multiplies_assign(vli_vector_cpu<VliType> & v_a, vli_vector_cpu<VliType> const& v_b )
{
    typedef typename VliType::size_type size_type;

    assert( v_a.size() == v_b.size() );

    #pragma omp parallel for private(i)
    for(size_type i = 0; i < v_b.size(); ++i)
        multiplication_classic_cpu( v_a[i], v_b[i]);
}
    
template<class Vli, int Order>

    
    
    
    
} //namespace detail
} //namespace vli

#endif // VLI_VECTOR_CPU_FUNCTION_HOOKS_HPP
