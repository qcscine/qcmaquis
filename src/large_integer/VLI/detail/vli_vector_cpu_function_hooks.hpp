#ifndef  VLI_VECTOR_CPU_FUNCTION_HOOKS_HPP 
#define  VLI_VECTOR_CPU_FUNCTION_HOOKS_HPP

#include "kernels_cpu.hpp"

namespace vli {

    /**
    template forward declaration 
    */
    template <class VliType>
    class vli_vector;

namespace detail {

template <class VliType>
void plus_assign(vli_vector<VliType> & v_a, vli_vector<VliType> const& v_b )
{
    assert( v_a.size() == v_b.size() );
    
    #pragma omp parallel for private(i)
    for(std::size_t i = 0; i < v_b.size(); ++i)
        addition_classic_cpu( v_a[i], v_b[i]);
}

template <class VliType>
void multiplies_assign(vli_vector<VliType> & v_a, VliType const& vli )
{
    #pragma omp parallel for private(i)
    for(std::size_t i = 0; i < v_a.size(); ++i)
        multiplication_classic_cpu( v_a[i], vli);
}

template <class VliType>
void entrywise_multiplies_assign(vli_vector<VliType> & v_a, vli_vector<VliType> const& v_b )
{
    assert( v_a.size() == v_b.size() );

    #pragma omp parallel for private(i)
    for(std::size_t i = 0; i < v_b.size(); ++i)
        multiplication_classic_cpu( v_a[i], v_b[i]);
}

} //namespace detail
} //namespace vli

#endif // VLI_VECTOR_CPU_FUNCTION_HOOKS_HPP
