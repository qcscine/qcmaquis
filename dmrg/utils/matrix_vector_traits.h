
#ifndef __ALPS_MATRIX_VECTOR_TRAITS_HPP__
#define __ALPS_MATRIX_VECTOR_TRAITS_HPP__

#include <complex>

namespace blas {
    
    template<class M>
    struct associated_diagonal_matrix { };
    
    template<class M>
    struct associated_vector { };
    
    template<class M>
    struct associated_real_diagonal_matrix { };
    
    template<class M>
    struct associated_real_vector { };
    
    namespace detail {
        template<class T> struct real_type { typedef T type; };
        template<class T> struct real_type<std::complex<T> > { typedef T type; };
    }
}

#endif
