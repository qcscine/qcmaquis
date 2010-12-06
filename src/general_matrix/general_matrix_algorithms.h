#ifndef __ALPS_GENERAL_MATRIX_ALGORITHMS_HPP__
#define __ALPS_GENERAL_MATRIX_ALGORITHMS_HPP__

#include "matrix_concept_check.hpp"
#include "diagonal_matrix.h"

#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

namespace blas
{
    namespace detail {
        template<typename T> struct sv_type { typedef T type; };
        template<typename T>
        struct sv_type<std::complex<T> > { typedef T type; };
    }
    
    template<typename T, class MemoryBlock>
    void svd(general_matrix<T, MemoryBlock> & M,
             general_matrix<T, MemoryBlock> & U,
             general_matrix<T, MemoryBlock>& V,
             typename associated_diagonal_matrix<general_matrix<T, MemoryBlock> >::type & S)
    {
        BOOST_CONCEPT_ASSERT((blas::Matrix<general_matrix<T, MemoryBlock> >));
        typename general_matrix<T, MemoryBlock>::size_type k = std::min(num_rows(M), num_columns(M));
        resize(U, num_rows(M), k);
        resize(V, k, num_columns(M));
        
        std::vector<typename detail::sv_type<T>::type> S_(k);
        boost::numeric::bindings::lapack::gesdd('S', M, S_, U, V);
        
        S = typename associated_diagonal_matrix<general_matrix<T, MemoryBlock> >::type(S_);
    }
    
    template<typename T, class MemoryBlock>
    void qr(general_matrix<T, MemoryBlock> & M,
            general_matrix<T, MemoryBlock> & Q,
            general_matrix<T, MemoryBlock> & R)
    {
        /* implement thin QR decomposition, i.e. for a (m,n) matrix, where m >= n, the result should be
         Q: (m,n)
         R: (n,n) */
    }
    
    template<typename T, class MemoryBlock>
    general_matrix<T, MemoryBlock> conjugate(general_matrix<T, MemoryBlock> M)
    {
        M.inplace_conjugate();
        return M;
    }
    
} /* namespace blas */

#endif
