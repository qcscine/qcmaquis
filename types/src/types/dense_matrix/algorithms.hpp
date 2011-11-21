#ifndef __ALPS_DENSE_MATRIX_ALGO_ALGORITHMS_HPP__
#define __ALPS_DENSE_MATRIX_ALGO_ALGORITHMS_HPP__
#include <vector>
#include <stdexcept>
#include "types/dense_matrix/matrix_concept_check.hpp"
#include "types/dense_matrix/diagonal_matrix.h"
#include "utils/function_objects.h"
#include "types/dense_matrix/matrix_algorithms.hpp"


#include <boost/numeric/bindings/lapack/driver/gesvd.hpp>
#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>
#include <boost/numeric/bindings/lapack/driver/syevd.hpp>
#include <boost/numeric/bindings/lapack/driver/heevd.hpp>
#include <boost/numeric/bindings/lapack/computational/geqrf.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include <utils/timings.h>

#include <alps/numeric/real.hpp>
#include <alps/numeric/imag.hpp>

// forward declaration for nested specialization, be cautious of the namespace

namespace maquis {
    namespace types {
        template<class T>
        class diagonal_matrix; 
    }
}

namespace maquis {
    namespace types {
        namespace algorithms {
      
        namespace detail {
            template<typename T> struct sv_type { typedef T type; };
            template<typename T>
            struct sv_type<std::complex<T> > { typedef T type; };
        }

        template<typename T, class MemoryBlock, class DiagMatrix>
        void svd(dense_matrix<T, MemoryBlock> M,
                 dense_matrix<T, MemoryBlock> & U,
                 dense_matrix<T, MemoryBlock>& V,
                 DiagMatrix & S)
        {
            BOOST_CONCEPT_ASSERT((maquis::types::Matrix<dense_matrix<T, MemoryBlock> >));
            typename dense_matrix<T, MemoryBlock>::size_type k = std::min(num_rows(M), num_cols(M));
            resize(U, num_rows(M), k);
            resize(V, k, num_cols(M));
            typename associated_vector<dense_matrix<typename detail::sv_type<T>::type, MemoryBlock> >::type S_(k);
            int info = boost::numeric::bindings::lapack::gesvd('S', 'S', M, S_, U, V);
            if (info != 0)
                throw std::runtime_error("Error in SVD!");

            S = DiagMatrix(S_);
        }
        
        template<typename T, class MemoryBlock>
        void qr(dense_matrix<T, MemoryBlock> M,
                dense_matrix<T, MemoryBlock> & Q,
                dense_matrix<T, MemoryBlock> & R)
        {
            typename dense_matrix<T, MemoryBlock>::size_type k = std::min(num_rows(M), num_cols(M));
    
            typename associated_vector<dense_matrix<typename detail::sv_type<T>::type, MemoryBlock> >::type tau(k);
            
            int info = 0; //boost::numeric::bindings::lapack::geqrf(M, tau);
            if (info != 0)
                throw std::runtime_error("Error in geqrf");
            
            resize(Q, num_rows(M), k);
            resize(R, k, num_cols(M));
            
            // get R
            std::fill(elements(R).first, elements(R).second, 0);
            for (std::size_t c = 0; c < num_cols(M); ++c)
                for (std::size_t r = 0; r <= c; ++r)
                    R(r, c) = M(r, c);
            
            // get Q from householder reflections in M
            std::fill(elements(Q).first, elements(Q).second, 0);
            
        }
        
        template<typename T, class MemoryBlock>
        dense_matrix<T, MemoryBlock> exp (dense_matrix<T, MemoryBlock> M, T const & alpha=1)
        {
            dense_matrix<T, MemoryBlock> N, tmp;
            typename associated_real_vector<dense_matrix<T, MemoryBlock> >::type Sv(num_rows(M));
            
            heev(M, N, Sv);
            
            typename associated_diagonal_matrix<dense_matrix<T, MemoryBlock> >::type S(Sv);
            S = exp(alpha*S);
            gemm(N, S, tmp);
            gemm(tmp, maquis::types::conjugate(transpose(N)), M);
            
            return M;
        }
  
        template<typename T, class MemoryBlock, class Generator>
        void generate(dense_matrix<T, MemoryBlock>& m, Generator g)
        {
           std::generate(elements(m).first, elements(m).second, g);
        }
    
        template<typename T, class MemoryBlock>
        void heev(dense_matrix<T, MemoryBlock> M,
                  dense_matrix<T, MemoryBlock> & evecs,
                  typename associated_real_vector<dense_matrix<T, MemoryBlock> >::type & evals) 
        {
            assert(num_rows(M) == num_cols(M));
            assert(evals.size() == num_rows(M));
#ifndef NDEBUG
            using utils::conj;
            for (int i = 0; i < num_rows(M); ++i)
                for (int j = 0; j < num_cols(M); ++j)
                    assert( abs( M(i,j) - conj(M(j,i)) ) < 1e-10 );
#endif
            
            boost::numeric::bindings::lapack::heevd('V', M, evals);
            // to be consistent with the SVD, I reorder in decreasing order
            std::reverse(evals.begin(), evals.end());
            // and the same with the matrix
            evecs.resize(num_rows(M), num_cols(M));
            for (std::size_t c = 0; c < num_cols(M); ++c)
            		std::copy(column(M, c).first, column(M, c).second,
                          column(evecs, num_cols(M)-1-c).first);
        }
        
        template<typename T, class MemoryBlock>
        void heev(dense_matrix<T, MemoryBlock> M,
                  dense_matrix<T, MemoryBlock> & evecs,
                  typename associated_diagonal_matrix<dense_matrix<T, MemoryBlock> >::type & evals)
        {
            assert(num_rows(M) == num_cols(M));
            typename associated_real_vector<dense_matrix<T, MemoryBlock> >::type evals_(num_rows(M));
            heev(M, evecs, evals_);
            evals = typename associated_diagonal_matrix<dense_matrix<T, MemoryBlock> >::type(evals_);
        }

        template<typename T, class MemoryBlock, class ThirdArgument>
        void syev(dense_matrix<T, MemoryBlock> M,
                  dense_matrix<T, MemoryBlock> & evecs,
                  ThirdArgument & evals)
        {
            heev(M, evecs, evals);
        }
        /*
        * Some block_matrix algorithms necessitate nested specialization due to ambient scheduler
        * the algos are full rewritten or partly with subset specialization 
        * an alternative implementation is presented inside p_dense_matrix/algorithms/algorithms.hpp
        */
        } // end namespace algorithms
    } // end namspace types
} //end namespace maquis

#endif
