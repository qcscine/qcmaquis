#ifndef __MAQUIS_MATRIX_ALGORITHMS_HPP__
#define __MAQUIS_MATRIX_ALGORITHMS_HPP__

#include "types/utils/traits.hpp"

namespace maquis {
    namespace types {
        
        
        /* NUMERIC MATRIX CONCEPT */

        namespace detail {
            template<typename T> struct sv_type { typedef T type; };
            template<typename T>
            struct sv_type<std::complex<T> > { typedef T type; };
        }

        template <class Matrix>
        void svd(Matrix M, Matrix & U, Matrix& V,
                 typename associated_real_diagonal_matrix<Matrix>::type& S);

        template <class Matrix>
        void qr(Matrix M, Matrix & Q, Matrix & R);
        
        template <class Matrix>
        Matrix exp (Matrix M, typename Matrix::value_type const & alpha=1);
  
        template <class Matrix, class Generator>
        void generate(Matrix& m, Generator g);
    
        template <class Matrix>
        void heev(Matrix M, Matrix & evecs,
                  typename associated_real_vector<Matrix>::type & evals);
        
        template <class Matrix>
        void heev(Matrix M, Matrix & evecs,
                  typename associated_diagonal_matrix<Matrix>::type & evals);

        template<class Matrix, class ThirdArgument>
        void syev(Matrix M, Matrix & evecs, ThirdArgument & evals);
        
        
        /* DMRG MATRIX CONCEPT */
        
        template <class Matrix>
        void reshape_r2l(Matrix& left, const Matrix& right,
                         size_t left_offset, size_t right_offset, 
                         size_t sdim, size_t ldim, size_t rdim);
        
        template <class Matrix>
        void reshape_l2r(const Matrix& left, Matrix& right,
                         size_t left_offset, size_t right_offset, 
                         size_t sdim, size_t ldim, size_t rdim);
        
        template <class Matrix>
        void lb_tensor_mpo(Matrix& out, const Matrix& in, const Matrix& alfa,
                           size_t out_offset, size_t in_offset, 
                           size_t sdim1, size_t sdim2, size_t ldim, size_t rdim);
        
        template <class Matrix>
        void rb_tensor_mpo(Matrix& out, const Matrix& in, const Matrix& alfa,
                           size_t out_offset, size_t in_offset, 
                           size_t sdim1, size_t sdim2, size_t ldim, size_t rdim);
         
        template <class Matrix>
        void scalar_norm(const Matrix& M, typename Matrix::value_type& ret);
        
        template <class Matrix>
        void scalar_norm(Matrix & M1, Matrix & M2, typename Matrix::value_type & ret);// not const due to nullcut        
        template <class DiagMatrix>
        void bond_renyi_entropies(const DiagMatrix & M, typename associated_real_vector<DiagMatrix>::type& sv);
        
        template <class Matrix>
        void left_right_boundary_init(Matrix & M);
        
    } // end namspace types
} //end namespace maquis

#endif
