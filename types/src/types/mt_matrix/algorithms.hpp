#ifndef MT_MATRIX_ALGO_ALGORITHMS_HPP__
#define MT_MATRIX_ALGO_ALGORITHMS_HPP__
#include <vector>
#include <stdexcept>

#include "types/mt_matrix/mt_matrix.h"
#include "types/dense_matrix/algorithms.hpp"

namespace maquis {
    namespace types {
        namespace algorithms {
            
            template<typename T, class DiagMatrix>
            void svd(mt_matrix<T> const & M,
                     mt_matrix<T> & U,
                     mt_matrix<T> & V,
                     DiagMatrix & S)
            {
                M.wait(); U.wait(); V.wait();
                svd(M.data_, U.data_, V.data_, S);
            }
            
            template<typename T>
            void qr(mt_matrix<T> const & M,
                    mt_matrix<T> & Q,
                    mt_matrix<T> & R)
            {
                M.wait(); Q.wait(); R.wait();
                qr(M.data_, Q.data_, R.data_);
            }
            
            template<typename T>
            mt_matrix<T> exp (mt_matrix<T> M, T const & alpha)
            {
                M.wait();
                mt_matrix<T> ret;
                ret.data_ = exp(M.data_, alpha);
                return ret;
            }
            
            template<typename T, class Generator>
            void generate(mt_matrix<T> & m, Generator g)
            {
                m.generate(g);
            }
            
            template<typename T>
            void heev(mt_matrix<T> M,
                      mt_matrix<T> & evecs,
                      typename associated_real_vector<dense_matrix<T> >::type & evals) 
            {
                M.wait(); evecs.wait();
                heev(M.data_, evecs.data_, evals);
            }
            
            template<typename T>
            void heev(mt_matrix<T> M,
                      mt_matrix<T> & evecs,
                      typename associated_diagonal_matrix<dense_matrix<T> >::type & evals)
            {
                M.wait(); evecs.wait();
                heev(M.data_, evecs.data_, evals);
            }
            
            template<typename T>
            void syev(mt_matrix<T> const & M,
                      mt_matrix<T> & vecs,
                      typename maquis::types::diagonal_matrix<double> & S)
            {
                M.wait();
                vecs.wait();
                syev(M.data_, vecs.data_, S);
            }
            
        } // end namespace algorithms
    } // end namspace types
} //end namespace maquis

#endif
