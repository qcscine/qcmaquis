#ifndef AMBIENT_NUMERIC_MATRIX_ALGORITHMS_EXPERIMENTAL
#define AMBIENT_NUMERIC_MATRIX_ALGORITHMS_EXPERIMENTAL

#include "ambient/numeric/matrix/matrix.h"
#include "ambient/numeric/kernels/kernels.hpp"
#include "ambient/numeric/kernels/experimental.hpp"
#include "ambient/numeric/matrix/diagonal_matrix.hpp"

#define size_type       typename matrix<T,A>::size_type
#define real_type       typename matrix<T,A>::real_type
#define scalar_type     typename matrix<T,A>::scalar_type
#define difference_type typename matrix<T,A>::difference_type

namespace ambient { namespace numeric {

    template<int alfa, int beta, class MatrixViewA, class MatrixViewB, class MatrixViewC>
    inline void gemv(const MatrixViewA& a, size_t aoffset, 
                     const MatrixViewB& b, size_t boffset, 
                           MatrixViewC& c, size_t coffset, 
                           size_t m, size_t n)
    {
        kernels::gemv<alfa,beta,MatrixViewA,MatrixViewB,MatrixViewC>::spawn<complexity::N2>(a, aoffset, 
                                                                                            b, boffset, 
                                                                                            c, coffset, 
                                                                                            m, n); 
    }

    template<int ADD, class MA, class MB, class MC, class MF>
    inline void gemv_scale(const MA& a, size_t aoffset, 
                           const MB& b, size_t boffset, 
                                 MC& c, size_t coffset, 
                           const MF& f, size_t foffset, 
                           size_t m, size_t n)
    {
        kernels::gemv_scale<ADD,MA,MB,MC,MF>::spawn<complexity::N2>(a, aoffset, 
                                                                    b, boffset, 
                                                                    c, coffset, 
                                                                    f, foffset, 
                                                                    m, n); 
    }

    template<typename T, class A>
    inline void gbbrd(matrix<T,A>& a, diagonal_matrix<T>& d, diagonal_matrix<T>& e, matrix<T,A>& u, matrix<T,A>& v){
        kernels::gbbrd<T>::spawn<complexity::N3>(a, d, e, u, v);
    }

    template<typename T, class A>
    inline void gebrd(matrix<T,A>& a, diagonal_matrix<T>& d, diagonal_matrix<T>& e, matrix<T,A>& u, matrix<T,A>& v){
        kernels::gebrd<T>::spawn<complexity::N3>(a, d, e, u, v);
    }

    template<typename T, class A>
    inline void bdsqr(diagonal_matrix<T>& d, diagonal_matrix<T>& e, matrix<T,A>& u, matrix<T,A>& v){
        kernels::bdsqr<T>::spawn<complexity::N3>(d, e, u, v);
    }

    template<typename T, class A>
    inline void gebd2(matrix<T,A>& a, diagonal_matrix<T>& d, diagonal_matrix<T>& e, diagonal_matrix<T>& tq, diagonal_matrix<T>& tp){
        kernels::gebd2<T>::spawn<complexity::N3>(a, d, e, tq, tp);
    }

    template<PLASMA_enum TR, typename T, class A>
    inline void larfg(matrix<T,A>& a, diagonal_matrix<T>& t, diagonal_matrix<T>& d, size_t k){
        kernels::larfg<T,TR>::spawn<complexity::N3>(a, t, d, k);
    }

    template<typename T, class A>
    inline void labrd_update_col(matrix<T,A>& say, const matrix<T,A>& sax, 
                           matrix<T,A>& sy, const matrix<T,A>& sx, 
                           diagonal_matrix<T>& tq, 
                           diagonal_matrix<T>& d, 
                           int i)
    {
        kernels::labrd_update_col<T>::spawn<complexity::N3>(say, sax, sy, sx, tq, d, i);
    }

    template<typename T, class A>
    inline void labrd_reduce_col(matrix<T,A>& say, const matrix<T,A>& sax, 
                           matrix<T,A>& sy, const matrix<T,A>& sx, 
                           int i)
    {
        kernels::labrd_reduce_col<T>::spawn<complexity::N3>(say, sax, sy, sx, i);
    }

    template<typename T, class A>
    inline void labrd_update_row(const matrix<T,A>& say, matrix<T,A>& sax, 
                           const matrix<T,A>& sy, matrix<T,A>& sx, 
                           diagonal_matrix<T>& tp, 
                           diagonal_matrix<T>& e, 
                           int i)
    {
        kernels::labrd_update_row<T>::spawn<complexity::N3>(say, sax, sy, sx, tp, e, i);
    }

    template<typename T, class A>
    inline void labrd_reduce_row(const matrix<T,A>& say, matrix<T,A>& sax, 
                           const matrix<T,A>& sy, matrix<T,A>& sx, 
                           int i)
    {
        kernels::labrd_reduce_row<T>::spawn<complexity::N3>(say, sax, sy, sx, i);
    }

    template<PLASMA_enum UL, size_t OFF, typename T, class A>
    inline void laset2(matrix<T,A>& a, const T& alfa = 0.0){
        kernels::laset2<T,UL,OFF>::spawn<complexity::N2>(a, alfa);
    }

    template<PLASMA_enum UL, typename T, class A>
    inline void copy_band(const matrix<T,A>& src, matrix<T,A>& dst, size_t dj){
        kernels::copy_band<T,UL>::spawn<complexity::N2>(src, dst, dj);
    }

    template <int alfa, typename T, class A>
    inline void add_vectors(matrix<T,A>& lhs, size_t loffset, const matrix<T,A>& rhs, size_t roffset, size_t size){ 
        kernels::add_vectors<alfa, T>::spawn<complexity::N2>(lhs, loffset, rhs, roffset, size); 
    }

    template<typename T, class A>
    inline void sqrt_inplace(matrix<T,A>& a){
        kernels::template sqrt_inplace<T>::template spawn<complexity::N>(a);
    }

    template<typename T, class A>
    inline void norm_vector(const matrix<T,A>& a, matrix<T,A>& b){ 
        kernels::template norm_vector<T>::template spawn<complexity::N>(a, b);
    }

    template<typename T, class A> 
    inline double max_vector(const matrix<T,A>& a){ 
        ambient::numeric::future<double> ret(0.0);
        kernels::template max_vector<T>::template spawn<complexity::N>(a, ret);
        return (double)ret; 
    }

    template<class Matrix>
    inline void fill_gaussian(Matrix& a){
        kernels::template init_gaussian<typename Matrix::value_type>::template spawn<complexity::N2>(a);
    }

} }

#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type
#endif
