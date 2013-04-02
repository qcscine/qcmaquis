#ifndef AMBIENT_NUMERIC_MATRIX_ALGORITHMS_EXPERIMENTAL
#define AMBIENT_NUMERIC_MATRIX_ALGORITHMS_EXPERIMENTAL

#include "ambient/numeric/matrix/matrix.h"
#include "ambient/numeric/kernels/kernels.hpp"
#include "ambient/numeric/matrix/diagonal_matrix.hpp"

#define size_type       typename matrix<T>::size_type
#define real_type       typename matrix<T>::real_type
#define scalar_type     typename matrix<T>::scalar_type
#define difference_type typename matrix<T>::difference_type

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

    template<typename T>
    inline void gbbrd(matrix<T>& a, diagonal_matrix<T>& d, diagonal_matrix<T>& e, matrix<T>& u, matrix<T>& v){
        kernels::gbbrd<T>::spawn<complexity::N3>(a, d, e, u, v);
    }

    template<typename T>
    inline void gebrd(matrix<T>& a, diagonal_matrix<T>& d, diagonal_matrix<T>& e, matrix<T>& u, matrix<T>& v){
        kernels::gebrd<T>::spawn<complexity::N3>(a, d, e, u, v);
    }

    template<typename T>
    inline void bdsqr(diagonal_matrix<T>& d, diagonal_matrix<T>& e, matrix<T>& u, matrix<T>& v){
        kernels::bdsqr<T>::spawn<complexity::N3>(d, e, u, v);
    }

    template<typename T>
    inline void gebd2(matrix<T>& a, diagonal_matrix<T>& d, diagonal_matrix<T>& e, diagonal_matrix<T>& tq, diagonal_matrix<T>& tp){
        kernels::gebd2<T>::spawn<complexity::N3>(a, d, e, tq, tp);
    }

    template<PLASMA_enum TR, typename T>
    inline void larfg(matrix<T>& a, diagonal_matrix<T>& t, diagonal_matrix<T>& d, size_t k){
        kernels::larfg<T,TR>::spawn<complexity::N3>(a, t, d, k);
    }

    template<typename T>
    inline void labrd_update_col(matrix<T>& say, const matrix<T>& sax, 
                           matrix<T>& sy, const matrix<T>& sx, 
                           diagonal_matrix<T>& tq, 
                           diagonal_matrix<T>& d, 
                           int i)
    {
        kernels::labrd_update_col<T>::spawn<complexity::N3>(say, sax, sy, sx, tq, d, i);
    }

    template<typename T>
    inline void labrd_reduce_col(matrix<T>& say, const matrix<T>& sax, 
                           matrix<T>& sy, const matrix<T>& sx, 
                           int i)
    {
        kernels::labrd_reduce_col<T>::spawn<complexity::N3>(say, sax, sy, sx, i);
    }

    template<typename T>
    inline void labrd_update_row(const matrix<T>& say, matrix<T>& sax, 
                           const matrix<T>& sy, matrix<T>& sx, 
                           diagonal_matrix<T>& tp, 
                           diagonal_matrix<T>& e, 
                           int i)
    {
        kernels::labrd_update_row<T>::spawn<complexity::N3>(say, sax, sy, sx, tp, e, i);
    }

    template<typename T>
    inline void labrd_reduce_row(const matrix<T>& say, matrix<T>& sax, 
                           const matrix<T>& sy, matrix<T>& sx, 
                           int i)
    {
        kernels::labrd_reduce_row<T>::spawn<complexity::N3>(say, sax, sy, sx, i);
    }

    template<PLASMA_enum UL, size_t OFF, typename T>
    inline void laset2(matrix<T>& a, const T& alfa = 0.0){
        kernels::laset2<T,UL,OFF>::spawn<complexity::N2>(a, alfa);
    }

    template<PLASMA_enum UL, typename T>
    inline void copy_band(const matrix<T>& src, matrix<T>& dst, size_t dj){
        kernels::copy_band<T,UL>::spawn<complexity::N2>(src, dst, dj);
    }

    template <int alfa, typename T>
    inline void add_vectors(matrix<T>& lhs, size_t loffset, const matrix<T>& rhs, size_t roffset, size_t size){ 
        kernels::add_vectors<alfa, T>::spawn<complexity::N2>(lhs, loffset, rhs, roffset, size); 
    }

} }

#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type
#endif
