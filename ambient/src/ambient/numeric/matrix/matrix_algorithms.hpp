#ifndef __AMBIENT_NUMERIC_MATRIX_ALGORITHMS_HPP__
#define __AMBIENT_NUMERIC_MATRIX_ALGORITHMS_HPP__

#include "ambient/numeric/matrix/matrix.h"
#include "ambient/numeric/matrix/kernels/utils.hpp"
#include "ambient/numeric/matrix/kernels/kernels.hpp"
#include "ambient/numeric/matrix/diagonal_matrix.hpp"

#define size_type       typename matrix<T>::size_type
#define real_type       typename matrix<T>::real_type
#define scalar_type     typename matrix<T>::scalar_type
#define difference_type typename matrix<T>::difference_type

// Functions used by the ground state algorithm: //
//
// norm_square
// overlap
// transpose
// transpose_inplace
// trace
// add_inplace
// sub_inplace
// resize
// mul_inplace (scalar)
// div_inplace
// copy
// gemm (matrix-matrix)
// gemm (matrix-diagonal)
// svd
// heev

namespace ambient { namespace numeric {

    template<typename T>
    bool is_hermitian(const matrix<T>& a)
    {
        /*if(num_rows(a) != num_cols(a))
            return false;
        for (size_t i=0; i < num_rows(a); ++i)
            for(size_t j=0; j < num_cols(a); ++j)
                if ( a(i,j) != conj(a(j,i)) )
                    return false;*/
        return false;
    }


    template<int alfa, int beta, class MatrixViewA, class MatrixViewB, class MatrixViewC>
    inline void gemv(const MatrixViewA& a, size_t aoffset, 
                     const MatrixViewB& b, size_t boffset, 
                           MatrixViewC& c, size_t coffset, 
                           size_t m, size_t n)
    {
        kernels::gemv<alfa,beta,MatrixViewA,MatrixViewB,MatrixViewC>::spawn(a, aoffset, 
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
        kernels::gemv_scale<ADD,MA,MB,MC,MF>::spawn(a, aoffset, 
                                                    b, boffset, 
                                                    c, coffset, 
                                                    f, foffset, 
                                                    m, n); 
    }

    template<class MatrixViewA, class MatrixViewB, typename T>
    inline void gemm(const MatrixViewA& a, const MatrixViewB& b, matrix<T>& c){
        kernels::gemm<MatrixViewA,MatrixViewB,T>::spawn(a, b, c); 
    }

    template<class MatrixViewA, typename T, typename D>
    inline void gemm(const MatrixViewA& a, const diagonal_matrix<D>& b, matrix<T>& c){ 
        kernels::gemm_diagonal_rhs<MatrixViewA,T,D>::spawn(a, b, c); 
    }

    template<class MatrixViewB, typename T, typename D>
    inline void gemm(const diagonal_matrix<D>& a, const MatrixViewB& b, matrix<T>& c){ 
        kernels::gemm_diagonal_lhs<MatrixViewB,T,D>::spawn(a, b, c); 
    }

    template<typename T>
    inline void scale(matrix<T>& a, size_t ai, size_t aj, const diagonal_matrix<T>& alfa, size_t alfai, size_t alfaj){
        kernels::scale_offset<T>::spawn(a, ai, aj, alfa, alfai);
    }

    template<typename T>
    inline void gbbrd(matrix<T>& a, diagonal_matrix<T>& d, diagonal_matrix<T>& e, matrix<T>& u, matrix<T>& v){
        kernels::gbbrd<T>::spawn(a, d, e, u, v);
    }

    template<typename T>
    inline void gebrd(matrix<T>& a, diagonal_matrix<T>& d, diagonal_matrix<T>& e, matrix<T>& u, matrix<T>& v){
        kernels::gebrd<T>::spawn(a, d, e, u, v);
    }

    template<typename T>
    inline void bdsqr(diagonal_matrix<T>& d, diagonal_matrix<T>& e, matrix<T>& u, matrix<T>& v){
        kernels::bdsqr<T>::spawn(d, e, u, v);
    }

    template<typename T>
    inline void gebd2(matrix<T>& a, diagonal_matrix<T>& d, diagonal_matrix<T>& e, diagonal_matrix<T>& tq, diagonal_matrix<T>& tp){
        kernels::gebd2<T>::spawn(a, d, e, tq, tp);
    }

    template<PLASMA_enum TR, typename T>
    inline void larfg(matrix<T>& a, diagonal_matrix<T>& t, diagonal_matrix<T>& d, size_t k){
        kernels::larfg<T,TR>::spawn(a, t, d, k);
    }

    template<typename T>
    inline void labrd_update_col(matrix<T>& say, const matrix<T>& sax, 
                           matrix<T>& sy, const matrix<T>& sx, 
                           diagonal_matrix<T>& tq, 
                           diagonal_matrix<T>& d, 
                           int i)
    {
        kernels::labrd_update_col<T>::spawn(say, sax, sy, sx, tq, d, i);
    }

    template<typename T>
    inline void labrd_reduce_col(matrix<T>& say, const matrix<T>& sax, 
                           matrix<T>& sy, const matrix<T>& sx, 
                           int i)
    {
        kernels::labrd_reduce_col<T>::spawn(say, sax, sy, sx, i);
    }

    template<typename T>
    inline void labrd_update_row(const matrix<T>& say, matrix<T>& sax, 
                           const matrix<T>& sy, matrix<T>& sx, 
                           diagonal_matrix<T>& tp, 
                           diagonal_matrix<T>& e, 
                           int i)
    {
        kernels::labrd_update_row<T>::spawn(say, sax, sy, sx, tp, e, i);
    }

    template<typename T>
    inline void labrd_reduce_row(const matrix<T>& say, matrix<T>& sax, 
                           const matrix<T>& sy, matrix<T>& sx, 
                           int i)
    {
        kernels::labrd_reduce_row<T>::spawn(say, sax, sy, sx, i);
    }

    template<typename T>
    inline void svd(matrix<T>& a, matrix<T>& u, matrix<T>& vt, diagonal_matrix<double>& s){
        int m = num_rows(a);
        int n = num_cols(a);
        int k = std::min(m,n);
        u.resize(m, k);
        vt.resize(k, n);
        s.resize(k, k);
        kernels::svd<T>::spawn(a, u, vt, s);
    }

    template<typename T>
    inline void heev(matrix<T> a, matrix<T>& evecs, diagonal_matrix<double>& evals){
        assert(num_rows(a) == num_cols(a) && num_rows(evals) == num_rows(a));
        kernels::heev<T>::spawn(a, evals); // destoys U triangle of M
        evecs.swap(a);
    }

    template<typename T>
    inline void syev(const matrix<T>& a, matrix<T>& evecs, diagonal_matrix<double>& evals){
        heev(a, evecs, evals);
    }

    template<typename T>
    inline void geqrt(matrix<T>& a, matrix<T>& t){
        kernels::geqrt<T>::spawn(a, t);
    }

    template<PLASMA_enum TR, typename T>
    inline void ormqr(size_t k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c){
        kernels::ormqr<T,TR>::spawn(k, a, t, c);
    }

    template<typename T>
    inline void tsqrt(matrix<T>& a1, matrix<T>& a2, matrix<T>& t){
        kernels::tsqrt<T>::spawn(a1, a2, t);
    }

    template<PLASMA_enum TR, typename T>
    inline void tsmqr(size_t k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t){
        kernels::tsmqr<T,TR>::spawn(k, a1, a2, v, t);
    }

    template<typename T>
    inline void gelqt(matrix<T>& a, matrix<T>& t){
        kernels::gelqt<T>::spawn(a, t);
    }

    template<PLASMA_enum TR, typename T>
    inline void ormlq(size_t k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c){
        kernels::ormlq<T,TR>::spawn(k, a, t, c);
    }

    template<typename T>
    inline void tslqt(matrix<T>& a1, matrix<T>& a2, matrix<T>& t){
        kernels::tslqt<T>::spawn(a1, a2, t);
    }

    template<PLASMA_enum TR, typename T>
    inline void tsmlq(size_t k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t){
        kernels::tsmlq<T,TR>::spawn(k, a1, a2, v, t);
    }

    template<PLASMA_enum UL, size_t OFF, typename T>
    inline void laset2(matrix<T>& a, const T& alfa = 0.0){
        kernels::laset2<T,UL,OFF>::spawn(a, alfa);
    }

    template<typename T> 
    inline void qr(matrix<T> a, matrix<T>& q, matrix<T>& r){ 
        int m = num_rows(a);            
        int n = num_cols(a);            
        int k = std::min(m,n); 
        resize(q, m, k);  
        resize(r, k, n); 
        kernels::qr<T>::spawn(a, q, r);  
    }

    template<typename T> 
    inline void lq(matrix<T> a, matrix<T>& l, matrix<T>& q){ 
        int m = num_rows(a);            
        int n = num_cols(a);            
        int k = std::min(m,n); 
        resize(l, m, k); 
        resize(q, k, n);  
        kernels::lq<T>::spawn(a, l, q);  
    } 

    template<typename T>
    inline matrix<T> exp(const matrix<T>& a, const T& alfa = 1.){
        printf("ERROR: NOT TESTED (EXP)\n");
        diagonal_matrix<double> evals(a.num_rows());
        matrix<T> evecs = matrix<T>();
        heev(a, evecs, evals);
        return (evecs * exp(evals, alfa))*conj(transpose(evecs));
    }

    template<typename T> inline void resize(matrix<T>& a, size_t m, size_t n){ 
        assert(a.num_rows() != 0 && a.num_cols() != 0);
        if(a.num_rows() == m && a.num_cols() == n) return;
        matrix<T> resized(m, n);
        if(!a.impl->weak())
            kernels::resize<T>::spawn(resized, a, std::min(m, a.num_rows()), std::min(n, a.num_cols()));
        a.swap(resized);
    }

    template<typename T> inline scalar_type trace(const matrix<T>& a){ 
        scalar_type trace(0.);
        kernels::trace<T>::spawn(a, trace);
        return trace;
    }

    template <typename T>
    inline real_type norm_square(const matrix<T>& a){ 
        real_type norm(0.); 
        kernels::scalar_norm<T>::spawn(a, norm); 
        return norm; 
    }

    template <typename T>
    inline scalar_type overlap(const matrix<T>& a, const matrix<T>& b){ 
        scalar_type overlap(0.); 
        kernels::overlap<T>::spawn(a, b, overlap); 
        return overlap; 
    }
        
    template<typename T>
    inline void swap(matrix<T>& x, matrix<T>& y){ 
        x.swap(y);                     
    }

    template<typename T>
    inline const matrix<T>& conj(const matrix<T>& a){
        //a.conj();
        return a;
    }

    template<typename T>
    inline void conj_inplace(matrix<T>& a){
        // gs (doubles)
        // does nothing for now
    }

    template<typename T>
    inline void transpose_inplace(matrix<T>& a){
        matrix<T> t(a.num_cols(), a.num_rows());
        kernels::transpose_out<T>::spawn(a, t);
        a.swap(t);
    }

    template<typename T>
    inline transpose_view<matrix<T> > transpose(const matrix<T>& a){
        return transpose_view<matrix<T> >(a);
    }

    template<typename T>
    inline void adjoint_inplace(matrix<T>& a){
        transpose_inplace(a);
    }

    template<typename T>
    inline transpose_view<matrix<T> > adjoint(const matrix<T>& a){
        return transpose_view<matrix<T> >(a);
    }

    template<class Matrix>
    inline void fill_identity(Matrix& a){
        kernels::init_identity<typename Matrix::value_type>::spawn(a);
    }

    template<class Matrix>
    inline void fill_random(Matrix& a){
        kernels::init_random<typename Matrix::value_type>::spawn(a);
    }

    template<class Matrix>
    inline void fill_value(Matrix& a, typename Matrix::value_type value){
        if(value == 0.) return; // matrices are 0 by default
        kernels::init_value<typename Matrix::value_type>::spawn(a, value);
    }

    template<typename T>
    inline void remove_rows(matrix<T>& a, size_type i, difference_type k){
        assert( i+k <= a.num_rows() );
        assert(false); printf("ERROR: NOT TESTED (REMOVE ROWS)\n");
        //kernels::remove_rows<T>::spawn(a, i, k);
        resize(a, a.num_rows()-k, a.num_cols());
    }

    template<typename T>
    inline void remove_cols(matrix<T>& a, size_type j, difference_type k){
        assert( j+k <= a.num_cols() );
        assert(false); printf("ERROR: NOT TESTED (REMOVE COLS)\n");
        //kernels::remove_cols<T>::spawn(a, j, k);
        resize(a, a.num_rows(), a.num_cols()-k);
    }

    template <int alfa, typename T>
    inline void add_vectors(matrix<T>& lhs, size_t loffset, const matrix<T>& rhs, size_t roffset, size_t size){ 
        kernels::add_vectors<alfa, T>::spawn(lhs, loffset, rhs, roffset, size); 
    }

    template <typename T>
    inline void add_inplace(matrix<T>& lhs, const matrix<T>& rhs){ 
        kernels::add<T>::spawn(lhs, rhs); 
    }

    template <typename T>
    inline void sub_inplace(matrix<T>& lhs, const matrix<T>& rhs){ 
        kernels::sub<T>::spawn(lhs, rhs); 
    }

    template <typename T>
    inline void mul_inplace(matrix<T>& a, const matrix<T>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM INPLACE)\n");
        //kernels::gemm_inplace<T>::spawn(a, rhs);
    }

    template <typename T>
    inline void mul_inplace(matrix<T>& a, const diagonal_matrix<T>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM DIAG INPLACE)\n");
        //kernels::gemm_diag_inplace<T>::spawn(a, rhs);
    }

    template <typename T>
    inline void mul_inplace(matrix<T>& a, const scalar_type& rhs) { 
        kernels::scale<T>::spawn(a, rhs); 
    }

    template <typename T>
    inline void div_inplace(matrix<T>& a, const scalar_type& rhs) { 
        kernels::scale_inverse<T>::spawn(a, rhs); 
    }

    template <typename T>
    inline void save(const matrix<T>& a, const size_t& tag) { 
        kernels::save<T>::spawn(a, tag); 
    }

    template <typename T>
    inline void load(matrix<T>& a, const size_t& tag) { 
        kernels::load<T>::spawn(a, tag); 
    }

    template<typename T>
    inline void copy(const matrix<T>& src, matrix<T>& dst){
        ambient::fuse(src.impl, dst.impl);
    }

    template<typename T>
    inline void copy_block(const matrix<T>& src, size_t si, size_t sj, 
                           matrix<T>& dst, size_t di, size_t dj, 
                           size_t m, size_t n)
    {
        kernels::copy_block<T>::spawn(src, si, sj, dst, di, dj, m, n); 
    }

    template<PLASMA_enum UL, typename T>
    inline void copy_band(const matrix<T>& src, matrix<T>& dst, size_t dj){
        kernels::copy_band<T,UL>::spawn(src, dst, dj);
    }

    template<typename T>
    inline void copy_rt(const matrix<T>& src, matrix<T>& dst){ 
        kernels::copy_rt<T>::spawn(src, dst);
    }

    template<typename T>
    inline void copy_lt(const matrix<T>& src, matrix<T>& dst){ 
        kernels::copy_lt<T>::spawn(src, dst);
    }

    template<typename T>
    inline void copy_block_s(const matrix<T>& src, size_t si, size_t sj, 
                             matrix<T>& dst, size_t di, size_t dj, 
                             const matrix<T>& alfa, size_t ai, size_t aj,
                             size_t m, size_t n)
    { 
        kernels::copy_block_s<T>::spawn(src, si, sj, dst, di, dj, alfa, ai, aj, m, n);
    }

    template<typename T>
    inline void copy_block_sa(const matrix<T>& src, size_t si, size_t sj, 
                              matrix<T>& dst, size_t di, size_t dj, 
                              const matrix<T>& alfa, size_t ai, size_t aj,
                              size_t m, size_t n)
    { 
        kernels::copy_block_sa<T>::spawn(src, si, sj, dst, di, dj, alfa, ai, aj, m, n);
    }

    template<typename T>
    inline void destroy(matrix<T>& a){
        ambient::destroy(a.impl); 
    }

    template<typename T>
    bool operator == (const matrix<T>& a, const matrix<T>& b){
        if(num_cols(a) != num_cols(b) || num_rows(a) != num_rows(b)){
            printf("Sizes are different: %lu x %lu against %lu x %lu\n",
                    num_cols(a), num_rows(a), num_cols(b), num_rows(b));
            return false;
        }
        ambient::future<bool> ret(true);
        kernels::validation<T>::spawn(a, b, ret);
        return (bool)ret;
    }

    template<typename T>
    bool operator == (matrix<T> a, const transpose_view<matrix<T> >& b){
        if(num_cols(a) != num_cols(b) || num_rows(a) != num_rows(b)){
            printf("Sizes are different: %lu x %lu against %lu x %lu\n",
                    num_cols(a), num_rows(a), num_cols(b), num_rows(b));
            return false;
        }
        transpose_inplace(a);
        ambient::future<bool> ret(true);
        kernels::validation<T>::spawn(a, b, ret);
        return (bool)ret;
    }

    template<typename T>
    bool operator == (const diagonal_matrix<T>& a, const diagonal_matrix<T>& b){
        if(num_rows(a) != num_rows(b)){
            printf("Sizes are different: %lu x %lu against %lu x %lu\n",
                    num_cols(a), num_rows(a), num_cols(b), num_rows(b));
            return false;
        }
        ambient::future<bool> ret(true);
        kernels::validation<T>::spawn(a, b, ret);
        return (bool)ret;
    }

    template <typename T>
    std::ostream& operator << (std::ostream& o, const matrix<T>& a){
        for(size_type i=0; i< a.num_rows(); ++i){
            for(size_type j=0; j < a.num_cols(); ++j){
                //if(a(i,j) < 10e-16) ambient::cout << 0 << " ";
                //else 
                if(a(i,j)*a(i,j) < 10e-40) printf("     ");
                else if(a(i,j) < 0) printf("%.2f ", -1*a(i,j)); //ambient::cout << a(i,j) << " ";
                else printf("%.2f ", a(i,j)); //ambient::cout << a(i,j) << " ";
            }
            ambient::cout << std::endl << "\n";
        }
        return o;
    }

    template <typename T> 
    inline matrix<T> operator + (matrix<T> lhs, const matrix<T>& rhs){
        return (lhs += rhs);
    }

    template <typename T> 
    inline matrix<T> operator - (matrix<T> lhs, const matrix<T>& rhs){
        return (lhs -= rhs);
    }

    template <typename T>
    inline const matrix<T> operator * (matrix<T> lhs, const matrix<T>& rhs){
        return (lhs *= rhs);
    }

    template<typename T, typename T2> 
    inline const matrix<T> operator * (matrix<T> lhs, const T2& rhs){
        return (lhs *= rhs);
    }

    template<typename T, typename T2> 
    inline const matrix<T> operator * (const T2& lhs, matrix<T> rhs){
        return (rhs *= lhs);
    }

    template<typename T> inline size_type num_rows(const matrix<T>& a){
        return a.num_rows();
    }

    template<typename T> inline size_type num_cols(const matrix<T>& a){
        return a.num_cols();
    }

    template<typename T> inline size_type num_rows(const transpose_view< matrix<T> >& a){
        return a.impl->spec.dim.x;
    }

    template<typename T> inline size_type num_cols(const transpose_view< matrix<T> >& a){
        return a.impl->spec.dim.y;
    }

} }

#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type
#endif
