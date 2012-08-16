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
    bool is_hermitian(matrix<T> const& m)
    {
        /*if(num_rows(m) != num_cols(m))
            return false;
        for (size_t i=0; i<num_rows(m); ++i)
            for(size_t j=0; j<num_cols(m); ++j)
                if ( m(i,j) != conj(m(j,i)) )
                    return false;*/
        return false;
    }


    template<class MatrixViewA, class MatrixViewB, typename T>
    inline void gemm(const MatrixViewA& a, const MatrixViewB& b, matrix<T>& c){
        kernels::gemm_general<MatrixViewA,MatrixViewB,T>::spawn(a, b, c); 
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
    inline void svd(matrix<T> a, matrix<T>& u, matrix<T>& vt, diagonal_matrix<double>& s){
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
    inline matrix<T> exp(const matrix<T>& m, const T& alfa = 1.){
        printf("ERROR: NOT TESTED (EXP)\n");
        diagonal_matrix<double> evals(m.num_rows());
        matrix<T> evecs = matrix<T>();
        heev(m, evecs, evals);
        return (evecs * exp(evals, alfa))*conj(transpose(evecs));
    }

    template<typename T> inline void resize(matrix<T>& m, size_t rows, size_t cols){ 
        assert(m.num_rows() != 0 && m.num_cols() != 0);
        if(m.num_rows() == rows && m.num_cols() == cols) return;
        matrix<T> resized(rows, cols);
        if(!m.impl->weak())
            kernels::resize<T>::spawn(resized, m, std::min(rows, m.num_rows()), std::min(cols, m.num_cols()));
        m.swap(resized);
    }

    template<typename T> inline scalar_type trace(const matrix<T>& m){ 
        scalar_type trace(0.);
        kernels::trace<T>::spawn(m, trace);
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
    inline const matrix<T>& conj(const matrix<T>& m){
        //m.conj();
        return m;
    }

    template<typename T>
    inline void conj_inplace(matrix<T>& m){
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
    inline matrix<T> transpose(const matrix<T>& a){
        matrix<T> t(a.num_cols(), a.num_rows());
        kernels::transpose_out<T>::spawn(a, t);
        return t;
    }

    template<typename T>
    inline void adjoint_inplace(matrix<T>& m){
        transpose_inplace(m);
    }

    template<typename T>
    inline const matrix<T>& adjoint(const matrix<T>& m){
        return transpose(m);
    }

    template<typename T, class G>
    inline void generate(matrix<T>& m, G g){ 
        fill_random(m);  // G isn't used
    }

    template<typename T>
    inline void fill_identity(matrix<T>& m){
        kernels::init_identity<T>::spawn(m);
    }

    template<typename T>
    inline void fill_random(matrix<T>& m){
        kernels::init_random<T>::spawn(m);
    }

    template<typename T>
    inline void fill_value(matrix<T>& m, T value){
        if(value == 0.) return; // matrices are 0 by default
        kernels::init_value<T>::spawn(m, value);
    }

    template<typename T>
    inline void remove_rows(matrix<T>& m, size_type i, difference_type k){
        assert( i+k <= m->num_rows() );
        assert(false); printf("ERROR: NOT TESTED (REMOVE ROWS)\n");
        //kernels::remove_rows<T>::spawn(m, i, k);
        resize(m, m.num_rows()-k, m.num_cols());
    }

    template<typename T>
    inline void remove_cols(matrix<T>& m, size_type j, difference_type k){
        assert( j+k <= m->num_cols() );
        assert(false); printf("ERROR: NOT TESTED (REMOVE COLS)\n");
        //kernels::remove_cols<T>::spawn(m, j, k);
        resize(m, m.num_rows(), m.num_cols()-k);
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
    inline void mul_inplace(matrix<T>& m, const matrix<T>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM INPLACE)\n");
        //kernels::gemm_inplace<T>::spawn(m, rhs);
    }

    template <typename T>
    inline void mul_inplace(matrix<T>& m, const diagonal_matrix<T>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM DIAG INPLACE)\n");
        //kernels::gemm_diag_inplace<T>::spawn(m, rhs);
    }

    template <typename T>
    inline void mul_inplace(matrix<T>& a, const scalar_type& rhs) { 
        kernels::scale<T>::spawn(a, rhs); 
    }

    template <typename T>
    inline void div_inplace(matrix<T>& a, const scalar_type& rhs) { 
        kernels::scale_inverse<T>::spawn(a, rhs); 
    }

    template<typename T>
    inline void copy(matrix<T>& ac, const matrix<T>& a){ 
        kernels::copy<T>::spawn(ac, a); 
    }

    template<typename T>
    inline void destroy(matrix<T>& a){
        ambient::destroy(a.impl); 
    }

    template<typename T>
    bool operator == (const matrix<T>& a, const matrix<T>& b){
        ambient::future<bool> ret(true);
        kernels::validation<T>::spawn(a, b, ret); 
        return (bool)ret;
    }

    template <typename T>
    std::ostream& operator << (std::ostream& o, matrix<T> const& m){
        for(size_type i=0; i< m.num_rows(); ++i){
            for(size_type j=0; j < m.num_cols(); ++j)
                ambient::cout << m(i,j) << " ";
            ambient::cout << std::endl;
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

    template<typename T> inline size_type num_rows(const matrix<T>& m){ 
        return m.num_rows(); 
    }

    template<typename T> inline size_type num_cols(const matrix<T>& m){
        return m.num_cols(); 
    }

    template<typename T> inline size_type num_rows(const transpose_view< matrix<T> >& m){ 
        return m.impl->spec.dim.x; 
    }

    template<typename T> inline size_type num_cols(const transpose_view< matrix<T> >& m){ 
        return m.impl->spec.dim.y; 
    }

} }

#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type 
#endif
