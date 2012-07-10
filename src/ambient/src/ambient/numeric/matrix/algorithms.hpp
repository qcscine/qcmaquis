#ifndef __AMBIENT_NUMERIC_MATRIX_ALGORITHMS_HPP__
#define __AMBIENT_NUMERIC_MATRIX_ALGORITHMS_HPP__

#include "ambient/numeric/matrix/kernels/utils.hpp"
#include "ambient/numeric/matrix/kernels/kernels.hpp"
#include "ambient/numeric/matrix/kernels/atomics.hpp"

#define ATOMIC(condition, kernel, ...) assert(condition); ambient::push< kernels::kernel ## _atomic<T> >(__VA_ARGS__);

#define size_type       typename matrix<T>::size_type
#define real_type       typename matrix<T>::real_type
#define scalar_type     typename matrix<T>::scalar_type
#define difference_type typename matrix<T>::difference_type

namespace ambient { namespace numeric {

    // {{{ matrix advanced algorithms

    template <typename T>
    inline real_type norm_square(const matrix<T>& a){ real_type norm(0); ATOMIC(a.atomic(), scalar_norm, a, norm); return norm; } // gs

    template <typename T>
    inline scalar_type overlap(const matrix<T>& a, const matrix<T>& b){ scalar_type overlap(0); ATOMIC(a.atomic(), overlap, a, b, overlap); return overlap; } // gs

    // }}}

    // {{{ matrix operators (free functions)

    template <typename T>
    std::ostream& operator << (std::ostream& o, matrix<T> const& m){
        for(size_type i=0; i< m.num_rows(); ++i){
            for(size_type j=0; j < m.num_cols(); ++j)
                ambient::cout << m(i,j) << " ";
            ambient::cout << std::endl;
        }
        return o;
    }

    template <typename T> inline matrix<T> operator + (matrix<T> lhs, const matrix<T>& rhs){ return (lhs += rhs); }
    template <typename T> inline matrix<T> operator - (matrix<T> lhs, const matrix<T>& rhs){ return (lhs -= rhs); }
    template <typename T> inline const matrix<T> operator * (matrix<T> lhs, const matrix<T>& rhs){ return (lhs *= rhs); }
    template<typename T, typename T2> inline const matrix<T> operator * (matrix<T> lhs, const T2& rhs){ return (lhs *= rhs); }
    template<typename T, typename T2> inline const matrix<T> operator * (const T2& lhs, matrix<T> rhs){ return (rhs *= lhs); }
    template<typename T> inline size_type num_rows(const matrix<T>& m){ return m.num_rows(); }
    template<typename T> inline size_type num_cols(const matrix<T>& m){ return m.num_cols(); }
    template<typename T> inline void resize(matrix<T>& m, size_t rows, size_t cols){ m.resize(rows, cols); }
    template<typename T> inline scalar_type trace(const matrix<T>& m){ return m.trace(); }

    template<typename T>
    inline void transpose_inplace(matrix<T>& a){
        // gs
        matrix<T> t(a.num_cols(), a.num_rows());
        ATOMIC(a.atomic(), transpose_out, a, t);
        a.swap(t);
    }

    template<typename T>
    inline matrix<T> transpose(const matrix<T>& a){
        // gs
        matrix<T> t(a.num_cols(), a.num_rows());
        ATOMIC(a.atomic(), transpose_out, a, t);
        return t;
    }

    template<typename T>
    inline void conj_inplace(matrix<T>& m){
        //m.conj();
    }

    template<typename T>
    inline const matrix<T>& conj(const matrix<T>& m){
        //m.conj();
        return m;
    }

    template<typename T>
    inline void adjoint_inplace(matrix<T>& m){
        transpose_inplace(m);
        //m.conj();
    }

    template<typename T>
    inline const matrix<T>& adjoint(const matrix<T>& m){
        return transpose(m);
        //m.conj();
    }

    template<typename T>
    inline matrix<T> exp(const matrix<T>& m, const T& alfa = 1.){
        printf("Error: Not checked <- EXP\n");
        typename alps::numeric::associated_real_diagonal_matrix< matrix<T> >::type evals(m.num_rows());
        matrix<T> evecs = matrix<T>();
        heev(m, evecs, evals);
        return (evecs * exp(evals, alfa))*conj(transpose(evecs));
    }

    // {{{ strassen matrix multiplication algorithm

    // {{{ strassen multiplication supplementary functions
    /*
    template<typename T>
    void gemm_strassen_gad(const matrix<T>& a, size_t ai, size_t aj, 
                           matrix<T>& r, size_t n)
    {
        ambient::push(ambient::gemm_strassen_gad_l<T>,  // output reduced matrix
                      ambient::gemm_strassen_gad_c<T>,  // a11 + a12,  a12 - a22
                      a, ai, aj, r, n);                 // a21 - a11,  a22 + a21
    }

    template<typename T>
    void gemm_strassen_dad(const matrix<T>& a, size_t ai, size_t aj,
                           const matrix<T>& b, size_t bi, size_t bj, 
                           matrix<T>& r, size_t n)
    {
        ambient::push(ambient::gemm_strassen_dad_l<T>,  // output reduced matrix
                      ambient::gemm_strassen_dad_c<T>,  // a11 + a22,  b11 + b22
                      a, ai, aj, b, bi, bj, r, n);

    }

    template<typename T>
    void gemm_strassen_pluseqs(const matrix<T>& a, size_t ai, size_t aj,
                               const matrix<T>& b, size_t bi, size_t bj, 
                                     matrix<T>& c, size_t ci, size_t cj, 
                                     size_t n)
    {
        ambient::push(ambient::add_sum_submx_l<T>, 
                      ambient::add_sum_submx_c<T>, 
                      a, ai, aj, b, bi, bj, c, ci, cj, n);
    }

    template<typename T>
    void gemm_strassen_pluseqd(const matrix<T>& a, size_t ai, size_t aj,
                               const matrix<T>& b, size_t bi, size_t bj, 
                                     matrix<T>& c, size_t ci, size_t cj, 
                                     size_t n)
    {
        ambient::push(ambient::add_dif_submx_l<T>, 
                      ambient::add_dif_submx_c<T>, 
                      a, ai, aj, b, bi, bj, c, ci, cj, n);
    }
    // }}}

    template<typename T>
    void gemm_strassen(const matrix<T>& a, size_t ai, size_t aj,
                       const matrix<T>& b, size_t bi, size_t bj,
                             matrix<T>& c, size_t ci, size_t cj, 
                             size_t n)
    {

        if(n > 128){ 
            // strassen algorithm recursion
            matrix<T> m1(n/2, n/2);
            matrix<T> m2(n/2, n/2);
            matrix<T> m3(n/2, n/2);
            matrix<T> m4(n/2, n/2);
            matrix<T> m5(n/2, n/2);
            matrix<T> m6(n/2, n/2);
            matrix<T> m7(n/2, n/2);

            matrix<T> ar(n,n);
            matrix<T> br(n,n);
            matrix<T> dr(n/2,n);

            gemm_strassen_gad(a, ai, aj, ar, n);
            gemm_strassen_gad(b, bi, bj, br, n);
            gemm_strassen_dad(a, ai, aj, b, bi, bj, dr, n);

            gemm_strassen( dr , 0      , 0      , dr , 0      , n/2    , m1, 0, 0, n/2 );
            gemm_strassen( ar , n/2    , n/2    , b  , bi     , bj     , m2, 0, 0, n/2 );
            gemm_strassen( a  , ai     , aj     , br , 0      , n/2    , m3, 0, 0, n/2 );
            gemm_strassen( a  , ai+n/2 , aj+n/2 , br , n/2    , 0      , m4, 0, 0, n/2 );
            gemm_strassen( ar , 0      , 0      , b  , bi+n/2 , bj+n/2 , m5, 0, 0, n/2 );
            gemm_strassen( ar , n/2    , 0      , br , 0      , 0      , m6, 0, 0, n/2 );
            gemm_strassen( ar , 0      , n/2    , br , n/2    , n/2    , m7, 0, 0, n/2 );

            gemm_strassen_pluseqs( m2, 0, 0, m4, 0, 0, c , ci+n/2 , cj     , n/2 );
            gemm_strassen_pluseqs( m3, 0, 0, m5, 0, 0, c , ci     , cj+n/2 , n/2 );
            gemm_strassen_pluseqd( m1, 0, 0, m5, 0, 0, m4, 0      , 0      , n/2 );
            gemm_strassen_pluseqd( m1, 0, 0, m2, 0, 0, m3, 0      , 0      , n/2 );
            gemm_strassen_pluseqs( m4, 0, 0, m7, 0, 0, c , ci     , cj     , n/2 );
            gemm_strassen_pluseqs( m3, 0, 0, m6, 0, 0, c , ci+n/2 , cj+n/2 , n/2 );
        }else{
            // standard gemm (submatrices)
            ambient::push(ambient::gemm_submx_l<T>, ambient::gemm_submx_c<T>,
                          a, ai, aj, b, bi, bj, c, ci, cj, n);
        }
    }

    template<typename T>
    void gemm_strassen(const matrix<T>& a, const matrix<T>& b, matrix<T>& c){
        size_t n  = num_cols(a);
        assert(n == num_rows(a));
        assert(n == num_cols(b));
        assert(n == num_rows(b));
        c.resize(n, n);
        gemm_strassen(a, 0, 0, b, 0, 0, c, 0, 0, n);
    }*/
    // }}}

    template<class Tag1, class Tag2, typename T>
    inline void gemm(const matrix<T>& a, const matrix<T>& b, matrix<T>& c){ assert(num_cols(a) == num_rows(b)); ambient::push< kernels::gemm_general_atomic<Tag1,Tag2,T> >(a, b, c); } // gs

    template<class Tag1, class Tag2, typename T, typename D>
    inline void gemm(const matrix<T>& a, const diagonal_matrix<D>& b, matrix<T>& c){ assert(num_cols(a) == num_rows(b)); ambient::push< kernels::gemm_diagonal_rhs<Tag1,Tag2,T,D> >(a, b, c); } // gs

    template<class Tag1, class Tag2, typename T, typename D>
    inline void gemm(const diagonal_matrix<D>& a, const matrix<T>& b, matrix<T>& c){ assert(num_cols(a) == num_rows(b)); ambient::push< kernels::gemm_diagonal_lhs<Tag1,Tag2,T,D> >(a, b, c); } // gs

    template<typename T>
    inline void svd(matrix<T> a, matrix<T>& u, matrix<T>& vt, typename alps::numeric::associated_real_diagonal_matrix<matrix<T> >::type& s){ // gs
        int m = num_rows(a);
        int n = num_cols(a);
        int k = std::min(m,n);
        u.resize(m, k);
        vt.resize(k, n);
        s.resize(k, k);
        ATOMIC(a.atomic(), svd, a, u, vt, s);
    }

    template<typename T>
    inline void heev(matrix<T> a, matrix<T>& evecs, typename alps::numeric::associated_real_diagonal_matrix< matrix<T> >::type& evals){ // gs
        assert(num_rows(a) == num_cols(a) && num_rows(evals) == num_rows(a));
        ATOMIC(a.atomic(), heev, a, evals); // destoys U triangle of M
        evecs.swap(a);
    }

    template<typename T>
    inline void syev(const matrix<T>& a, matrix<T>& evecs, typename alps::numeric::associated_real_diagonal_matrix< matrix<T> >::type& evals){
        heev(a, evecs, evals);
    }

    template<typename T>
    inline void qr(const matrix<T>& m, matrix<T>& q, matrix<T>& r){
        assert(false); printf("Error: Not implemented <- QR\n");
    }

    template<typename T, class G>
    inline void generate(matrix<T>& m, G g){ m.fill_random(); } // warning: G isn't used
    // }}}

// {{{ implementation-specific type-nested algorithms //

namespace algorithms {

    template<typename T>
    inline void fill_identity(matrix_impl<T>& m){
        ATOMIC(m.atomic(), init_identity, m);
    }

    template<typename T>
    inline void fill_random(matrix_impl<T>& m){
        ATOMIC(m.atomic(), init_random, m);
    }

    template<typename T>
    inline void fill_value(matrix_impl<T>& m, T value){
        if(value == 0.) return; // matrices are 0s by default
        ATOMIC(m.atomic(), init_value, m, value);
    }

    template<typename T>
    inline void remove_rows(matrix_impl<T>& m, size_type i, difference_type k){
        assert(false); printf("Error: Not checked <- REMOVE ROWS\n");
        ambient::push< kernels::remove_rows<T> >(m, i, k);
    }

    template<typename T>
    inline void remove_cols(matrix_impl<T>& m, size_type j, difference_type k){
        assert(false); printf("Error: Not checked <- REMOVE COLS\n");
        ambient::push< kernels::remove_cols<T> >(m, j, k);
    }

    template<typename T>
    inline void resize(matrix_impl<T>& r, size_type rows, size_type cols, matrix_impl<T>& m, size_type orows, size_type ocols){ // gs
        ATOMIC(m.atomic(), resize, r, m, std::min(rows, orows), std::min(cols, ocols));
    }

    template<typename T>
    inline void conj_inplace(matrix_impl<T>& m){
        // gs (doubles)
        // does nothing for now
    }

    template<typename T>
    inline void transpose_inplace(matrix_impl<T>& m){
        assert(false); printf("Error: Not checked <- INPLACE TRANSPOSE\n");
        assert(m.num_rows() == m.num_cols()); // current limitation
        ambient::push< kernels::transpose<T> >(m);
    }

    template <typename T>
    inline scalar_type trace(const matrix_impl<T>& m){
        // gs
        scalar_type trace(0);
        ATOMIC(m.atomic(), trace, m, trace);
        return trace;
    }

    template <typename T>
    inline void add_inplace(matrix_impl<T>& m, const matrix_impl<T>& rhs){ ATOMIC(m.atomic(), add, m, rhs); } // gs

    template <typename T>
    inline void sub_inplace(matrix_impl<T>& m, const matrix_impl<T>& rhs){ ATOMIC(m.atomic(), sub, m, rhs); } // gs

    template <typename T>
    inline void scale_inplace(matrix_impl<T>& a, const scalar_type& rhs) { ATOMIC(a.atomic(), scale, a, rhs); } // gs

    template <typename T>
    inline void scale_inverse_inplace(matrix_impl<T>& a, const scalar_type& rhs) { ATOMIC(a.atomic(), scale_inverse, a, rhs); } // gs

    template<typename T>
    inline void copy(matrix_impl<T>& ac, const matrix_impl<T>& a){ ATOMIC(a.atomic(), copy, ac, a); } // gs

    template <typename T>
    inline void gemm_inplace(matrix_impl<T>& m, const matrix_impl<T>& rhs) {
        assert(false); printf("Error: Not checked <- GEMM INPLACE\n");
        ambient::push< kernels::gemm_inplace<T> >(m, rhs);
    }

    template <typename T>
    inline void gemm_diag_inplace(matrix_impl<T>& m, const matrix_impl<T>& rhs_diag) {
        assert(false); printf("Error: Not checked <- GEMM DIAG INPLACE\n");
        //ambient::push< kernels::gemm_diag_inplace<T> >(m, rhs);
    }
} 
// }}} end of implementation specific type-nested algorithms //

} } // namespace ambient::numeric

#undef ATOMIC
#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type 
#endif
