#ifndef __MAQUIS_TYPES_P_DENSE_MATRIX_ALGORITHMS_HPP__
#define __MAQUIS_TYPES_P_DENSE_MATRIX_ALGORITHMS_HPP__

#include "types/p_dense_matrix/kernels/r_kernels.hpp"
#include "types/p_dense_matrix/kernels/l_kernels.hpp"
#include "types/p_dense_matrix/kernels/c_kernels.hpp"

namespace maquis { namespace types { namespace algorithms {

    #define size_type typename p_dense_matrix_impl<T>::size_type
    #define difference_type typename p_dense_matrix_impl<T>::difference_type

// {{{ implementation-specific type-nested algorithms //

    template<typename T>
    void clear(p_dense_matrix_impl<T>& m){
        ambient::push(ambient::clear_l, ambient::clear_c, m);
    }

    template<typename T>
    void resize(p_dense_matrix_impl<T>& m, size_type rows, size_type cols){
        ambient::push(ambient::resize_l<T>, ambient::resize_c<T>, m, rows, cols);
    }

    template<typename T>
    void remove_rows(p_dense_matrix_impl<T>& m, size_type i, difference_type k){
        ambient::push(ambient::remove_rows_l<T>, ambient::remove_rows_c<T>, m, i, k);
    }

    template<typename T>
    void remove_cols(p_dense_matrix_impl<T>& m, size_type j, difference_type k){
        ambient::push(ambient::remove_cols_l<T>, ambient::remove_cols_c<T>, m, j, k);
    }

    template<typename T>
    void inplace_conjugate(p_dense_matrix_impl<T>& m){
        // TODO: does nothing for now
    }

    template<typename T>
    void inplace_transpose(p_dense_matrix_impl<T>& m){
        ambient::push(ambient::transpose_l<T>, ambient::transpose_c<T>, m);
    }

    template <typename T>
    T trace(p_dense_matrix_impl<T>& m) {
        T trace; // stack memory
        T* out = &trace;
        ambient::push(ambient::trace_l<T>, ambient::trace_c<T>, m, out);
        ambient::playout();
        return trace;
    }

    template <typename T>
    void add_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs) {
        ambient::push(ambient::mem_bound_l<T>, ambient::add_c<T>, m, rhs);
    }

    template <typename T>
    void sub_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs) {
        ambient::push(ambient::mem_bound_l<T>, ambient::sub_c<T>, m, rhs);
    }

    template <typename T>
    void gemm_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs) {
        ambient::push(ambient::gemm_inplace_l<T>, ambient::gemm_inplace_c<T>, m, rhs);
    }

    template <typename T, typename S>
    void scale_inplace(p_dense_matrix_impl<T>& m, const S& rhs) {
        ambient::push(ambient::scale_l<T>, ambient::scale_c<T>, m, rhs);
    }

    template <typename T>
    void cpy(p_dense_matrix_impl<T>& dst, const p_dense_matrix_impl<T>& src){
        ambient::push(ambient::copy_l<p_dense_matrix_impl<T> >, ambient::copy_c<p_dense_matrix_impl<T> >, dst, src);
    }

// }}} end of implementation specific type-nested algorithms //

    template<typename T>
    p_dense_matrix<T> conjugate(p_dense_matrix<T> m){
        m.conjugate();
        return m;
    }

    template<typename T>
    p_dense_matrix<T> transpose(const p_dense_matrix<T> m){
        m.transpose();
        return m;
    }


    // {{{ strassen multiplication supplementary functions

    template<typename T>
    void gemm_strassen_gad(const p_dense_matrix<T>& a, size_t ai, size_t aj, 
                           p_dense_matrix<T>& r, size_t n)
    {
        ambient::push(ambient::gemm_strassen_gad_l<T>,  // output reduced matrix
                      ambient::gemm_strassen_gad_c<T>,  // a11 + a12,  a12 - a22
                      a, ai, aj, r, n);                 // a21 - a11,  a22 + a21
    }

    template<typename T>
    void gemm_strassen_dad(const p_dense_matrix<T>& a, size_t ai, size_t aj,
                           const p_dense_matrix<T>& b, size_t bi, size_t bj, 
                           p_dense_matrix<T>& r, size_t n)
    {
        ambient::push(ambient::gemm_strassen_dad_l<T>,  // output reduced matrix
                      ambient::gemm_strassen_dad_c<T>,  // a11 + a22,  b11 + b22
                      a, ai, aj, b, bi, bj, r, n);

    }

    template<typename T>
    void gemm_strassen_pluseqs(const p_dense_matrix<T>& a, size_t ai, size_t aj,
                               const p_dense_matrix<T>& b, size_t bi, size_t bj, 
                                     p_dense_matrix<T>& c, size_t ci, size_t cj, 
                                     size_t n)
    {
        ambient::push(ambient::add_sum_submx_l<T>, 
                      ambient::add_sum_submx_c<T>, 
                      a, ai, aj, b, bi, bj, c, ci, cj, n);
    }

    template<typename T>
    void gemm_strassen_pluseqd(const p_dense_matrix<T>& a, size_t ai, size_t aj,
                               const p_dense_matrix<T>& b, size_t bi, size_t bj, 
                                     p_dense_matrix<T>& c, size_t ci, size_t cj, 
                                     size_t n)
    {
        ambient::push(ambient::add_dif_submx_l<T>, 
                      ambient::add_dif_submx_c<T>, 
                      a, ai, aj, b, bi, bj, c, ci, cj, n);
    }

    // }}}

    template<typename T>
    void gemm_strassen(const p_dense_matrix<T>& a, size_t ai, size_t aj,
                       const p_dense_matrix<T>& b, size_t bi, size_t bj,
                             p_dense_matrix<T>& c, size_t ci, size_t cj, 
                             size_t n)
    {

        if(n > 128){ 
            //printf("Performing strassen with size %d\n", n);
            // strassen algorithm recursion
            p_dense_matrix<T> m1(n/2, n/2);
            p_dense_matrix<T> m2(n/2, n/2);
            p_dense_matrix<T> m3(n/2, n/2);
            p_dense_matrix<T> m4(n/2, n/2);
            p_dense_matrix<T> m5(n/2, n/2);
            p_dense_matrix<T> m6(n/2, n/2);
            p_dense_matrix<T> m7(n/2, n/2);

            p_dense_matrix<T> ar(n,n);
            p_dense_matrix<T> br(n,n);
            p_dense_matrix<T> dr(n/2,n);

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
    void gemm_strassen(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c){
        size_t n  = num_cols(a);
        assert(n == num_rows(a));
        assert(n == num_cols(b));
        assert(n == num_rows(b));
        c.resize(n, n);
        gemm_strassen(a, 0, 0, b, 0, 0, c, 0, 0, n);
    }

    template<typename T>
    void gemm(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c){
        assert(num_cols(a) == num_rows(b));
        c.resize(a.num_rows(), b.num_cols());
        ambient::push(ambient::gemm_l<T>, ambient::gemm_c<T>, a, b, c);
    }

    template<typename T>
    void gemm(const p_dense_matrix<T>& a, const p_diagonal_matrix<T>& b, p_dense_matrix<T>& c){
        assert(num_cols(a) == num_rows(b));
        c.resize(a.num_rows(), b.num_cols());
        ambient::push(ambient::gemm_diagonal_rhs_l<T>, ambient::gemm_diagonal_rhs_c<T>, a, b, c);
    }

    template<typename T>
    void gemm(const p_diagonal_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c){
        assert(num_cols(a) == num_rows(b));
        c.resize(a.num_rows(), b.num_cols());
        ambient::push(ambient::gemm_diagonal_lhs_l<T>, ambient::gemm_diagonal_lhs_c<T>, a, b, c);
    }

    template<typename T>
    void svd(const p_dense_matrix<T>& a, p_dense_matrix<T>& u, p_dense_matrix<T>& vt,
             typename associated_diagonal_matrix<p_dense_matrix<double> >::type& s)
    {
        int m = num_rows(a);
        int n = num_cols(a);
        int k = std::min(m,n);
        u.resize(m, k);
        vt.resize(k, n);
        s.resize(k, k);
        ambient::push(ambient::svd_l<T>, ambient::svd_c<T>, a, m, n, u, vt, s);
    }

    template<typename T>
    void heev(p_dense_matrix<T> a, p_dense_matrix<T>& evecs,
              typename associated_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    {
        assert(num_rows(a) == num_cols(a));
        assert(num_rows(evals) == num_rows(a));
        int m = num_rows(a);
        evecs.resize(m, m);
        ambient::push(ambient::heev_l, ambient::heev_c, a, m, evals); // destoys U triangle of M
        evecs = a;
    }

    template<typename T>
    void syev(p_dense_matrix<T> a, p_dense_matrix<T>& evecs,
              typename associated_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    {
        heev(a, evecs, evals);
    }

    template<typename T>
    void qr(p_dense_matrix<T> m, p_dense_matrix<T> & q, p_dense_matrix<T> & r){
        assert(false);
    }

    template<typename T, class G>
    void generate(p_dense_matrix<T>& m, G g){
        m.generate();
    }

    // garbage: std replacements //
    template<typename T>
    void copy(typename associated_vector<T>::type& sc, typename associated_diagonal_matrix<T>::type& s){
    // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector 
        ambient::push(ambient::associated_copy_l<T>, ambient::copy_c<T>, sc, s);
    }

    template<typename T>
    void copy_sqr_gt(std::vector<typename T::value_type>& sc, typename associated_diagonal_matrix<T>::type& s, const double& prec){
    // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector
        std::vector<typename T::value_type>* sc_ptr = &sc;
        ambient::push(ambient::push_back_sqr_gt_l, ambient::push_back_sqr_gt_c, sc_ptr, s, prec);
    }

    template<typename T>
    void copy_after(std::vector<typename T::value_type>& sc, const size_type pos, typename associated_diagonal_matrix<T>::type& s){
    // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector
        std::vector<typename T::value_type>* sc_ptr = &sc;  
        ambient::push(ambient::copy_after_std_l, ambient::copy_after_std_c, sc_ptr, pos, s);
    }

    template<typename T>
    void copy_after(typename associated_vector<T>::type& sc, const size_type pos, typename associated_diagonal_matrix<T>::type& s){
    // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector 
        ambient::push(ambient::copy_after_l<T>, ambient::copy_after_c<T>, sc, pos, s);
    }

    #undef size_type
    #undef difference_type 

} } } // namespace maquis::types::algorithms

#endif
