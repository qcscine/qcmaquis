#ifndef __MAQUIS_TYPES_P_DENSE_MATRIX_ALGORITHMS_HPP__
#define __MAQUIS_TYPES_P_DENSE_MATRIX_ALGORITHMS_HPP__

#include "types/p_dense_matrix/kernels/l_kernels.hpp"
#include "types/p_dense_matrix/kernels/c_kernels.hpp"

//#define AMBIENT_SERIAL_CHECK

#ifdef AMBIENT_SERIAL_CHECK
#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/utils/bindings.hpp"
#endif

#define size_type       typename p_dense_matrix<T>::size_type
#define scalar_type     typename p_dense_matrix<T>::scalar_type
#define difference_type typename p_dense_matrix<T>::difference_type

namespace maquis { namespace types {

    // {{{ p_dense_matrix advanced algorithms

    template <typename T>
    inline void reshape_r2l(p_dense_matrix<T>& left, const p_dense_matrix<T>& right,
                            size_t left_offset, size_t right_offset, 
                            size_t sdim, size_t ldim, size_t rdim)
    { // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sl = maquis::traits::matrix_cast<dense_matrix<T> >(left);
        dense_matrix<T> sr = maquis::traits::matrix_cast<dense_matrix<T> >(right);
        reshape_r2l(sl, sr, left_offset, right_offset, sdim, ldim, rdim);
#endif
        ambient::push(ambient::reshape_r2l_l<T>, ambient::reshape_r2l_c<T>, left, right, 
                      left_offset, right_offset, sdim, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
        if(sl == left){} else printf("--------------------- RESHAPE R2L WAS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void reshape_l2r(const p_dense_matrix<T>& left, p_dense_matrix<T>& right,
                            size_t left_offset, size_t right_offset, 
                            size_t sdim, size_t ldim, size_t rdim)
    { // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sl = maquis::traits::matrix_cast<dense_matrix<T> >(left);
        dense_matrix<T> sr = maquis::traits::matrix_cast<dense_matrix<T> >(right);
        reshape_l2r(sl, sr, left_offset, right_offset, sdim, ldim, rdim);
#endif
        ambient::push(ambient::reshape_l2r_l<T>, ambient::reshape_l2r_c<T>, left, right, 
                      left_offset, right_offset, sdim, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
        if(sr == right){} else printf("--------------------- RESHAPE L2R WAS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void lb_tensor_mpo(p_dense_matrix<T>& out, const p_dense_matrix<T>& in, const p_dense_matrix<T>& alfa,
                              size_t out_offset, size_t in_offset, 
                              size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    { // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sout = maquis::traits::matrix_cast<dense_matrix<T> >(out);
        dense_matrix<T> sin = maquis::traits::matrix_cast<dense_matrix<T> >(in);
        dense_matrix<T> salfa = maquis::traits::matrix_cast<dense_matrix<T> >(alfa);
        lb_tensor_mpo(sout, sin, salfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#endif
        ambient::push(ambient::lb_tensor_mpo_l<T>, ambient::lb_tensor_mpo_c<T>,
                      out, in, alfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
        if(sout == out){} else printf("--------------------- LB TENSOR MPO WAS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void rb_tensor_mpo(p_dense_matrix<T>& out, const p_dense_matrix<T>& in, const p_dense_matrix<T>& alfa,
                              size_t out_offset, size_t in_offset, 
                              size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    { // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sout = maquis::traits::matrix_cast<dense_matrix<T> >(out);
        dense_matrix<T> sin = maquis::traits::matrix_cast<dense_matrix<T> >(in);
        dense_matrix<T> salfa = maquis::traits::matrix_cast<dense_matrix<T> >(alfa);
        rb_tensor_mpo(sout, sin, salfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#endif
        ambient::push(ambient::rb_tensor_mpo_l<T>, ambient::rb_tensor_mpo_c<T>,
                      out, in, alfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
        if(sout == out){} else printf("--------------------- RB TENSOR MPO WAS INCORRECT!\n");
#endif
    }
     
    template <typename T>
    inline void scalar_norm(const p_dense_matrix<T>& a, scalar_type& ret){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> s = maquis::traits::matrix_cast<dense_matrix<T> >(a);
        typename dense_matrix<T>::value_type sret((T)ret);
        scalar_norm(s, sret);
#endif
        size_t m = num_rows(a);
        size_t n = num_cols(a);
        ambient::push(ambient::scalar_norm_l<T>, ambient::scalar_norm_c<T>, a, m, n, ret);
#ifdef AMBIENT_SERIAL_CHECK
        if(std::abs((T)ret - sret) > 0.01) printf("--------------------- SCALAR NORM IS INCORRECT (%.2f vs %.2f)\n", (T)ret, sret); 
#endif
    }

    template <typename T>
    inline void scalar_norm(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b, scalar_type& ret){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sa = maquis::traits::matrix_cast<dense_matrix<T> >(a);
        dense_matrix<T> sb = maquis::traits::matrix_cast<dense_matrix<T> >(b);
        T sret((T)ret); scalar_norm(sa, sb, sret);
#endif
        size_t m = num_rows(a);
        size_t n = num_cols(a);
        ambient::push(ambient::scalar_overlap_l<T>, ambient::scalar_overlap_c<T>, a, b, m, n, ret);
#ifdef AMBIENT_SERIAL_CHECK
        if(std::abs((T)ret - sret) > 0.01) printf("--------------------- SCALAR NORM OVERLAP IS INCORRECT (%.2f vs %.2f)\n", (T)ret, sret); 
#endif
    }

    template <typename T>
    inline void bond_renyi_entropies(const p_diagonal_matrix<T>& m, typename associated_real_vector<p_dense_matrix<T> >::type& sv){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        diagonal_matrix<T> sm = maquis::traits::matrix_cast<diagonal_matrix<T> >(m);
        typename associated_real_vector<p_dense_matrix<T> >::type ssv(sv);
        bond_renyi_entropies(sm, ssv);
#endif
        std::vector<T>* sc_ptr = &sv;
        ambient::push(ambient::push_back_sqr_gt_l<T>, ambient::push_back_sqr_gt_c<T>, m, sc_ptr);
#ifdef AMBIENT_SERIAL_CHECK
        if(sv != ssv) printf("--------------------- BOND RENYI ENTROPIES WAS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void left_right_boundary_init(p_dense_matrix<T> & a){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sa = maquis::traits::matrix_cast<dense_matrix<T> >(a);
        left_right_boundary_init(sa);
#endif
        a.fill_value(1.0);
#ifdef AMBIENT_SERIAL_CHECK
        if(sa == a){}else printf("--------------------- INCORRECT LEFT RIGHT BOUNDARY INIT!\n");
#endif
    }

    // }}}

    // {{{ p_dense_matrix operators (free functions)
    template <typename T>
    inline p_dense_matrix<T> operator + (p_dense_matrix<T> lhs, const p_dense_matrix<T>& rhs){
        return (lhs += rhs); 
    }

    template <typename T>
    inline p_dense_matrix<T> operator - (p_dense_matrix<T> lhs, const p_dense_matrix<T>& rhs){ 
        return (lhs -= rhs); 
    }

    template<typename T>
    inline const p_dense_matrix<T> operator * (p_dense_matrix<T> lhs, const p_dense_matrix<T>& rhs){ 
        return (lhs *= rhs); 
    }

    template<typename T, typename T2>
    inline const p_dense_matrix<T> operator * (p_dense_matrix<T> lhs, const T2& rhs){ 
        return (lhs *= rhs); 
    }

    template<typename T, typename T2>
    inline const p_dense_matrix<T> operator * (const T2& lhs, p_dense_matrix<T> rhs){ 
        return (rhs *= lhs); 
    }

    template <typename T>
    std::ostream& operator << (std::ostream& o, p_dense_matrix<T> const& m){
        for(size_type i=0; i< m.num_rows(); ++i){
            for(size_type j=0; j < m.num_cols(); ++j)
                maquis::cout << m(i,j) << " ";
            maquis::cout << std::endl;
        }
        return o;
    }

    template<typename T>
    inline size_type num_rows(const p_dense_matrix<T>& m){
        return m.num_rows();
    }

    template<typename T>
    inline size_type num_cols(const p_dense_matrix<T>& m){
        return m.num_cols();
    }

    template<typename T>
    inline void resize(p_dense_matrix<T>& m, size_t rows, size_t cols){
        m.resize(rows, cols);
    }

    template<typename T>
    inline scalar_type trace(const p_dense_matrix<T>& m){
        return m.trace();
    }

    template<typename T>
    inline p_dense_matrix<T> transpose(const p_dense_matrix<T>& m){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sm = maquis::traits::matrix_cast<dense_matrix<T> >(m);
#endif
        p_dense_matrix<T> t(m.num_cols(), m.num_rows());
        ambient::push(ambient::transpose_out_l<T>, ambient::transpose_out_c<T>, m, t);
        // p_dense_matrix<T> t(m); // alternative
        // t.transpose();
#ifdef AMBIENT_SERIAL_CHECK
        if(transpose(sm) == t){}else printf("--------------------- INCORRECT TRANSPOSE!\n");
#endif
        return t;
    }

    template<typename T>
    inline p_dense_matrix<T> conjugate(p_dense_matrix<T> m){
        m.conjugate();
        return m;
    }

    template<typename T>
    p_dense_matrix<T> exp(p_dense_matrix<T> m, T const & alfa = 1.){
        //printf("exp\n");
        typename associated_real_diagonal_matrix< p_dense_matrix<T> >::type evals(m.num_rows());
        p_dense_matrix<T> evecs;
        heev(m, evecs, evals);
        p_dense_matrix<T> e = evecs * exp(evals, alfa);
        e *= conjugate(transpose(evecs));
        return e;
    }

    // {{{ strassen matrix multiplication algorithm

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
    // }}}

    template<typename T>
    inline void gemm(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c){
        // gs
        if(num_cols(a) != num_rows(b)) maquis::cout << "Incorrect GEMM operation\n";
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sa = maquis::traits::matrix_cast<dense_matrix<T> >(a);
        dense_matrix<T> sb = maquis::traits::matrix_cast<dense_matrix<T> >(b);
        dense_matrix<T> sc = maquis::traits::matrix_cast<dense_matrix<T> >(c);
        gemm(sa,sb,sc);
#endif
        c.resize(a.num_rows(), b.num_cols());
        ambient::push(ambient::gemm_l<T>, ambient::gemm_c<T>, a, b, c);
#ifdef AMBIENT_SERIAL_CHECK
        if(sc == c){} else printf("--------------------- GEMM WAS INCORRECT!\n");
#endif
    }

    template<typename T, typename D>
    inline void gemm(const p_dense_matrix<T>& a, const p_diagonal_matrix<D>& b, p_dense_matrix<T>& c){
        // gs
        if(num_cols(a) != num_rows(b)) maquis::cout << "Incorrect GEMM operation\n";
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sa = maquis::traits::matrix_cast<dense_matrix<T> >(a);
        diagonal_matrix<D> sb = maquis::traits::matrix_cast<diagonal_matrix<D> >(b);
        dense_matrix<T> sc = maquis::traits::matrix_cast<dense_matrix<T> >(c);
        gemm(sa, sb, sc);
#endif
        size_t m = a.num_rows();
        size_t n = b.num_cols();
        size_t k = a.num_cols();
        c.resize(m, n);
        ambient::push(ambient::gemm_diagonal_rhs_l<T,D>, ambient::gemm_diagonal_rhs_c<T,D>, a, b, c, m, n, k);
#ifdef AMBIENT_SERIAL_CHECK
        if(sc == c){} else printf("--------------------- GEMM PDP WAS INCORRECT!\n");
#endif
    }

    template<typename T, typename D>
    inline void gemm(const p_diagonal_matrix<D>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c){
        // gs
        if(num_cols(a) != num_rows(b)) maquis::cout << "Incorrect GEMM operation\n";
#ifdef AMBIENT_SERIAL_CHECK
        diagonal_matrix<D> sa = maquis::traits::matrix_cast<diagonal_matrix<D> >(a);
        dense_matrix<T> sb = maquis::traits::matrix_cast<dense_matrix<T> >(b);
        dense_matrix<T> sc = maquis::traits::matrix_cast<dense_matrix<T> >(c);
        gemm(sa, sb, sc);
#endif
        size_t m = a.num_rows();
        size_t n = b.num_cols();
        size_t k = a.num_cols();
        c.resize(m, n);
        ambient::push(ambient::gemm_diagonal_lhs_l<T,D>, ambient::gemm_diagonal_lhs_c<T,D>, a, b, c, m, n, k);
#ifdef AMBIENT_SERIAL_CHECK
        if(sc == c){} else printf("--------------------- GEMM DPP WAS INCORRECT!\n");
#endif
    }

    template<typename T>
    inline void svd(const p_dense_matrix<T>& a, p_dense_matrix<T>& u, p_dense_matrix<T>& vt,
                    typename associated_real_diagonal_matrix<p_dense_matrix<T> >::type& s)
    { // gs
        int m = num_rows(a);
        int n = num_cols(a);
        int k = std::min(m,n);
#ifdef AMBIENT_SERIAL_CHECK // warning: we don't check the correctness (out can vary)
        dense_matrix<T> sa = maquis::traits::matrix_cast<dense_matrix<T> >(a);
        dense_matrix<T> su(m,k);
        dense_matrix<T> sv(k,n);
        typename associated_real_diagonal_matrix< dense_matrix<T> >::type ss(k,k);
        svd(sa,su,sv,ss);
        u  = maquis::traits::matrix_cast<p_dense_matrix<T> >(su);
        vt = maquis::traits::matrix_cast<p_dense_matrix<T> >(sv);
        s  = maquis::traits::matrix_cast<typename associated_real_diagonal_matrix< p_dense_matrix<T> >::type>(ss);
#else
        u.resize(m, k);
        vt.resize(k, n);
        s.resize(k, k);
        ambient::push(ambient::svd_l<T>, ambient::svd_c<T>, a, m, n, u, vt, s);
#endif
    }

    template<typename T>
    inline void heev(p_dense_matrix<T> a, p_dense_matrix<T>& evecs,
                     typename associated_real_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    { // gs
        if(num_rows(a) != num_cols(a) || num_rows(evals) != num_rows(a)) maquis::cout << "Incorrect HEEV operation\n";
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sa = maquis::traits::matrix_cast<dense_matrix<T> >(a);
        dense_matrix<T> sevecs = maquis::traits::matrix_cast<dense_matrix<T> >(evecs);
        typename associated_real_diagonal_matrix< dense_matrix<T> >::type sevals = 
            maquis::traits::matrix_cast<typename associated_real_diagonal_matrix< dense_matrix<T> >::type >(evals);
        heev(sa, sevecs, sevals.get_values());
#endif
        size_t m = num_rows(a);
        evecs.resize(m, m);
        ambient::push(ambient::heev_l<double>, ambient::heev_c<double>, a, m, evals); // destoys U triangle of M
        evecs = a;
#ifdef AMBIENT_SERIAL_CHECK
        if(sevecs == evecs){}else printf("--------------------- HEEV WAS INCORRECT!\n");
#endif
    }

    template<typename T>
    inline void syev(p_dense_matrix<T> a, p_dense_matrix<T>& evecs,
                     typename associated_real_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    {
        heev(a, evecs, evals);
    }

    template<typename T>
    inline void qr(p_dense_matrix<T> m, p_dense_matrix<T> & q, p_dense_matrix<T> & r){
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- QR");
    }

    template<typename T, class G>
    inline void generate(p_dense_matrix<T>& m, G g){ // warning: G isn't used
        m.fill_random();
    }
    // }}}

// {{{ implementation-specific type-nested algorithms //
namespace algorithms {

    template<typename T>
    inline void resize(p_dense_matrix_impl<T>& m, size_type rows, size_type cols, size_type orows, size_type ocols){ 
        // gs
        ambient::push(ambient::resize_l<T>, ambient::resize_c<T>, m, rows, cols, orows, ocols);
    }

    template<typename T>
    inline void remove_rows(p_dense_matrix_impl<T>& m, size_type i, difference_type k){
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- REMOVE ROWS");
        ambient::push(ambient::remove_rows_l<T>, ambient::remove_rows_c<T>, m, i, k);
    }

    template<typename T>
    inline void remove_cols(p_dense_matrix_impl<T>& m, size_type j, difference_type k){
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- REMOVE COLS");
        ambient::push(ambient::remove_cols_l<T>, ambient::remove_cols_c<T>, m, j, k);
    }

    template<typename T>
    inline void fill_identity(p_dense_matrix_impl<T>& a){
        size_t m = a.num_rows();
        size_t n = a.num_cols();
        ambient::push(ambient::init_identity_l<T>, ambient::init_identity_c<T>, a, m, n);
    }

    template<typename T>
    inline void fill_random(p_dense_matrix_impl<T>& a){
        size_t m = a.num_rows();
        size_t n = a.num_cols();
        ambient::push(ambient::init_random_l<T>, ambient::init_random_c<T>, a, m, n);
    }

    template<typename T>
    inline void fill_value(p_dense_matrix_impl<T>& a, T value){
        if(value == 0.) return; // matrices are 0s by default
        size_t m = a.num_rows();
        size_t n = a.num_cols();
        ambient::push(ambient::init_value_l<T>, ambient::init_value_c<T>, a, m, n, value);
    }

    template<typename T>
    inline void inplace_conjugate(p_dense_matrix_impl<T>& m){
        // gs (doubles)
        //printf("inplace conjugate\n");
        // TODO: does nothing for now
    }

    template<typename T>
    inline void inplace_transpose(p_dense_matrix_impl<T>& m){
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- INPLACE TRANSPOSE");
        assert(m.num_rows() == m.num_cols()); // current limitation
        ambient::push(ambient::transpose_l<T>, ambient::transpose_c<T>, m);
    }

    template <typename T>
    inline scalar_type trace(const p_dense_matrix_impl<T>& a) {
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sm = maquis::traits::matrix_cast<dense_matrix<T> >(a);
#endif
        scalar_type trace;
        size_t n = std::min(a.num_rows(), a.num_cols());
        ambient::push(ambient::trace_l<T>, ambient::trace_c<T>, a, n, trace);
#ifdef AMBIENT_SERIAL_CHECK
        if(maquis::types::trace(sm) != (T)trace) printf("--------------------- TRACE IS INCORRECT!\n");
#endif
        return trace;
    }

    template <typename T>
    inline void add_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs) {
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sm = maquis::traits::matrix_cast<dense_matrix<T> >(m);
        dense_matrix<T> srhs = maquis::traits::matrix_cast<dense_matrix<T> >(rhs);
        sm += srhs;
#endif
        ambient::push(ambient::mem_bound_l<T>, ambient::add_c<T>, m, rhs);
#ifdef AMBIENT_SERIAL_CHECK
        if(sm == m){}else printf("--------------------- ADD INPLACE IS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void sub_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs) {
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> sm = maquis::traits::matrix_cast<dense_matrix<T> >(m);
        dense_matrix<T> srhs = maquis::traits::matrix_cast<dense_matrix<T> >(rhs);
        sm -= srhs;
#endif
        ambient::push(ambient::mem_bound_l<T>, ambient::sub_c<T>, m, rhs);
#ifdef AMBIENT_SERIAL_CHECK
        if(sm == m){}else printf("--------------------- SUB INPLACE IS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void gemm_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs) {
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- GEMM INPLACE");
        ambient::push(ambient::gemm_inplace_l<T>, ambient::gemm_inplace_c<T>, m, rhs);
    }

    template <typename T>
    inline void gemm_diag_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs_diag) {
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- GEMM DIAG INPLACE");
        //ambient::push(ambient::gemm_diag_inplace_l<T>, ambient::gemm_diag_inplace_c<T>, m, rhs);
    }

    template <typename T>
    inline void scale_inplace(p_dense_matrix_impl<T>& a, const scalar_type& rhs) {
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        dense_matrix<T> s = maquis::traits::matrix_cast<dense_matrix<T> >(a);
        s *= (T)rhs;
#endif
        size_t m = a.num_rows();
        size_t n = a.num_cols();
        ambient::push(ambient::scale_l<T>, ambient::scale_c<T>, a, m, n, rhs);
#ifdef AMBIENT_SERIAL_CHECK
        if(s == a){} else printf("--------------------- SCALE WAS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void cpy(p_dense_matrix_impl<T>& dst, const p_dense_matrix_impl<T>& src){
        // gs
        ambient::push(static_cast<typename p_dense_matrix_impl<T>::copy_t>(&ambient::copy_l), 
                      static_cast<typename p_dense_matrix_impl<T>::copy_t>(&ambient::copy_c), dst, src);
    }
} 
// }}} end of implementation specific type-nested algorithms //

} } // namespace maquis::types

#undef size_type
#undef scalar_type
#undef difference_type 
#endif
