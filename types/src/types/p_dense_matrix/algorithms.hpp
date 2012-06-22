#ifndef __MAQUIS_TYPES_P_DENSE_MATRIX_ALGORITHMS_HPP__
#define __MAQUIS_TYPES_P_DENSE_MATRIX_ALGORITHMS_HPP__

#include "types/p_dense_matrix/kernels/utils.hpp"
#include "types/p_dense_matrix/kernels/kernels.hpp"
#include "types/p_dense_matrix/kernels/atomics.hpp"

//#define AMBIENT_SERIAL_CHECK
#define USE_ATOMIC(condition, kernel, ...) assert(condition); ambient::push< ambient::kernel ## _atomic<T> >(__VA_ARGS__);

#ifdef AMBIENT_SERIAL_CHECK
#include "alps/numeric/matrix/matrix.hpp"
#include "alps/numeric/matrix/matrix_blas.hpp"
#include "alps/numeric/matrix/matrix_interface.hpp"
#include "alps/numeric/matrix/resizable_matrix_interface.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "types/utils/bindings.hpp"
#endif

#define size_type       typename p_dense_matrix<T>::size_type
#define scalar_type     typename p_dense_matrix<T>::scalar_type
#define difference_type typename p_dense_matrix<T>::difference_type

namespace maquis { namespace types {

    // {{{ p_dense_matrix advanced algorithms

    template <typename T>
    inline void reshape_l2r(const p_dense_matrix<T>& left, p_dense_matrix<T>& right,
                            size_t left_offset, size_t right_offset, 
                            size_t sdim, size_t ldim, size_t rdim)
    { // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sl = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(left);
        alps::numeric::matrix<T> sr = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(right);
        reshape_l2r(sl, sr, left_offset, right_offset, sdim, ldim, rdim);
#endif
        USE_ATOMIC(left.atomic() && right.atomic(), reshape_l2r, 
                   left, right, left_offset, right_offset, sdim, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
        if(sr == right){} else printf("--------------------- RESHAPE L2R WAS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void reshape_r2l(p_dense_matrix<T>& left, const p_dense_matrix<T>& right,
                            size_t left_offset, size_t right_offset, 
                            size_t sdim, size_t ldim, size_t rdim)
    { // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sl = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(left);
        alps::numeric::matrix<T> sr = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(right);
        reshape_r2l(sl, sr, left_offset, right_offset, sdim, ldim, rdim);
#endif
        USE_ATOMIC(left.atomic() && right.atomic(), reshape_r2l, 
                   left, right, left_offset, right_offset, sdim, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
        if(sl == left){} else printf("--------------------- RESHAPE R2L WAS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void lb_tensor_mpo(p_dense_matrix<T>& out, const p_dense_matrix<T>& in, const p_dense_matrix<T>& alfa,
                              size_t out_offset, size_t in_offset, 
                              size_t sdim1, size_t sdim2, size_t ldim, size_t rdim)
    { // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sout = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(out);
        alps::numeric::matrix<T> sin = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(in);
        alps::numeric::matrix<T> salfa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(alfa);
        lb_tensor_mpo(sout, sin, salfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#endif
        USE_ATOMIC(out.atomic() && in.atomic() && alfa.atomic(), lb_tensor_mpo, 
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
        alps::numeric::matrix<T> sout = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(out);
        alps::numeric::matrix<T> sin = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(in);
        alps::numeric::matrix<T> salfa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(alfa);
        rb_tensor_mpo(sout, sin, salfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#endif
        USE_ATOMIC(out.atomic() && in.atomic() && alfa.atomic(), rb_tensor_mpo, 
                   out, in, alfa, out_offset, in_offset, sdim1, sdim2, ldim, rdim);
#ifdef AMBIENT_SERIAL_CHECK
        if(sout == out){} else printf("--------------------- RB TENSOR MPO WAS INCORRECT!\n");
#endif
    }
     
    template <typename T>
    inline void scalar_norm(const p_dense_matrix<T>& a, scalar_type& ret){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> s = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
        typename alps::numeric::matrix<T>::value_type sret((T)ret);
        scalar_norm(s, sret);
#endif
        USE_ATOMIC(a.atomic(), scalar_norm, a, ret);
#ifdef AMBIENT_SERIAL_CHECK
        if(std::abs((T)ret - sret) > 0.01) printf("--------------------- SCALAR NORM IS INCORRECT (%.2f vs %.2f)\n", (T)ret, sret); 
#endif
    }

    template <typename T>
    inline void scalar_norm(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b, scalar_type& ret){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
        alps::numeric::matrix<T> sb = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(b);
        T sret((T)ret); scalar_norm(sa, sb, sret);
#endif
        USE_ATOMIC(a.atomic(), scalar_overlap, a, b, ret);
#ifdef AMBIENT_SERIAL_CHECK
        if(std::abs((T)ret - sret) > 0.01) printf("--------------------- SCALAR NORM OVERLAP IS INCORRECT (%.2f vs %.2f)\n", (T)ret, sret); 
#endif
    }

    template <typename T>
    inline void bond_renyi_entropies(const p_diagonal_matrix<T>& m, typename alps::numeric::associated_real_vector<p_dense_matrix<T> >::type& sv){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        diagonal_matrix<T> sm = maquis::traits::matrix_cast<diagonal_matrix<T> >(m);
        typename alps::numeric::associated_real_vector<p_dense_matrix<T> >::type ssv(sv);
        bond_renyi_entropies(sm, ssv);
#endif
        std::vector<T>* sc_ptr = &sv;
        ambient::push< ambient::push_back_sqr_gt<T> >(m, sc_ptr);
        ambient::playout();
#ifdef AMBIENT_SERIAL_CHECK
        if(sv != ssv) printf("--------------------- BOND RENYI ENTROPIES WAS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void left_right_boundary_init(p_dense_matrix<T>& a){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
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
    inline p_dense_matrix<T> operator + (const p_dense_matrix<T>& lhs, const p_dense_matrix<T>& rhs){
        return (lhs.clone() += rhs); 
    }

    template <typename T>
    inline p_dense_matrix<T> operator - (const p_dense_matrix<T>& lhs, const p_dense_matrix<T>& rhs){ 
        return (lhs.clone() -= rhs); 
    }

    template<typename T>
    inline const p_dense_matrix<T> operator * (const p_dense_matrix<T>& lhs, const p_dense_matrix<T>& rhs){ 
        return (lhs.clone() *= rhs); 
    }

    template<typename T, typename T2>
    inline const p_dense_matrix<T> operator * (const p_dense_matrix<T>& lhs, const T2& rhs){ 
        return (lhs.clone() *= rhs); 
    }

    template<typename T, typename T2>
    inline const p_dense_matrix<T> operator * (const T2& lhs, const p_dense_matrix<T>& rhs){ 
        return (rhs.clone() *= lhs); 
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
    inline const p_dense_matrix<T>& transpose(const p_dense_matrix<T>& a){
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sm = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
#endif
        p_dense_matrix<T> t(a.num_cols(), a.num_rows());
        USE_ATOMIC(a.atomic(), transpose_out, a, t);
        const_cast<p_dense_matrix<T>&>(a).swap(t); // somewhat inplace :)
        // p_dense_matrix<T> t(m); // alternative
        // t.transpose();
#ifdef AMBIENT_SERIAL_CHECK
        if(transpose(sm) == t){}else printf("--------------------- INCORRECT TRANSPOSE!\n");
#endif
        return a;
    }

    template<typename T>
    inline const p_dense_matrix<T>& conjugate(const p_dense_matrix<T>& m){
        //m.conjugate();
        return m;
    }

    template<typename T>
    inline p_dense_matrix<T> exp(const p_dense_matrix<T>& m, const T& alfa = 1.){
        printf("----------------------------------- exp is actually used!\n\n");
        typename alps::numeric::associated_real_diagonal_matrix< p_dense_matrix<T> >::type evals(m.num_rows());
        p_dense_matrix<T> evecs;
        heev(m, evecs, evals);
        return (evecs * exp(evals, alfa))*conjugate(transpose(evecs));
    }

    // {{{ strassen matrix multiplication algorithm

    // {{{ strassen multiplication supplementary functions
    /*
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
    }*/
    // }}}

    template<typename T>
    inline void gemm(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c){
        // gs
        assert(num_cols(a) == num_rows(b));
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
        alps::numeric::matrix<T> sb = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(b);
        alps::numeric::matrix<T> sc = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(c);
        gemm(sa,sb,sc);
#endif
        p_dense_matrix<T> result(a.num_rows(), b.num_cols());
        USE_ATOMIC(a.atomic() && b.atomic(), gemm_general, a, b, result);
        c.swap(result);
#ifdef AMBIENT_SERIAL_CHECK
        if(sc == c){} else printf("--------------------- GEMM WAS INCORRECT!\n");
#endif
    }

    template<typename T, typename D>
    inline void gemm(const p_dense_matrix<T>& a, const p_diagonal_matrix<D>& b, p_dense_matrix<T>& c){
        // gs
        assert(num_cols(a) == num_rows(b));
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
        diagonal_matrix<D> sb = maquis::traits::matrix_cast<diagonal_matrix<D> >(b);
        alps::numeric::matrix<T> sc = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(c);
        gemm(sa, sb, sc);
#endif
        size_t m = a.num_rows();
        size_t n = b.num_cols();
        size_t k = a.num_cols();
        p_dense_matrix<T> result(m, n);
        ambient::push< ambient::gemm_diagonal_rhs<T,D> >(a, b, result, m, n, k);
        c.swap(result);
#ifdef AMBIENT_SERIAL_CHECK
        if(sc == c){} else printf("--------------------- GEMM PDP WAS INCORRECT!\n");
#endif
    }

    template<typename T, typename D>
    inline void gemm(const p_diagonal_matrix<D>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c){
        // gs
        assert(num_cols(a) == num_rows(b));
#ifdef AMBIENT_SERIAL_CHECK
        diagonal_matrix<D> sa = maquis::traits::matrix_cast<diagonal_matrix<D> >(a);
        alps::numeric::matrix<T> sb = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(b);
        alps::numeric::matrix<T> sc = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(c);
        gemm(sa, sb, sc);
#endif
        size_t m = a.num_rows();
        size_t n = b.num_cols();
        size_t k = a.num_cols();
        p_dense_matrix<T> result(m, n);
        ambient::push< ambient::gemm_diagonal_lhs<T,D> >(a, b, result, m, n, k);
        c.swap(result);
#ifdef AMBIENT_SERIAL_CHECK
        if(sc == c){} else printf("--------------------- GEMM DPP WAS INCORRECT!\n");
#endif
    }

    template<typename T>
    inline void svd(const p_dense_matrix<T>& a, p_dense_matrix<T>& u, p_dense_matrix<T>& vt,
                    typename alps::numeric::associated_real_diagonal_matrix<p_dense_matrix<T> >::type& s)
    { // gs
        int m = num_rows(a);
        int n = num_cols(a);
        int k = std::min(m,n);
#ifdef AMBIENT_SERIAL_CHECK // warning: we don't check the correctness (out can vary)
        alps::numeric::matrix<T> sa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
        alps::numeric::matrix<T> su(m,k);
        alps::numeric::matrix<T> sv(k,n);
        typename alps::numeric::associated_real_diagonal_matrix< alps::numeric::matrix<T> >::type ss(k,k);
        svd(sa,su,sv,ss);
        u  = maquis::traits::matrix_cast<p_dense_matrix<T> >(su);
        vt = maquis::traits::matrix_cast<p_dense_matrix<T> >(sv);
        s  = maquis::traits::matrix_cast<typename alps::numeric::associated_real_diagonal_matrix< p_dense_matrix<T> >::type>(ss);
#else
        u.resize(m, k);
        vt.resize(k, n);
        s.resize(k, k);
        p_dense_matrix<T> ac(a);
        USE_ATOMIC(a.atomic(), svd, ac, u, vt, s);
#endif
    }

    template<typename T>
    inline void heev(const p_dense_matrix<T>& a, p_dense_matrix<T>& evecs,
                     typename alps::numeric::associated_real_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    { // gs
        assert(num_rows(a) == num_cols(a) && num_rows(evals) == num_rows(a));
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sa = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
        alps::numeric::matrix<T> sevecs = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(evecs);
        typename alps::numeric::associated_real_diagonal_matrix< alps::numeric::matrix<T> >::type sevals = 
            maquis::traits::matrix_cast<typename alps::numeric::associated_real_diagonal_matrix< alps::numeric::matrix<T> >::type >(evals);
        heev(sa, sevecs, sevals.get_values());
#endif
        p_dense_matrix<T> ac(a);
        USE_ATOMIC(a.atomic(), heev, ac, evals); // destoys U triangle of M
        evecs.swap(ac);
#ifdef AMBIENT_SERIAL_CHECK
        if(sevecs == evecs){}else printf("--------------------- HEEV WAS INCORRECT!\n");
#endif
    }

    template<typename T>
    inline void syev(const p_dense_matrix<T>& a, p_dense_matrix<T>& evecs,
                     typename alps::numeric::associated_real_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    {
        heev(a, evecs, evals);
    }

    template<typename T>
    inline void qr(const p_dense_matrix<T>& m, p_dense_matrix<T>& q, p_dense_matrix<T>& r){
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
    inline void resize(p_dense_matrix_impl<T>& r, size_type rows, size_type cols, p_dense_matrix_impl<T>& m, size_type orows, size_type ocols){ 
        // gs
        USE_ATOMIC(m.atomic(), resize, r, m, std::min(rows, orows), std::min(cols, ocols));
    }

    template<typename T>
    inline void remove_rows(p_dense_matrix_impl<T>& m, size_type i, difference_type k){
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- REMOVE ROWS");
        ambient::push< ambient::remove_rows<T> >(m, i, k);
    }

    template<typename T>
    inline void remove_cols(p_dense_matrix_impl<T>& m, size_type j, difference_type k){
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- REMOVE COLS");
        ambient::push< ambient::remove_cols<T> >(m, j, k);
    }

    template<typename T>
    inline void fill_identity(p_dense_matrix_impl<T>& a){
        size_t m = a.num_rows();
        size_t n = a.num_cols();
        ambient::push< ambient::init_identity<T> >(a, m, n);
    }

    template<typename T>
    inline void fill_random(p_dense_matrix_impl<T>& a){
        size_t m = a.num_rows();
        size_t n = a.num_cols();
        ambient::push< ambient::init_random<T> >(a, m, n);
    }

    template<typename T>
    inline void fill_value(p_dense_matrix_impl<T>& a, T value){
        if(value == 0.) return; // matrices are 0s by default
        size_t m = a.num_rows();
        size_t n = a.num_cols();
        ambient::push< ambient::init_value<T> >(a, m, n, value);
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
        ambient::push< ambient::transpose<T> >(m);
    }

    template <typename T>
    inline scalar_type trace(const p_dense_matrix_impl<T>& a) {
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sm = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
#endif
        scalar_type trace;
        size_t n = std::min(a.num_rows(), a.num_cols());
        ambient::push< ambient::trace<T> >(a, n, trace);
#ifdef AMBIENT_SERIAL_CHECK
        if(trace(sm) != (T)trace) printf("--------------------- TRACE IS INCORRECT!\n");
#endif
        return trace;
    }

    template <typename T>
    inline void add_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs) {
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sm = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(m);
        alps::numeric::matrix<T> srhs = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(rhs);
        sm += srhs;
#endif
        USE_ATOMIC(m.atomic(), add, m, rhs);
#ifdef AMBIENT_SERIAL_CHECK
        if(sm == m){}else printf("--------------------- ADD INPLACE IS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void sub_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs) {
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> sm = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(m);
        alps::numeric::matrix<T> srhs = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(rhs);
        sm -= srhs;
#endif
        USE_ATOMIC(m.atomic(), sub, m, rhs);
#ifdef AMBIENT_SERIAL_CHECK
        if(sm == m){}else printf("--------------------- SUB INPLACE IS INCORRECT!\n");
#endif
    }

    template <typename T>
    inline void gemm_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs) {
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- GEMM INPLACE");
        ambient::push< ambient::gemm_inplace<T> >(m, rhs);
    }

    template <typename T>
    inline void gemm_diag_inplace(p_dense_matrix_impl<T>& m, const p_dense_matrix_impl<T>& rhs_diag) {
        assert(false); printf("NOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO <- GEMM DIAG INPLACE");
        //ambient::push< ambient::gemm_diag_inplace<T> >(m, rhs);
    }

    template <typename T>
    inline void scale_inplace(p_dense_matrix_impl<T>& a, const scalar_type& rhs) {
        // gs
#ifdef AMBIENT_SERIAL_CHECK
        alps::numeric::matrix<T> s = maquis::traits::matrix_cast<alps::numeric::matrix<T> >(a);
        s *= (T)rhs;
#endif
        USE_ATOMIC(a.atomic(), scale, a, rhs);
#ifdef AMBIENT_SERIAL_CHECK
        if(s == a){} else printf("--------------------- SCALE WAS INCORRECT!\n");
#endif
    }

    template<typename T>
    inline void copy(p_dense_matrix_impl<T>& ac, const p_dense_matrix_impl<T>& a){
        USE_ATOMIC(a.atomic(), copy, ac, a);
    }
} 
// }}} end of implementation specific type-nested algorithms //

} } // namespace maquis::types

#undef size_type
#undef scalar_type
#undef difference_type 
#endif
