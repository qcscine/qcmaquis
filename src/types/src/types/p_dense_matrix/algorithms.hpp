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

    template<typename T>
    void gemm(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c){
        assert(num_cols(a) == num_rows(b));
        c.resize(a.num_rows(), b.num_cols());
        ambient::push(ambient::gemm_l<T>, ambient::gemm_c<T>, *a.impl, *b.impl, *c.impl);
    }

    template<typename T>
    void gemm(const p_dense_matrix<T>& a, const p_diagonal_matrix<T>& b, p_dense_matrix<T>& c){
        assert(num_cols(a) == num_rows(b));
        c.resize(a.num_rows(), b.num_cols());
        ambient::push(ambient::gemm_diagonal_rhs_l<T>, ambient::gemm_diagonal_rhs_c<T>, *a.impl, *b.get_data().impl, *c.impl);
    }

    template<typename T>
    void gemm(const p_diagonal_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c){
        assert(num_cols(a) == num_rows(b));
        c.resize(a.num_rows(), b.num_cols());
        ambient::push(ambient::gemm_diagonal_lhs_l<T>, ambient::gemm_diagonal_lhs_c<T>, *a.get_data().impl, *b.impl, *c.impl);
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
        ambient::push(ambient::svd_l<T>, ambient::svd_c<T>, *a.impl, m, n, *u.impl, *vt.impl, *s.get_data().impl);
    }

    template<typename T>
    void heev(p_dense_matrix<T> a, p_dense_matrix<T>& evecs,
              typename associated_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    {
        assert(num_rows(a) == num_cols(a));
        assert(num_rows(evals) == num_rows(a));
        int m = num_rows(a);
        evecs.resize(m, m);
        ambient::push(ambient::heev_l, ambient::heev_c, *a.impl, m, *evals.get_data().impl); // destoys U triangle of M
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
        ambient::push(ambient::associated_copy_l<T>, ambient::copy_c<T>, *sc.get_data().impl, *s.get_data().impl);
    }

    template<typename T>
    void copy_sqr_gt(std::vector<typename T::value_type>& sc, typename associated_diagonal_matrix<T>::type& s, const double& prec){
    // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector
        std::vector<typename T::value_type>* sc_ptr = &sc;
        ambient::push(ambient::push_back_sqr_gt_l, ambient::push_back_sqr_gt_c, sc_ptr, *s.get_data().impl, prec);
    }

    template<typename T>
    void copy_after(std::vector<typename T::value_type>& sc, const size_type pos, typename associated_diagonal_matrix<T>::type& s){
    // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector
        std::vector<typename T::value_type>* sc_ptr = &sc;  
        ambient::push(ambient::copy_after_std_l, ambient::copy_after_std_c, sc_ptr, pos, *s.get_data().impl);
    }

    template<typename T>
    void copy_after(typename associated_vector<T>::type& sc, const size_type pos, typename associated_diagonal_matrix<T>::type& s){
    // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector 
        ambient::push(ambient::copy_after_l<T>, ambient::copy_after_c<T>, *sc.get_data().impl, pos, *s.get_data().impl);
    }

    #undef size_type
    #undef difference_type 

} } } // namespace maquis::types::algorithms

#endif
