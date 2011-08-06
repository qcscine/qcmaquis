#ifndef __ALPS_P_DENSE_MATRIX_ALGORITHMS_HPP__
#define __ALPS_P_DENSE_MATRIX_ALGORITHMS_HPP__

#include "p_dense_matrix/concept/matrix_concept_check.hpp"

#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>
#include <boost/numeric/bindings/std/vector.hpp>

#include "p_dense_matrix/p_dense_matrix.h"
#include "utils/timings.h"

namespace blas
{
    template<typename T>
    p_dense_matrix<T> transpose(const p_dense_matrix<T>& m)
    {   // TODO: perhaps this could return a proxy object
        //printf("transpose: %d %d\n", m.num_rows(), m.num_cols());
        p_dense_matrix<T> mt(num_cols(m), num_rows(m));
        ambient::push(ambient::transpose_l, ambient::transpose_c, mt, m);
        return mt;
    }

    template<typename T>
    const T trace(const p_dense_matrix<T>& m)
    {
        assert(num_rows(m) == num_cols(m));
        T* tr = (T*)calloc(1,sizeof(T));
        if(m.is_abstract()) m.touch();
        ambient::push(ambient::trace_l, ambient::trace_c, m, tr);
        ambient::playout(); // to remove in future
        return *tr;
    }
        
    template<typename T>
    p_dense_matrix<T> conjugate(p_dense_matrix<T> m)
    {
        m.inplace_conjugate();
        return m;
    }
        
    template<typename T>
    void pblas_gemm(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c)
    {
        //c.resize(a.num_rows(), b.num_cols());
        c.set_init(ambient::null_i<T>);
        ambient::push(ambient::gemm_l_scalapack, ambient::gemm_c_scalapack, a, b, c);
    }

    template<typename T>
    void gemm(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c)
    {
        //c.resize(a.num_rows(), b.num_cols());
        c.set_init(ambient::null_i<T>);
        //printf("gemm: %d %d\n", c.num_rows(), c.num_cols());
        ambient::push(ambient::gemm_l, ambient::gemm_c, a, b, c);
    }

    template<typename T>
    void gemm(const p_dense_matrix<T>& a, const p_diagonal_matrix<T>& b, p_dense_matrix<T>& c)
    {
        assert(num_cols(a) == num_rows(b));
        c.resize(a.num_rows(), b.num_cols());
        //printf("gemm: %d %d\n", c.num_rows(), c.num_cols());
        ambient::push(ambient::gemm_diagonal_rhs_l, ambient::gemm_diagonal_rhs_c, a, b.get_data(), c);
    }
    
    template<typename T>
    void gemm(const p_diagonal_matrix<T>& a, const p_dense_matrix<T>& b, p_dense_matrix<T>& c)
    {
        assert(num_cols(a) == num_rows(b));
        c.resize(a.num_rows(), b.num_cols());
        //printf("gemm: %d %d\n", c.num_rows(), c.num_cols());
        ambient::push(ambient::gemm_diagonal_lhs_l,ambient::gemm_diagonal_lhs_c, a.get_data(), b, c);
    }
    
    template<typename T>
    void svd(const p_dense_matrix<T>& a,
                   p_dense_matrix<T>& u,
                   p_dense_matrix<T>& vt,
             typename associated_diagonal_matrix<p_dense_matrix<T> >::type& s)
    {
        int m = num_rows(a);
        int n = num_cols(a);
        int k = std::min(m,n);
        u.resize(m, k);
        vt.resize(k, n);
        s.resize(k, k);
        ambient::push(ambient::svd_l_scalapack, ambient::svd_c_scalapack, a, m, n, u, vt, s.get_data());
    }

    template<typename T>
    void heev(p_dense_matrix<T> a,
              p_dense_matrix<T>& evecs,
              typename associated_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    {
        assert(num_rows(a) == num_cols(a));
        assert(num_rows(evals) == num_rows(a));
        int m = num_rows(a);

        evecs.resize(m, m);
        ambient::push(ambient::heev_l_scalapack, ambient::heev_c_scalapack, a, m, evals.get_data()); // destoys U triangle of M
        evecs = a;
    }
    
    template<typename T>
    void syev(p_dense_matrix<T> a,
              p_dense_matrix<T>& evecs,
              typename associated_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    {
        heev(a, evecs, evals);
    }
    
    template<typename T>
    void qr(p_dense_matrix<T> m,
            p_dense_matrix<T> & q,
            p_dense_matrix<T> & r)
    {
        assert(false);
        /* implement thin QR decomposition, i.e. for a (m,n) matrix, where m >= n, the result should be
         Q: (m,n)
         R: (n,n) */
    }
 
    template<typename T>
    void validation(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b)
    {
        ambient::push(ambient::validation_l, ambient::validation_c, a, b);
    }
  
    template<typename T>
    void associated_validation(const p_dense_matrix<T>& a, const p_dense_matrix<T>& b )
    {
        ambient::push(ambient::associated_validation_l, ambient::associated_validation_c, a, b);
    }
 
    template<typename T, class G>
    void generate(p_dense_matrix<T>& a, G g)
    {
        assert(a.is_abstract());
        a.set_init(ambient::random_i<T>);
        a.decouple_copy();
        //a.touch(); // to purge (redunant)
    }
 
    template<typename T>
    void copy(typename associated_vector<T>::type& sc, typename associated_diagonal_matrix<T>::type& s)
    { // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector 
        ambient::push(ambient::associated_copy_l, ambient::copy_c, sc.get_data(), s.get_data());
    }

    template<typename T>
    void copy_sqr_gt(std::vector<typename T::value_type>& sc, typename associated_diagonal_matrix<T>::type& s, const double& prec)
    { // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector
        std::vector<typename T::value_type>* sc_ptr = &sc;
        ambient::push(ambient::push_back_sqr_gt_l, ambient::push_back_sqr_gt_c, sc_ptr, s.get_data(), prec);
        ambient::playout(); // the sc reference points to the deleted object if not to issue playout here
    }

    template<typename T>
    void copy_after(std::vector<typename T::value_type>& sc, const size_t pos, typename associated_diagonal_matrix<T>::type& s)
    { // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector
        std::vector<typename T::value_type>* sc_ptr = &sc;  
        ambient::push(ambient::copy_after_std_l, ambient::copy_after_std_c, sc_ptr, pos, s.get_data());
    }
 
    template<typename T>
    void copy_after(typename associated_vector<T>::type& sc, const size_t pos, typename associated_diagonal_matrix<T>::type& s)
    { // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector 
        ambient::push(ambient::copy_after_l, ambient::copy_after_c, sc.get_data(), pos, s.get_data());
    }

    template<class T>
    void find_if(typename associated_vector<T>::type& s, const double& value, size_t* out_value)
    {
        ambient::push(ambient::associated_find_if_l, ambient::associated_find_if_c, s.get_data(), value, out_value);
    }
 
    void variable_free(void* a)
    { 
        ambient::push(ambient::variable_free_l, ambient::variable_free_c, a);
    }
 
} /* namespace blas */

#endif
