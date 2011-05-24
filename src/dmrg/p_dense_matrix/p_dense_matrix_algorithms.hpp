#ifndef __ALPS_DENSE_MATRIX_ALGORITHMS_HPP__
#define __ALPS_DENSE_MATRIX_ALGORITHMS_HPP__

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
        assert(false);
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

    template<class T>
    void reverse(typename associated_diagonal_matrix<T>::type & s)
    { // reverse only the first col
        size_t num_rows = s.num_rows();
        ambient::push(ambient::associated_reverse_l, ambient::associated_reverse_c, s.get_data(), num_rows);
    }
    
    template<typename T>
    void svd(const p_dense_matrix<T>& a,
                   p_dense_matrix<T>& u,
                   p_dense_matrix<T>& v,
             typename associated_diagonal_matrix<p_dense_matrix<T> >::type& s)
    {
        BOOST_CONCEPT_ASSERT((blas::Matrix<p_dense_matrix<T> >));
        typename p_dense_matrix<T>::size_type k = std::min(num_rows(a), num_cols(a));
        u.resize(num_rows(a), k);
        v.resize(k, num_cols(a));
        s.resize(k, k);
        //printf("svd: %d %d; %d %d; %d %d\n", a.num_rows(), a.num_cols(), u.num_rows(), u.num_cols(), v.num_rows(), v.num_cols());
        ambient::push(ambient::svd_l_scalapack, ambient::svd_c_scalapack, a, u, v, s.get_data());
    }

    template<typename T>
    void syev(p_dense_matrix<T> m,
              p_dense_matrix<T>& evecs,
              typename associated_diagonal_matrix< p_dense_matrix<T> >::type& evals)
    {
        assert(num_rows(m) == num_cols(m));
        assert(num_rows(evals) == num_rows(m));

        evecs.resize(num_rows(m), num_cols(m));
        ambient::push(ambient::syev_l_scalapack, ambient::syev_c_scalapack, m, evals.get_data(), evecs); // destoys U triangle of M
        reverse< p_dense_matrix<T> >(evals);          
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
        a.set_init(ambient::random_i<T>);
    }
 
    template<typename T>
    void copy(typename associated_vector<T>::type& sc, typename associated_diagonal_matrix<T>::type& s)
    { // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector 
        ambient::push(ambient::associated_copy_l, ambient::copy_c, sc.get_data(), s.get_data());
    }
 
    template<typename T>
    void copy_after(typename associated_vector<T>::type& sc, const size_t pos, typename associated_diagonal_matrix<T>::type& s)
    { // this kernel copies only the first cols of the work group, only used with associated_diagonal_matrix and associated_vector 
        ambient::push(ambient::copy_after_l, ambient::copy_after_c, sc.get_data(), pos, s.get_data());
    }
    
    template<class T>
    void sort(typename associated_vector<T>::type& s)
    {
        ambient::playout(); // because of bugbug place...
        ambient::push(ambient::associated_sort_l, ambient::associated_sort_c, s.get_data());
        size_t grid_dim_y = get_grid_dim(s.get_data()).y; // bugbug
        size_t num = __a_ceil(grid_dim_y/2);
        for(size_t i=0; i < num ; i++){
            ambient::push(ambient::associated_oe_sort_l, ambient::associated_sort_o_c, s.get_data());
            ambient::push(ambient::associated_oe_sort_l, ambient::associated_sort_e_c, s.get_data());
        }
        ambient::push(ambient::associated_oe_sort_l, ambient::move_offset_c, s.get_data());
    }

    template<class T>
    void find_if(typename associated_vector<T>::type& s, const double& value, size_t* out_value)
    {
        ambient::push(ambient::associated_find_if_l, ambient::associated_find_if_c, s.get_data(), value, out_value);
    }

    template<class T>
    void accumulate(typename associated_vector<T>::type& s, const size_t* begin, double* out_value)
    {
        ambient::push(ambient::associated_accumulate_l, ambient::associated_accumulate_c, s.get_data(), begin, out_value);
    }

    template<class T>
    void max(typename associated_vector<T>::type& s, const double& evalscut, const size_t& mmax, double* out_value) 
    {
        ambient::push(ambient::associated_max_l, ambient::associated_max_c, s.get_data(), evalscut, mmax, out_value); 
    }
 
    void variable_free(void* a)
    { 
        ambient::push(ambient::variable_free_l, ambient::variable_free_c, a);
    }
 
} /* namespace blas */

#endif
