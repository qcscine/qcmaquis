#ifndef AMBIENT_NUMERIC_MATRIX_ALGORITHMS
#define AMBIENT_NUMERIC_MATRIX_ALGORITHMS

#include "ambient/numeric/matrix/matrix.h"
#include "ambient/numeric/kernels/kernels.hpp"
#include "ambient/numeric/matrix/diagonal_matrix.hpp"

#define size_type       typename matrix<T,A>::size_type
#define real_type       typename matrix<T,A>::real_type
#define scalar_type     typename matrix<T,A>::scalar_type
#define difference_type typename matrix<T,A>::difference_type

namespace ambient { namespace numeric {

    template<typename T, class A>
    bool is_hermitian(const matrix<T,A>& a)
    {
        if(num_rows(a) != num_cols(a))
            return false;
        for (size_t i=0; i < num_rows(a); ++i)
            for(size_t j=0; j < num_cols(a); ++j){
               if ( a(i,j) != ambient::numeric::kernels::helper_complex<T>::conj(a(j,i)))
                   return false;
            }
        return true;
    }

    template<typename T, class A>
    inline matrix<T,A> exp(const matrix<T,A>& a, const T& alfa = 1.){
        printf("ERROR: NOT TESTED (EXP)\n");
        assert( false );
        // TODO: - right eigenvalues/eigenvectors from geev
        //       - inverse of eigenvectors (getrf && getri)
        //       - return Nr * exp(S, alpha) * inverse(Nr);
        return matrix<T,A>();
    }

    template<typename T, class A>
    inline matrix<T,A> exp_hermitian(const matrix<T,A>& a, const T& alfa = 1.){
        printf("ERROR: NOT TESTED (EXP_HERMITIAN)\n");
        diagonal_matrix<T> evals(a.num_rows());
        matrix<T,A> evecs = matrix<T,A>();
        heev(a, evecs, evals);
        return (evecs * exp(evals, alfa))*conj(transpose(evecs));
    }

    template<typename T, class A>
    inline matrix<T,A> conj(const matrix<T,A>& a){
        matrix<T,A> b(a);
        conj_inplace(b);
        return b;
    }

    template<class A>
    inline void conj_inplace(matrix<std::complex<double>, A>& a){
        kernels::conj_inplace<std::complex<double>,A>::spawn<complexity::N>(a);
    }

    template<class A>
    inline void conj_inplace(matrix<double, A>& a){ }

    template<typename T, class A>
    inline void remove_rows(matrix<T,A>& a, size_type i, difference_type k){
        assert( i+k <= a.num_rows() );
        assert(false); printf("ERROR: NOT TESTED (REMOVE ROWS)\n");
        //kernels::remove_rows<T>::spawn<complexity::N2>(a, i, k);
        resize(a, a.num_rows()-k, a.num_cols());
    }

    template<typename T, class A>
    inline void remove_cols(matrix<T,A>& a, size_type j, difference_type k){
        assert( j+k <= a.num_cols() );
        assert(false); printf("ERROR: NOT TESTED (REMOVE COLS)\n");
        //kernels::remove_cols<T>::spawn<complexity::N2>(a, j, k);
        resize(a, a.num_rows(), a.num_cols()-k);
    }

    template <typename T, class A>
    inline void mul_inplace(matrix<T,A>& a, const matrix<T,A>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM INPLACE)\n");
        //kernels::gemm_inplace<T>::spawn(a, rhs);
    }

    template <typename T, class A>
    inline void mul_inplace(matrix<T,A>& a, const diagonal_matrix<T>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM DIAG INPLACE)\n");
        //kernels::gemm_diag_inplace<T>::spawn<complexity::N2>(a, rhs);
    }

    template <typename T, class A>
    inline void save(const matrix<T,A>& a, const size_t& tag) { 
        kernels::save<T>::spawn<complexity::N2>(a, tag); 
    }

    template <typename T, class A>
    inline void load(matrix<T,A>& a, const size_t& tag) { 
        kernels::load<T>::spawn<complexity::N2>(a, tag); 
    }

    template<class Matrix>
    inline void persist(const Matrix& a){
        ambient::persist(a.core);
    }

    template<class MatrixViewA, class MatrixViewB, typename T, class A>
    inline void gemm(const MatrixViewA& a, const MatrixViewB& b, matrix<T,A>& c){
        kernels::gemm<MatrixViewA,MatrixViewB,matrix<T,A>,T>::spawn<complexity::N3>(a, b, c); 
    }

    template<class MatrixViewA, typename T, typename D, class A>
    inline void gemm(const MatrixViewA& a, const diagonal_matrix<D>& b, matrix<T,A>& c){ 
        kernels::gemm_diagonal_rhs<MatrixViewA,T,D>::spawn<complexity::N3>(a, b, c); 
    }

    template<class MatrixViewB, typename T, typename D, class A>
    inline void gemm(const diagonal_matrix<D>& a, const MatrixViewB& b, matrix<T,A>& c){ 
        kernels::gemm_diagonal_lhs<MatrixViewB,T,D>::spawn<complexity::N3>(a, b, c); 
    }

    template<typename T, class A>
    inline void scale(matrix<T,A>& a, size_t ai, size_t aj, const diagonal_matrix<T>& alfa, size_t alfai, size_t alfaj){
        kernels::scale_offset<T>::spawn<complexity::N2>(a, ai, aj, alfa, alfai);
    }

    template<typename T, class A>
    inline void svd(matrix<T,A>& a, matrix<T,A>& u, matrix<T,A>& vt, diagonal_matrix<double>& s){
        int m = num_rows(a);
        int n = num_cols(a);
        int k = std::min(m,n);
        u.resize(m, k);
        vt.resize(k, n);
        s.resize(k, k);
        kernels::svd<T>::spawn<complexity::N3>(a, u, vt, s);
    }

    template<typename T, class A>
    inline void inverse(matrix<T,A>& a){
        kernels::inverse<T>::spawn<complexity::N3>(a);
    }

    template<typename T, class A>
    inline void geev(matrix<T,A> a, matrix<T,A>& lv, matrix<T,A>& rv, diagonal_matrix<T>& s){
        kernels::geev<T>::spawn<complexity::N3>(a, lv, rv, s); 
    }

    template<typename T, class A>
    inline void heev(matrix<T,A> a, matrix<T,A>& evecs, diagonal_matrix<double>& evals){
        assert(num_rows(a) == num_cols(a) && num_rows(evals) == num_rows(a));
        kernels::heev<T>::spawn<complexity::N3>(a, evals); // destoys U triangle of M
        evecs.swap(a);
    }

    template<typename T, class A>
    inline void syev(const matrix<T,A>& a, matrix<T,A>& evecs, diagonal_matrix<double>& evals){
        heev(a, evecs, evals); // should it be syev instead?
    }

    template<typename T, class A>
    inline void geqrt(matrix<T,A>& a, matrix<T,A>& t){
        kernels::geqrt<T>::spawn<complexity::N3>(a, t);
    }

    template<PLASMA_enum TR, typename T, class A>
    inline void ormqr(size_t k, const matrix<T,A>& a, const matrix<T,A>& t, matrix<T,A>& c){
        kernels::ormqr<T,TR>::spawn<complexity::N3>(k, a, t, c);
    }

    template<typename T, class A>
    inline void tsqrt(matrix<T,A>& a1, matrix<T,A>& a2, matrix<T,A>& t){
        kernels::tsqrt<T>::spawn<complexity::N3>(a1, a2, t);
    }

    template<PLASMA_enum TR, typename T, class A>
    inline void tsmqr(size_t k, matrix<T,A>& a1, matrix<T,A>& a2, const matrix<T,A>& v, const matrix<T,A>& t){
        kernels::tsmqr<T,TR>::spawn<complexity::N3>(k, a1, a2, v, t);
    }

    template<typename T, class A>
    inline void gelqt(matrix<T,A>& a, matrix<T,A>& t){
        kernels::gelqt<T>::spawn<complexity::N3>(a, t);
    }

    template<PLASMA_enum TR, typename T, class A>
    inline void ormlq(size_t k, const matrix<T,A>& a, const matrix<T,A>& t, matrix<T,A>& c){
        kernels::ormlq<T,TR>::spawn<complexity::N3>(k, a, t, c);
    }

    template<typename T, class A>
    inline void tslqt(matrix<T,A>& a1, matrix<T,A>& a2, matrix<T,A>& t){
        kernels::tslqt<T>::spawn<complexity::N3>(a1, a2, t);
    }

    template<PLASMA_enum TR, typename T, class A>
    inline void tsmlq(size_t k, matrix<T,A>& a1, matrix<T,A>& a2, const matrix<T,A>& v, const matrix<T,A>& t){
        kernels::tsmlq<T,TR>::spawn<complexity::N3>(k, a1, a2, v, t);
    }

    template<typename T, class A> 
    inline void resize(matrix<T,A>& a, size_t m, size_t n){ 
        assert(a.num_rows() != 0 && a.num_cols() != 0);
        if(a.num_rows() == m && a.num_cols() == n) return;
        matrix<T,A> resized(m, n);
        if(!a.core->weak())
            kernels::resize<T,A>::spawn<complexity::N2>(resized, a, std::min(m, a.num_rows()), std::min(n, a.num_cols()));
        a.swap(resized);
    }

    template<typename T, class A> 
    inline scalar_type trace(const matrix<T,A>& a){ 
        scalar_type trace(0.);
        kernels::trace<T,A>::spawn<complexity::N>(a, trace);
        return trace;
    }

    template <typename T, class A>
    inline real_type norm_square(const matrix<T,A>& a){ 
        real_type norm(0.); 
        kernels::scalar_norm<T>::spawn<complexity::N2>(a, norm); 
        return norm; 
    }

    template <typename T, class A>
    inline scalar_type overlap(const matrix<T,A>& a, const matrix<T,A>& b){ 
        scalar_type overlap(0.); 
        kernels::overlap<T>::spawn<complexity::N2>(a, b, overlap); 
        return overlap; 
    }
        
    template<typename T, class A>
    inline void swap(matrix<T,A>& x, matrix<T,A>& y){ 
        x.swap(y);                     
    }

    template<typename T, class A>
    inline void transpose_inplace(matrix<T,A>& a){
        matrix<T,A> t(a.num_cols(), a.num_rows());
        kernels::transpose_out<T,A>::spawn<complexity::N2>(a, t);
        a.swap(t);
    }

    template<typename T, class A>
    inline transpose_view<matrix<T,A> > transpose(const matrix<T,A>& a){
        return transpose_view<matrix<T,A> >(a);
    }

    template<class Matrix>
    inline void fill_identity(Matrix& a){
        kernels::init_identity<typename Matrix::value_type>::spawn<complexity::N2>(a);
    }

    template<class Matrix>
    inline void fill_random(Matrix& a){
        kernels::init_random<typename Matrix::value_type>::spawn<complexity::N2>(a);
        // ambient::sync(); // uncomment for reproduceability
    }

    template<class Matrix>
    inline void fill_value(Matrix& a, typename Matrix::value_type value){
        if(value == 0.) return; // matrices are 0 by default
        kernels::init_value<typename Matrix::value_type, typename Matrix::allocator_type>::spawn<complexity::N2>(a, value);
    }

    template <typename T, class A>
    inline void add_inplace(matrix<T,A>& lhs, const matrix<T,A>& rhs){ 
        kernels::add<T,A>::spawn<complexity::N2>(lhs, rhs); 
    }

    template <typename T, class A>
    inline void sub_inplace(matrix<T,A>& lhs, const matrix<T,A>& rhs){ 
        kernels::sub<T>::spawn<complexity::N2>(lhs, rhs); 
    }

    template <typename T, class A>
    inline void mul_inplace(matrix<T,A>& a, const scalar_type& rhs) { 
        kernels::scale<T>::spawn<complexity::N2>(a, rhs); 
    }

    template <typename T, class A>
    inline void div_inplace(matrix<T,A>& a, const scalar_type& rhs) { 
        kernels::scale_inverse<T>::spawn<complexity::N2>(a, rhs); 
    }

    template<typename T, class A>
    inline void copy(const matrix<T,A>& src, matrix<T,A>& dst){
        ambient::fuse(src.core, dst.core);
    }

    template<typename T, class A>
    inline void copy_block(const matrix<T,A>& src, size_t si, size_t sj, 
                           matrix<T,A>& dst, size_t di, size_t dj, 
                           size_t m, size_t n)
    {
        kernels::copy_block<T>::spawn<complexity::N2>(src, si, sj, dst, di, dj, m, n); 
    }

    template<typename T, class A>
    inline void copy_rt(const matrix<T,A>& src, matrix<T,A>& dst){ 
        kernels::copy_rt<T>::spawn<complexity::N2>(src, dst);
    }

    template<typename T, class A>
    inline void copy_lt(const matrix<T,A>& src, matrix<T,A>& dst){ 
        kernels::copy_lt<T>::spawn<complexity::N2>(src, dst);
    }

    template<typename T, class A>
    inline void copy_block_s(const matrix<T,A>& src, size_t si, size_t sj, 
                             matrix<T,A>& dst, size_t di, size_t dj, 
                             const matrix<T,A>& alfa, size_t ai, size_t aj,
                             size_t m, size_t n)
    { 
        kernels::copy_block_s<T>::spawn<complexity::N2>(src, si, sj, dst, di, dj, alfa, ai, aj, m, n);
    }

    template<class A1, class A2, class A3, typename T>
    inline void copy_block_sa(const matrix<T,A1>& src, size_t si, size_t sj, 
                              matrix<T,A2>& dst, size_t di, size_t dj, 
                              const matrix<T,A3>& alfa, size_t ai, size_t aj,
                              size_t m, size_t n)
    { 
        kernels::copy_block_sa<A1,A2,A3,T>::spawn<complexity::N2>(src, si, sj, dst, di, dj, alfa, ai, aj, m, n);
    }

    template<typename T, class A>
    bool operator == (const matrix<T,A>& a, const matrix<T,A>& b){
        if(num_cols(a) != num_cols(b) || num_rows(a) != num_rows(b)){
            printf("Sizes are different: %lu x %lu against %lu x %lu\n",
                    num_cols(a), num_rows(a), num_cols(b), num_rows(b));
            return false;
        }
        ambient::numeric::future<bool> ret(true);
        kernels::validation<T>::spawn<complexity::N2>(a, b, ret);
        return (bool)ret;
    }

    template<typename T, class A>
    bool operator == (matrix<T,A> a, const transpose_view<matrix<T,A> >& b){
        if(num_cols(a) != num_cols(b) || num_rows(a) != num_rows(b)){
            printf("Sizes are different: %lu x %lu against %lu x %lu\n",
                    num_cols(a), num_rows(a), num_cols(b), num_rows(b));
            return false;
        }
        transpose_inplace(a);
        ambient::numeric::future<bool> ret(true);
        kernels::validation<T>::spawn<complexity::N2>(a, b, ret);
        return (bool)ret;
    }

    template<typename T>
    bool operator == (const diagonal_matrix<T>& a, const diagonal_matrix<T>& b){
        if(num_rows(a) != num_rows(b)){
            printf("Sizes are different: %lu x %lu against %lu x %lu\n",
                    num_cols(a), num_rows(a), num_cols(b), num_rows(b));
            return false;
        }
        ambient::numeric::future<bool> ret(true);
        kernels::validation<T>::spawn<complexity::N2>(a, b, ret);
        return (bool)ret;
    }

    template <typename T, class A>
    std::ostream& operator << (std::ostream& o, const matrix<T,A>& a){
        std::cout.precision(2);
        std::cout.setf(std::ios::fixed, std::ios::floatfield);
        for(size_type i=0; i< a.num_rows(); ++i){
            for(size_type j=0; j < a.num_cols(); ++j){
                std::cout << a(i,j) << " ";
            }
            ambient::cout << "\n\n";
        }
        return o;
    }

    template <typename T, class A> 
    inline matrix<T,A> operator + (matrix<T,A> lhs, const matrix<T,A>& rhs){
        return (lhs += rhs);
    }

    template <typename T, class A> 
    inline matrix<T,A> operator - (matrix<T,A> lhs, const matrix<T,A>& rhs){
        return (lhs -= rhs);
    }

    template <typename T, class A>
    inline const matrix<T,A> operator * (matrix<T,A> lhs, const matrix<T,A>& rhs){
        return (lhs *= rhs);
    }

    template<typename T, typename T2, class A> 
    inline const matrix<T,A> operator * (matrix<T,A> lhs, const T2& rhs){
        return (lhs *= rhs);
    }

    template<typename T, typename T2, class A> 
    inline const matrix<T,A> operator * (const T2& lhs, matrix<T,A> rhs){
        return (rhs *= lhs);
    }

    template<typename T, class A> inline size_type num_rows(const matrix<T,A>& a){
        return a.num_rows();
    }

    template<typename T, class A> inline size_type num_cols(const matrix<T,A>& a){
        return a.num_cols();
    }

    template<typename T, class A> inline size_type num_rows(const transpose_view< matrix<T,A> >& a){
        return a.core->spec.dim.x;
    }

    template<typename T, class A> inline size_type num_cols(const transpose_view< matrix<T,A> >& a){
        return a.core->spec.dim.y;
    }

} }

#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type
#endif
