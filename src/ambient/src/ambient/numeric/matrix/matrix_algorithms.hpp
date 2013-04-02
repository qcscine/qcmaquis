#ifndef AMBIENT_NUMERIC_MATRIX_ALGORITHMS
#define AMBIENT_NUMERIC_MATRIX_ALGORITHMS

#include "ambient/numeric/matrix/matrix.h"
#include "ambient/numeric/kernels/kernels.hpp"
#include "ambient/numeric/matrix/diagonal_matrix.hpp"

#define size_type       typename matrix<T>::size_type
#define real_type       typename matrix<T>::real_type
#define scalar_type     typename matrix<T>::scalar_type
#define difference_type typename matrix<T>::difference_type

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

    template<typename T>
    inline matrix<T> exp(const matrix<T>& a, const T& alfa = 1.){
        printf("ERROR: NOT TESTED (EXP)\n");
        assert( false );
        // TODO: - right eigenvalues/eigenvectors from geev
        //       - inverse of eigenvectors (getrf && getri)
        //       - return Nr * exp(S, alpha) * inverse(Nr);
        return matrix<T>();
    }

    template<typename T>
    inline matrix<T> exp_hermitian(const matrix<T>& a, const T& alfa = 1.){
        printf("ERROR: NOT TESTED (EXP_HERMITIAN)\n");
        diagonal_matrix<T> evals(a.num_rows());
        matrix<T> evecs = matrix<T>();
        heev(a, evecs, evals);
        return (evecs * exp(evals, alfa))*conj(transpose(evecs));
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
    inline void adjoint_inplace(matrix<T>& a){
        transpose_inplace(a);
    }

    template<typename T>
    inline transpose_view<matrix<T> > adjoint(const matrix<T>& a){
        return transpose_view<matrix<T> >(a);
    }

    template<typename T>
    inline void remove_rows(matrix<T>& a, size_type i, difference_type k){
        assert( i+k <= a.num_rows() );
        assert(false); printf("ERROR: NOT TESTED (REMOVE ROWS)\n");
        //kernels::remove_rows<T>::spawn<complexity::N2>(a, i, k);
        resize(a, a.num_rows()-k, a.num_cols());
    }

    template<typename T>
    inline void remove_cols(matrix<T>& a, size_type j, difference_type k){
        assert( j+k <= a.num_cols() );
        assert(false); printf("ERROR: NOT TESTED (REMOVE COLS)\n");
        //kernels::remove_cols<T>::spawn<complexity::N2>(a, j, k);
        resize(a, a.num_rows(), a.num_cols()-k);
    }

    template <typename T>
    inline void mul_inplace(matrix<T>& a, const matrix<T>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM INPLACE)\n");
        //kernels::gemm_inplace<T>::spawn(a, rhs);
    }

    template <typename T>
    inline void mul_inplace(matrix<T>& a, const diagonal_matrix<T>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM DIAG INPLACE)\n");
        //kernels::gemm_diag_inplace<T>::spawn<complexity::N2>(a, rhs);
    }

    template <typename T>
    inline void save(const matrix<T>& a, const size_t& tag) { 
        kernels::save<T>::spawn<complexity::N2>(a, tag); 
    }

    template <typename T>
    inline void load(matrix<T>& a, const size_t& tag) { 
        kernels::load<T>::spawn<complexity::N2>(a, tag); 
    }

    template<class Matrix>
    inline void persist(const Matrix& a){
        ambient::persist(a.core);
    }

    template<class MatrixViewA, class MatrixViewB, typename T>
    inline void gemm(const MatrixViewA& a, const MatrixViewB& b, matrix<T>& c){
        kernels::gemm<MatrixViewA,MatrixViewB,T>::spawn<complexity::N3>(a, b, c); 
    }

    template<class MatrixViewA, typename T, typename D>
    inline void gemm(const MatrixViewA& a, const diagonal_matrix<D>& b, matrix<T>& c){ 
        kernels::gemm_diagonal_rhs<MatrixViewA,T,D>::spawn<complexity::N3>(a, b, c); 
    }

    template<class MatrixViewB, typename T, typename D>
    inline void gemm(const diagonal_matrix<D>& a, const MatrixViewB& b, matrix<T>& c){ 
        kernels::gemm_diagonal_lhs<MatrixViewB,T,D>::spawn<complexity::N3>(a, b, c); 
    }

    template<typename T>
    inline void scale(matrix<T>& a, size_t ai, size_t aj, const diagonal_matrix<T>& alfa, size_t alfai, size_t alfaj){
        kernels::scale_offset<T>::spawn<complexity::N2>(a, ai, aj, alfa, alfai);
    }

    template<typename T>
    inline void svd(matrix<T>& a, matrix<T>& u, matrix<T>& vt, diagonal_matrix<double>& s){
        int m = num_rows(a);
        int n = num_cols(a);
        int k = std::min(m,n);
        u.resize(m, k);
        vt.resize(k, n);
        s.resize(k, k);
        kernels::svd<T>::spawn<complexity::N3>(a, u, vt, s);
    }

    template<typename T>
    inline void heev(matrix<T> a, matrix<T>& evecs, diagonal_matrix<double>& evals){
        assert(num_rows(a) == num_cols(a) && num_rows(evals) == num_rows(a));
        kernels::heev<T>::spawn<complexity::N3>(a, evals); // destoys U triangle of M
        evecs.swap(a);
    }

    template<typename T>
    inline void syev(const matrix<T>& a, matrix<T>& evecs, diagonal_matrix<double>& evals){
        heev(a, evecs, evals); // should it be syev instead?
    }

    template<typename T>
    inline void geqrt(matrix<T>& a, matrix<T>& t){
        kernels::geqrt<T>::spawn<complexity::N3>(a, t);
    }

    template<PLASMA_enum TR, typename T>
    inline void ormqr(size_t k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c){
        kernels::ormqr<T,TR>::spawn<complexity::N3>(k, a, t, c);
    }

    template<typename T>
    inline void tsqrt(matrix<T>& a1, matrix<T>& a2, matrix<T>& t){
        kernels::tsqrt<T>::spawn<complexity::N3>(a1, a2, t);
    }

    template<PLASMA_enum TR, typename T>
    inline void tsmqr(size_t k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t){
        kernels::tsmqr<T,TR>::spawn<complexity::N3>(k, a1, a2, v, t);
    }

    template<typename T>
    inline void gelqt(matrix<T>& a, matrix<T>& t){
        kernels::gelqt<T>::spawn<complexity::N3>(a, t);
    }

    template<PLASMA_enum TR, typename T>
    inline void ormlq(size_t k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c){
        kernels::ormlq<T,TR>::spawn<complexity::N3>(k, a, t, c);
    }

    template<typename T>
    inline void tslqt(matrix<T>& a1, matrix<T>& a2, matrix<T>& t){
        kernels::tslqt<T>::spawn<complexity::N3>(a1, a2, t);
    }

    template<PLASMA_enum TR, typename T>
    inline void tsmlq(size_t k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t){
        kernels::tsmlq<T,TR>::spawn<complexity::N3>(k, a1, a2, v, t);
    }

    template<typename T> inline void resize(matrix<T>& a, size_t m, size_t n){ 
        assert(a.num_rows() != 0 && a.num_cols() != 0);
        if(a.num_rows() == m && a.num_cols() == n) return;
        matrix<T> resized(m, n);
        if(!a.core->weak())
            kernels::resize<T>::spawn<complexity::N2>(resized, a, std::min(m, a.num_rows()), std::min(n, a.num_cols()));
        a.swap(resized);
    }

    template<typename T> inline scalar_type trace(const matrix<T>& a){ 
        scalar_type trace(0.);
        kernels::trace<T>::spawn<complexity::N>(a, trace);
        return trace;
    }

    template <typename T>
    inline real_type norm_square(const matrix<T>& a){ 
        real_type norm(0.); 
        kernels::scalar_norm<T>::spawn<complexity::N2>(a, norm); 
        return norm; 
    }

    template <typename T>
    inline scalar_type overlap(const matrix<T>& a, const matrix<T>& b){ 
        scalar_type overlap(0.); 
        kernels::overlap<T>::spawn<complexity::N2>(a, b, overlap); 
        return overlap; 
    }
        
    template<typename T>
    inline void swap(matrix<T>& x, matrix<T>& y){ 
        x.swap(y);                     
    }

    template<typename T>
    inline void transpose_inplace(matrix<T>& a){
        matrix<T> t(a.num_cols(), a.num_rows());
        kernels::transpose_out<T>::spawn<complexity::N2>(a, t);
        a.swap(t);
    }

    template<typename T>
    inline transpose_view<matrix<T> > transpose(const matrix<T>& a){
        return transpose_view<matrix<T> >(a);
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
        kernels::init_value<typename Matrix::value_type>::spawn<complexity::N2>(a, value);
    }

    template <typename T>
    inline void add_inplace(matrix<T>& lhs, const matrix<T>& rhs){ 
        kernels::add<T>::spawn<complexity::N2>(lhs, rhs); 
    }

    template <typename T>
    inline void sub_inplace(matrix<T>& lhs, const matrix<T>& rhs){ 
        kernels::sub<T>::spawn<complexity::N2>(lhs, rhs); 
    }

    template <typename T>
    inline void mul_inplace(matrix<T>& a, const scalar_type& rhs) { 
        kernels::scale<T>::spawn<complexity::N2>(a, rhs); 
    }

    template <typename T>
    inline void div_inplace(matrix<T>& a, const scalar_type& rhs) { 
        kernels::scale_inverse<T>::spawn<complexity::N2>(a, rhs); 
    }

    template<typename T>
    inline void copy(const matrix<T>& src, matrix<T>& dst){
        ambient::fuse(src.core, dst.core);
    }

    template<typename T>
    inline void copy_block(const matrix<T>& src, size_t si, size_t sj, 
                           matrix<T>& dst, size_t di, size_t dj, 
                           size_t m, size_t n)
    {
        kernels::copy_block<T>::spawn<complexity::N2>(src, si, sj, dst, di, dj, m, n); 
    }

    template<typename T>
    inline void copy_rt(const matrix<T>& src, matrix<T>& dst){ 
        kernels::copy_rt<T>::spawn<complexity::N2>(src, dst);
    }

    template<typename T>
    inline void copy_lt(const matrix<T>& src, matrix<T>& dst){ 
        kernels::copy_lt<T>::spawn<complexity::N2>(src, dst);
    }

    template<typename T>
    inline void copy_block_s(const matrix<T>& src, size_t si, size_t sj, 
                             matrix<T>& dst, size_t di, size_t dj, 
                             const matrix<T>& alfa, size_t ai, size_t aj,
                             size_t m, size_t n)
    { 
        kernels::copy_block_s<T>::spawn<complexity::N2>(src, si, sj, dst, di, dj, alfa, ai, aj, m, n);
    }

    template<typename T>
    inline void copy_block_sa(const matrix<T>& src, size_t si, size_t sj, 
                              matrix<T>& dst, size_t di, size_t dj, 
                              const matrix<T>& alfa, size_t ai, size_t aj,
                              size_t m, size_t n)
    { 
        kernels::copy_block_sa<T>::spawn<complexity::N2>(src, si, sj, dst, di, dj, alfa, ai, aj, m, n);
    }

    template<typename T>
    bool operator == (const matrix<T>& a, const matrix<T>& b){
        if(num_cols(a) != num_cols(b) || num_rows(a) != num_rows(b)){
            printf("Sizes are different: %lu x %lu against %lu x %lu\n",
                    num_cols(a), num_rows(a), num_cols(b), num_rows(b));
            return false;
        }
        ambient::numeric::future<bool> ret(true);
        kernels::validation<T>::spawn<complexity::N2>(a, b, ret);
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

    template <typename T>
    std::ostream& operator << (std::ostream& o, const matrix<T>& a){
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
        return a.core->spec.dim.x;
    }

    template<typename T> inline size_type num_cols(const transpose_view< matrix<T> >& a){
        return a.core->spec.dim.y;
    }

} }

#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type
#endif
