#ifndef AMBIENT_NUMERIC_MATRIX_KERNELS
#define AMBIENT_NUMERIC_MATRIX_KERNELS

#include "ambient/numeric/kernels/math.hpp"
#include "ambient/numeric/kernels/utils.hpp"
#include "ambient/numeric/traits.hpp"
#include "ambient/utils/numeric.h"

namespace ambient { namespace numeric { namespace kernels {
    using ambient::unbound;
    using ambient::numeric::matrix;
    using ambient::numeric::traits::real_type;

    template<typename T> 
    struct geqrt : public kernel< geqrt<T> > 
    { static void c(matrix<T>& a, matrix<T>& t); };

    template<typename T, PLASMA_enum TR>
    struct ormqr : public kernel< ormqr<T,TR> > 
    { static void c(const size_t& k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c); };

    template<typename T>
    struct tsqrt : public kernel< tsqrt<T> > 
    { static void c(matrix<T>& a1, matrix<T>& a2, matrix<T>& t); };

    template<typename T, PLASMA_enum TR>
    struct tsmqr : public kernel< tsmqr<T,TR> > 
    { static void c(const size_t& k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t); };

    template<typename T>
    struct gelqt : public kernel< gelqt<T> > 
    { static void c(matrix<T>& a, matrix<T>& t); };

    template<typename T, PLASMA_enum TR>
    struct ormlq : public kernel< ormlq<T,TR> > 
    { static void c(const size_t& k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c); };

    template<typename T>
    struct tslqt : public kernel< tslqt<T> > 
    { static void c(matrix<T>& a1, matrix<T>& a2, matrix<T>& t); };

    template<typename T, PLASMA_enum TR>
    struct tsmlq : public kernel< tsmlq<T,TR> > 
    { static void c(const size_t& k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t); };

    template<class ViewA, class ViewB, class ViewC, typename T>
    struct gemm : public kernel< gemm<ViewA, ViewB, ViewC, T> > 
    { static void c(const matrix<T,typename ViewA::allocator_type>& a, const matrix<T,typename ViewB::allocator_type>& b, unbound< matrix<T,typename ViewC::allocator_type> >& c); };
        
    template<class ViewB, typename T, typename D>
    struct gemm_diagonal_lhs : public kernel< gemm_diagonal_lhs<ViewB,T,D> > 
    { static void c(const matrix<D>& a_diag, const matrix<T>& b, unbound< matrix<T> >& c); };
        
    template<typename T, typename D>
    struct gemm_diagonal_lhs<transpose_view<matrix<T> >,T,D> : public kernel< gemm_diagonal_lhs<transpose_view<matrix<T> >,T,D> > 
    { static void c(const matrix<D>& a_diag, const matrix<T>& b, unbound< matrix<T> >& c); };
        
    template<class ViewA, typename T, typename D>
    struct gemm_diagonal_rhs : public kernel< gemm_diagonal_rhs<ViewA,T,D> > 
    { static void c(const matrix<T>& a, const matrix<D>& b_diag, unbound< matrix<T> >& c); };

    template<typename T, typename D>
    struct gemm_diagonal_rhs<transpose_view<matrix<T> >,T,D> : public kernel< gemm_diagonal_rhs<transpose_view<matrix<T> >,T,D> > 
    { static void c(const matrix<T>& a, const matrix<D>& b_diag, unbound< matrix<T> >& c); };
        
    template<typename T, class A>
    struct trace : public kernel< trace<T,A> > 
    { static void c(const matrix<T,A>& a, future<T>& trace); };
        
    template<typename T>
    struct scalar_norm : public kernel< scalar_norm<T> > 
    { static void c(const matrix<T>& a, future<double>& norm); };
        
    template<typename T>
    struct overlap : public kernel< overlap<T> > 
    { static void c(const matrix<T>& a, const matrix<T>& b, future<T>& overlap); };
        
    template<typename T, class A>
    struct add : public kernel< add<T,A> > 
    { static void c(matrix<T,A>& a, const matrix<T,A>& b); };
        
    template<typename T>
    struct sub : public kernel< sub<T> > 
    { static void c(matrix<T>& a, const matrix<T>& b); };
        
    template<typename T>
    struct scale : public kernel< scale<T> > 
    { static void c(matrix<T>& a, const future<T>& t); };
        
    template<typename T>
    struct scale_offset : public kernel< scale_offset<T> > 
    { static void c(matrix<T>& a, const size_t& ai, const size_t& aj, const matrix<T>& alfa, const size_t& alfai); };
        
    template<typename T>
    struct scale_inverse : public kernel< scale_inverse<T> > 
    { static void c(matrix<T>& a, const future<T>& t); };
        
    template<typename T>
    struct sqrt_diagonal : public kernel< sqrt_diagonal<T> > 
    { static void c(matrix<T>& a); };
        
    template<typename T>
    struct exp_diagonal : public kernel< exp_diagonal<T> > 
    { static void c(matrix<T>& a, const T& alfa); };

    template<typename T, class A>
    struct transpose_out : public kernel< transpose_out<T,A> > 
    { static void c(const matrix<T,A>& a, unbound< matrix<T,A> >& t); };

    template<typename T, class A>
    struct conj_inplace : public kernel< conj_inplace<T,A> > 
    { static void c(matrix<T,A>& a); };

    template<typename T, class A>
    struct resize : public kernel< resize<T,A> > 
    { static void c(unbound< matrix<T,A> >& r, const matrix<T,A>& a, const size_t& m, const size_t& n); };
        
    template<typename T>
    struct init_identity : public kernel< init_identity<T> > 
    { static void c(unbound< matrix<T> >& a); };
        
    template<typename T, class A>
    struct init_value : public kernel< init_value<T,A> > 
    { static void c(unbound< matrix<T,A> >& a, const T& value); };
        
    template<typename T>
    struct round_square : public kernel< round_square<T> > 
    { static void c(const matrix<T>& a, std::vector<T>*& ac); };

    template<typename T>
    struct cast_to_vector : public kernel< cast_to_vector<T> > 
    { static void c(std::vector<T>*& ac, const matrix<T>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset); };
        
    template<typename T>
    struct cast_from_vector : public kernel< cast_from_vector<T> > 
    { static void c(const std::vector<T>*& ac, matrix<T>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset); };

    template<typename T1, typename T2>
    struct cast_from_vector_t : public kernel< cast_from_vector_t<T1,T2> > 
    { static void c(const std::vector<T1>*& ac, matrix<T2>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset); };

    template<typename T,typename D>
    struct cast_double_complex : public kernel< cast_double_complex<T, D> >
    {static void c(matrix<T>& a, const matrix<D>& b); };

    template<typename T, class A>
    struct save : public kernel< save<T,A> >
    { static void c(matrix<T,A> const& a, const size_t& tag); };

    template<typename T, class A>
    struct load : public kernel< load<T,A> >
    { static void c(matrix<T,A>& a, const size_t& tag); };
    
    template<typename T>
    struct svd : public kernel< svd<T> > 
    { static void c(const matrix<T>& a, unbound< matrix<T> >& u, unbound< matrix<T> >& vt, unbound< matrix<typename real_type<T>::type> >& s); };

    template<typename T>
    struct geev : public kernel< geev<T> > 
    { static void c(const matrix<T>& a, unbound< matrix<T> >& lv, unbound< matrix<T> >& rv, unbound< matrix<T> >& s); };

    template<typename T>
    struct inverse : public kernel< inverse<T> > 
    { static void c(matrix<T>& a); };

    template<typename T>
    struct heev : public kernel< heev<T> > 
    { static void c(matrix<T>& a, unbound< matrix<typename real_type<T>::type> >& w); };

    template<typename T>
    struct copy_rt : public kernel< copy_rt<T> > 
    { static void c(const matrix<T>& a, unbound< matrix<T> >& t); };

    template<typename T>
    struct copy_lt : public kernel< copy_lt<T> > 
    { static void c(const matrix<T>& a, unbound< matrix<T> >& t); };

    template<class A1, class A2, typename T>
    struct copy_block : public kernel< copy_block<A1,A2,T> > {
        static void c(const matrix<T,A1>& src, const size_t& si, const size_t& sj,
                      matrix<T,A2>& dst, const size_t& di, const size_t& dj, 
                      const size_t& m, const size_t& n);
    };

    template<typename T>
    struct copy_block_s : public kernel< copy_block_s<T> > {
        static void c(const matrix<T>& src, const size_t& si, const size_t& sj,
                      matrix<T>& dst, const size_t& di, const size_t& dj, 
                      const matrix<T>& alfa, const size_t& ai, const size_t& aj,
                      const size_t& m, const size_t& n);
    };

    template<class A1, class A2, class A3, typename T>
    struct copy_block_sa : public kernel< copy_block_sa<A1,A2,A3,T> > {
        static void c(const matrix<T,A1>& src, const size_t& si, const size_t& sj,
                      matrix<T,A2>& dst, const size_t& di, const size_t& dj, 
                      const matrix<T,A3>& alfa, const size_t& ai, const size_t& aj,
                      const size_t& m, const size_t& n);
    };
       
    template<typename T>
    struct init_random : public kernel< init_random<T> > {
        static void randomize(double& a);
        static void randomize(std::complex<double>& a);
        static void c(unbound< matrix<T> >& a);
    };

    template<typename T>
    struct validation : public kernel< validation<T> > {
        static double distance(const std::complex<double>& a, const std::complex<double>& b);
        static double magnitude(const std::complex<double>& a, const std::complex<double>& b);
        static double distance(double a, double b);
        static double magnitude(double a, double b);
        static void c(const matrix<T>& a, const matrix<T>& b, future<bool>& ret);
    };




////////////////////////////
// Kernels implementation //
////////////////////////////



    template<typename T>
    void geqrt<T>::c(matrix<T>& a, matrix<T>& t){
        T* tau  = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB>(); 
        T* work = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
        helper_plasma<T>::geqrt(a.num_rows(), a.num_cols(), PLASMA_IB,
                                (T*)revised(a), a.num_rows(),
                                (T*)updated(t), t.num_rows(),
                                tau, work);
    }

    template<typename T, PLASMA_enum TR>
    void ormqr<T,TR>::c(const size_t& k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c){
        T* work = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
        helper_plasma<T>::ormqr(PlasmaLeft, TR, c.num_rows(), c.num_cols(), k, PLASMA_IB,
                                (T*)current(a), a.num_rows(),
                                (T*)current(t), t.num_rows(),
                                (T*)revised(c), c.num_rows(),
                                 work, AMBIENT_IB);
    }

    template<typename T>
    void tsqrt<T>::c(matrix<T>& a1, matrix<T>& a2, matrix<T>& t){
        T* tau  = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB>();
        T* work = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
        helper_plasma<T>::tsqrt(a2.num_rows(), a2.num_cols(), PLASMA_IB,
                                (T*)revised(a1), a1.num_rows(),
                                (T*)revised(a2), a2.num_rows(),
                                (T*)updated(t), t.num_rows(),
                                tau, work);
    }

    template<typename T, PLASMA_enum TR>
    void tsmqr<T,TR>::c(const size_t& k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t){
        T* work = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
        helper_plasma<T>::tsmqr(PlasmaLeft, TR,
                                AMBIENT_IB, a1.num_cols(), a2.num_rows(), a2.num_cols(), k, PLASMA_IB,
                                (T*)revised(a1), a1.num_rows(),
                                (T*)revised(a2), a2.num_rows(),
                                (T*)current(v), v.num_rows(),
                                (T*)current(t), t.num_rows(),
                                work, PLASMA_IB);
    }

    template<typename T>
    void gelqt<T>::c(matrix<T>& a, matrix<T>& t){
        T* tau  = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB>();
        T* work = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
        helper_plasma<T>::gelqt(a.num_rows(), a.num_cols(), PLASMA_IB,
                                (T*)revised(a), a.num_rows(), 
                                (T*)updated(t),   t.num_rows(),
                                tau, work);
    }

    template<typename T, PLASMA_enum TR>
    void ormlq<T,TR>::c(const size_t& k, const matrix<T>& a, const matrix<T>& t, matrix<T>& c){
        T* work = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
        helper_plasma<T>::ormlq(PlasmaRight, TR,
                                c.num_rows(), c.num_cols(), k, PLASMA_IB,
                                (T*)current(a), a.num_rows(),
                                (T*)current(t), t.num_rows(),
                                (T*)revised(c), c.num_rows(),
                                work, AMBIENT_IB);
    }

    template<typename T>
    void tslqt<T>::c(matrix<T>& a1, matrix<T>& a2, matrix<T>& t){
        T* tau  = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB>();
        T* work = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
        helper_plasma<T>::tslqt(a2.num_rows(), a2.num_cols(), PLASMA_IB,
                                (T*)revised(a1), a1.num_rows(),
                                (T*)revised(a2), a2.num_rows(),
                                (T*)updated(t),     t.num_rows(),
                                tau, work);
    }

    template<typename T, PLASMA_enum TR>
    void tsmlq<T,TR>::c(const size_t& k, matrix<T>& a1, matrix<T>& a2, const matrix<T>& v, const matrix<T>& t){
        T* work = (T*)ambient::bulk.malloc<sizeof(T)*AMBIENT_IB*PLASMA_IB>();
        helper_plasma<T>::tsmlq(PlasmaRight, TR,
                                a1.num_rows(), AMBIENT_IB, a2.num_rows(), a2.num_cols(), k, PLASMA_IB,
                                (T*)revised(a1), a1.num_rows(),
                                (T*)revised(a2), a2.num_rows(),
                                (T*)current(v), v.num_rows(),
                                (T*)current(t), t.num_rows(),
                                work, AMBIENT_IB);
    }

    template<class ViewA, class ViewB, class ViewC, typename T>
    void gemm<ViewA, ViewB, ViewC, T>::c(const matrix<T,typename ViewA::allocator_type>& a, 
                                         const matrix<T,typename ViewB::allocator_type>& b, 
                                      unbound< matrix<T,typename ViewC::allocator_type> >& c){
        if(!raw(a).valid() || !raw(b).valid()){
            emptied(c);
            return;
        }
        T* ad = current(a);
        T* bd = current(b);
        T* cd = updated(c);
        int m = ViewA::rows(a);
        int k = ViewA::cols(a);
        int n = ViewB::cols(b);
        int lda = a.num_rows();
        int ldb = b.num_rows();
        int ldc = c.num_rows();
        static const double alfa(1.0); 
        static const double beta(0.0);
        helper_blas<T>::gemm(ViewA::code(), ViewB::code(), &m, &n, &k, &alfa, ad, &lda, bd, &ldb, &beta, cd, &ldc);
    }
        
    template<class ViewB, typename T, typename D>
    void gemm_diagonal_lhs<ViewB,T,D>::c(const matrix<D>& a_diag, const matrix<T>& b, unbound< matrix<T> >& c){
        int sizey = a_diag.num_rows();
        int size = b.num_cols();
        T* bd = current(b);
        T* cd = emptied(c);
        D* alfa = current(a_diag);

        for(int k = 0 ; k < sizey; k++)
    	    helper_blas<T>::axpy(&size, &alfa[k], &bd[k], &sizey, &cd[k], &sizey);
    }
        
    template<typename T, typename D>
    void gemm_diagonal_lhs<transpose_view<matrix<T> >,T,D>::c(const matrix<D>& a_diag, const matrix<T>& b, unbound< matrix<T> >& c){
        size_t sizex = b.num_cols();
        int size  = a_diag.num_rows();
        static const int ONE = 1;
        T* bd = current(b);
        T* cd = emptied(c);
        D* alfa = current(a_diag);
    
        for(int k = 0 ; k < sizex; k++)
    	    helper_blas<T>::axpy(&size, &alfa[k], &bd[k*size], &ONE, &cd[k], &size);// C - check carefully for TE a_diag double, b complex
    }
        
    template<class ViewA, typename T, typename D>
    void gemm_diagonal_rhs<ViewA,T,D>::c(const matrix<T>& a, const matrix<D>& b_diag, unbound< matrix<T> >& c){
        size_t sizex = b_diag.num_rows();
        int size = a.num_rows(); // for the case of complex
        static const int ONE = 1;
        T* ad = current(a);
        T* cd = emptied(c);
    	D* alfa = current(b_diag);
    
        for(int k = 0 ; k < sizex; k++)
    	    helper_blas<T>::axpy(&size, &alfa[k], &ad[k*size], &ONE, &cd[k*size], &ONE);
    }

    template<typename T, typename D>
    void gemm_diagonal_rhs<transpose_view<matrix<T> >,T,D>::c(const matrix<T>& a, const matrix<D>& b_diag, unbound< matrix<T> >& c){
        int sizey = b_diag.num_rows();
        int size = a.num_cols();
        static const int ONE = 1;
        T* ad = current(a);
        T* cd = emptied(c);
    	D* alfa = current(b_diag);
    
        for(int k = 0 ; k < sizey; k++)
    	   helper_blas<T>::axpy(&size, &alfa[k], &ad[k], &sizey, &cd[k*size], &ONE);// C - check carefully for TE b_diag double, b complex
    }

    template<typename T>
    void copy_rt<T>::c(const matrix<T>& a, unbound< matrix<T> >& t){
        T* ad  = current(a);
        T* td  = emptied(t);
        size_t sda = a.num_cols();
        size_t lda = a.num_rows();
        size_t ldt = t.num_rows();

        for(int j = 0; j < sda; ++j)
        for(int i = 0; i <= j && i < ldt; ++i)
        td[i+ldt*j] = ad[i+lda*j]; 
    }

    template<typename T>
    void copy_lt<T>::c(const matrix<T>& a, unbound< matrix<T> >& t){
        T* ad  = current(a);
        T* td  = emptied(t);
        size_t sdt = t.num_cols();
        size_t lda = a.num_rows();
        size_t ldt = t.num_rows();

        for(int j = 0; j < sdt; ++j)
        for(int i = j; i < lda; ++i)
        td[i+ldt*j] = ad[i+lda*j]; 
    }

    template<class A1, class A2, typename T>
    void copy_block<A1,A2,T>::c(const matrix<T,A1>& src, const size_t& si, const size_t& sj,
                                matrix<T,A2>& dst, const size_t& di, const size_t& dj, 
                                const size_t& m, const size_t& n)
    {
        T* sd = current(src);
        T* dd = m*n < ambient::square_dim(dst) ? revised(dst) : updated(dst);
        ambient::memptf<T, ambient::memcpy>(dd, dst.num_rows(), dim2(dj, di), 
                                            sd, src.num_rows(), dim2(sj, si), 
                                            dim2( n, m ));
    }

    template<typename T>
    void copy_block_s<T>::c(const matrix<T>& src, const size_t& si, const size_t& sj,
                            matrix<T>& dst, const size_t& di, const size_t& dj, 
                            const matrix<T>& alfa, const size_t& ai, const size_t& aj,
                            const size_t& m, const size_t& n)
    {
        T* sd = current(src);
        T* dd = m*n < ambient::square_dim(dst) ? revised(dst) : updated(dst);
        T factor = ((T*)current(alfa))[ai + aj*alfa.num_rows()];
        ambient::memptf<T, ambient::memscal>(dd, dst.num_rows(), dim2(dj, di), 
                                             sd, src.num_rows(), dim2(sj, si), 
                                             dim2( n, m ), factor);
    }

    template<class A1, class A2, class A3, typename T>
    void copy_block_sa<A1,A2,A3,T>::c(const matrix<T,A1>& src, const size_t& si, const size_t& sj,
                                      matrix<T,A2>& dst, const size_t& di, const size_t& dj, 
                                      const matrix<T,A3>& alfa, const size_t& ai, const size_t& aj,
                                      const size_t& m, const size_t& n)
    {
        T factor = ((T*)current(alfa))[ai + aj*alfa.num_rows()];
        ambient::memptf<T, ambient::memscala>(revised(dst), dst.num_rows(), dim2(dj, di), 
                                              current(src), src.num_rows(), dim2(sj, si), 
                                              dim2( n, m ), factor);
    }
        
    template<typename T, class A>
    void trace<T,A>::c(const matrix<T,A>& a, future<T>& trace){
        size_t m = a.num_rows();
        size_t n = a.num_cols();
        T* ad = current(a);
    
        size_t sizex = std::min(n,m);
        for(size_t jj = 0; jj < sizex; jj++)
            trace.get_naked() += ad[jj + jj*m];
    }
        
    template<typename T>
    void scalar_norm<T>::c(const matrix<T>& a, future<double>& norm){
        T* ad = current(a);
        norm.get_naked() = alps::numeric::real(ambient::dot(ad, ad, ambient::square_dim(a)));
    }
        
    template<typename T>
    void overlap<T>::c(const matrix<T>& a, const matrix<T>& b, future<T>& overlap){
        T* ad = current(a);
        T* bd = current(b);
        overlap.get_naked() = ambient::dot(ad, bd, ambient::square_dim(a));
    }

    template<typename T, class A>
    void add<T,A>::c(matrix<T,A>& a, const matrix<T,A>& b){
        T* ad = current(a);
        T* bd = current(b);
        T* ar = updated(a);

        int size = ambient::square_dim(a);
        for(int k = 0; k < size; k++) 
            ar[k] = ad[k] + bd[k];
    }

        
    template<typename T>
    void sub<T>::c(matrix<T>& a, const matrix<T>& b){
        T* ad = current(a);
        T* bd = current(b);
        T* ar = updated(a);

        int size = ambient::square_dim(a);
        for(int k = 0; k < size; k++) 
            ar[k] = ad[k] - bd[k];
    }
        
    template<typename T>
    void scale<T>::c(matrix<T>& a, const future<T>& t){
        T* ad = current(a);
        T* ar = updated(a);
        T factor = t.get_naked();
        int size = ambient::square_dim(a);
        for(int k=0; k < size; k++) 
            ar[k] = ad[k] * factor;
    }
        
    template<typename T>
    void scale_offset<T>::c(matrix<T>& a, const size_t& ai, const size_t& aj, const matrix<T>& alfa, const size_t& alfai){
        int m = num_rows(a);
        T* ad = &((T*)revised(a))[aj*m];
        T factor = ((T*)current(alfa))[alfai];
        for(int k = ai; k < m; k++) ad[k] *= factor;
    }
        
    template<typename T>
    void scale_inverse<T>::c(matrix<T>& a, const future<T>& t){
        T* ad = current(a);
        T* ar = updated(a);
        T factor = t.get_naked();
        int size = ambient::square_dim(a);
        for(int k=0; k < size; k++) 
            ar[k] = ad[k] / factor;
    }
        
    template<typename T>
    void sqrt_diagonal<T>::c(matrix<T>& a){
        size_t size = a.num_rows();
        T* ad = current(a);
        T* ar = updated(a);
        for(size_t i = 0; i < size; ++i) ar[i] = std::sqrt(ad[i]);
    }
        
    template<typename T>
    void exp_diagonal<T>::c(matrix<T>& a, const T& alfa){
        size_t size = a.num_rows();
        T* ad = current(a);
        T* ar = updated(a);
        for(size_t i = 0; i < size; ++i) ar[i] = std::exp(alfa*ad[i]);
    }

    template<typename T, class A>
    void transpose_out<T,A>::c(const matrix<T,A>& a, unbound< matrix<T,A> >& t){
        T* od = current(a);
        T* td = updated(t);
        int m = a.num_rows();
        int n = a.num_cols();

        for(int i = 0; i < m; i++){
            for(int j = 0; j < n; j++) *td++ = od[j*m];
            od++;
        }
    }

    template<typename T, class A>
    void conj_inplace<T,A>::c(matrix<T,A>& a){
        size_t size = a.num_rows()*a.num_cols();
        T* ad = current(a);
        T* ar = updated(a);
        for(int i=0; i < size; ++i)
            ar[i] = helper_complex<T>::conj(ad[i]);   
    }

    template<typename T, class A>
    void resize<T,A>::c(unbound< matrix<T,A> >& r, const matrix<T,A>& a, const size_t& m, const size_t& n){
        T* dd = m*n == ambient::square_dim(r) ? updated(r) : emptied(r);
        ambient::memptf<T, ambient::memcpy>(dd, r.num_rows(), dim2(0,0),
                                            current(a), a.num_rows(), dim2(0,0), dim2(n, m)); 
    }
        
    template<typename T>
    void init_identity<T>::c(unbound< matrix<T> >& a){
        size_t n = a.num_cols();
        size_t m = a.num_rows();
        T* ad = emptied(a);

        size_t sizex = std::min(m,n); // respecting borders
        for(size_t jj = 0; jj < sizex; ++jj) ad[jj + m*jj] = 1.;
    }
       
    template<typename T> 
    void init_random<T>::randomize(double& a){ 
        a = drand48();
    }

    template<typename T> 
    void init_random<T>::randomize(std::complex<double>& a){
        a.real(drand48());
        a.imag(drand48());
    }

    template<typename T>
    void init_random<T>::c(unbound< matrix<T> >& a){
        size_t size = ambient::square_dim(a);
        T* ad = updated(a);
        for(size_t i = 0; i < size; ++i) randomize(ad[i]);
    }
        
    template<typename T, class A>
    void init_value<T,A>::c(unbound< matrix<T,A> >& a, const T& value){
        size_t size = ambient::square_dim(a);
        T* ad = updated(a);
        for(size_t i = 0; i < size; ++i) ad[i] = value; // not a memset due to complex
    }
        
    template<typename T>
    void round_square<T>::c(const matrix<T>& a, std::vector<T>*& ac){
        T* ad = current(a);
        size_t sizey = a.num_rows();
        for(int i=0; i < sizey; i++){
            double v = std::abs(ad[i]);
            if(v > 1e-10) ac->push_back(v*v);
        }
    }

    template<typename T>
    void cast_to_vector<T>::c(std::vector<T>*& ac, const matrix<T>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset){
        T* ad = current(a);
        for(int j=0; j < n; ++j) std::memcpy((void*)&(*ac)[j*lda + offset],(void*)&ad[j*m], m*sizeof(T));  
    }
        
    template<typename T>
    void cast_from_vector<T>::c(const std::vector<T>*& ac, matrix<T>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset){
        T* ad = updated(a);
        for(int j=0; j < n; ++j) std::memcpy((void*)&ad[j*m],(void*)&(*ac)[offset + j*lda], m*sizeof(T));
    }

    template<typename T1, typename T2>
    void cast_from_vector_t<T1,T2>::c(const std::vector<T1>*& ac, matrix<T2>& a, const size_t& m, const size_t& n, const size_t& lda, const size_t& offset){
        T2* ad = updated(a);
        const T1* sd = &(*ac)[offset];
        for(int j=0; j < n; ++j) 
            for(int i=0; i < m; ++i)
                ad[j*m + i] = sd[j*lda + i];
    }

    template<typename T, typename D>
    void cast_double_complex<T,D>::c(matrix<T>& a, const matrix<D>& b){
        T* ad = updated(a);
        D* bd = current(b);
        size_t size = a.num_rows();
        for(size_t i = 0; i < size; ++i)
            ad[i] = helper_cast<T,D>::cast(bd[i]);
    };

    template<typename T, class A>
    void save<T,A>::c(matrix<T,A> const& a, const size_t& tag ){
        T* ad = (T*)current(a);
        ambient::fout.save(ad, tag, ambient::size(a)); 
    }

    template<typename T, class A>
    void load<T,A>::c(matrix<T,A>& a, const size_t& tag){
        T* ad = (T*)updated(a);
        ambient::fout.load(ad, tag, ambient::size(a)); 
    }

    template<typename T> 
    double validation<T>::distance(const std::complex<double>& a, const std::complex<double>& b){ 
        return fabs(std::norm(a) - std::norm(b));
    }
    template<typename T> 
    double validation<T>::magnitude(const std::complex<double>& a, const std::complex<double>& b){
        return std::max(fabs(std::norm(a)), fabs(std::norm(b)));
    }
    template<typename T> 
    double validation<T>::distance(double a, double b) { 
        return fabs(fabs(a) - fabs(b));    
    }
    template<typename T> 
    double validation<T>::magnitude(double a, double b){ 
        return std::max(fabs(a), fabs(b)); 
    }
    template<typename T>
    void validation<T>::c(const matrix<T>& a, const matrix<T>& b, future<bool>& ret){ // see paper for Reference Dongara 
        T* ad = current(a); 
        T* bd = current(b); 
        double epsilon = std::numeric_limits<double>::epsilon();
        int count = 0;
        size_t sizey = std::min(a.num_rows(), b.num_rows());
        size_t sizex = std::min(a.num_cols(), b.num_cols());
        
        std::cout.precision(16);
        std::cout.setf( std::ios::fixed, std:: ios::floatfield );

        for(size_t i=0; i < sizey; ++i){
            for(size_t j=0; j < sizex; ++j){
                T av = ad[i+j*a.num_rows()];
                T bv = bd[i+j*b.num_rows()];
                double d = distance(av, bv);
                double m = magnitude(av, bv);
                if(d > epsilon*256 && d/m > epsilon*256){ // || av*bv < 0 // 16 is recommended, 256 because MKL isn't bitwise stable
                    std::cout << i << " " << j << " : " << av << " " << bv << ", eps: " << d << "\n";
                    ret.get_naked() = false;
                    if(++count > 10) return;
                }

            }
        }
    }

    template<typename T>
    void svd<T>::c(const matrix<T>& a, unbound< matrix<T> >& u, unbound< matrix<T> >& vt, unbound< matrix<typename real_type<T>::type> >& s){
        int m = a.num_rows();
        int n = a.num_cols();
        int k = std::min(m,n);
        int info;
        int lwork = -1;
        T wkopt;
        T* ad  = current(a);
        T* ud  = updated(u);
        T* vtd = updated(vt);
        typename real_type<T>::type* sd  = updated(s);
        helper_lapack<T>::gesvd( "S", "S", &m, &n, ad, &m, sd, ud, &m, vtd, &k, &wkopt, &lwork, &info );
    }

    template<typename T>
    void geev<T>::c(const matrix<T>& a, unbound< matrix<T> >& lv, unbound< matrix<T> >& rv, unbound< matrix<T> >& s){
        int n = a.num_cols();
        int info;
        int lwork = -1;
        T wkopt;
        T* ad  = current(a);
        T* lvd = updated(lv);
        T* rvd = updated(rv);
        T* sd  = updated(s);
        helper_lapack<T>::geev("N", "V", &n, ad, &n, sd, lvd, &n, rvd, &n, &wkopt, &lwork, &info); 
    }

    template<typename T>
    void inverse<T>::c(matrix<T> & a){
        int info;
        int m = a.num_rows();
        int n = a.num_cols();
        T* ad = (T*)revised(a);    
        int* ipivd = new int[n];
        helper_lapack<T>::getrf(&m, &n, ad, &m, ipivd, &info);
        helper_lapack<T>::getri(&n, ad, &n, ipivd, &info);
        delete [] ipivd;
    }

    template<typename T>
    void heev<T>::c(matrix<T>& a, unbound< matrix<typename real_type<T>::type> >& w){
        int m = a.num_rows();
        int info, lwork = -1;
        T wkopt;
        T* work;
        T* ad = (T*)std::malloc(ambient::size(a));
        typename real_type<T>::type* wd = (typename real_type<T>::type*)std::malloc(ambient::size(w));
        std::memcpy(ad, (T*)current(a), ambient::size(a));
        std::memcpy(wd, (typename real_type<T>::type*)current(w), ambient::size(w));

        helper_lapack<T>::syev("V","U",&m,ad,&m,wd,&wkopt,&lwork,&info);

        typename real_type<T>::type s;
        for(int i=0; i < (int)(m/2); i++){ 
            s = wd[i];
            wd[i] = wd[m-i-1];
            wd[m-i-1] = s;
        } 
        // reversing eigenvectors
        size_t len = m*sizeof(T);
        work = (T*)std::malloc(len);
        for (int i=0; i < (int)(m/2); i++){ 
            std::memcpy(work, &ad[i*m], len);
            std::memcpy(&ad[i*m], &ad[(m-1-i)*m], len);
            std::memcpy(&ad[(m-1-i)*m], work, len);
        }
        std::memcpy((T*)updated(a), ad, ambient::size(a)); std::free(ad);
        std::memcpy((typename real_type<T>::type*)updated(w), wd, ambient::size(w));
        std::free(wd);
    }

} } }

#endif
