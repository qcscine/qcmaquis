/*
 * Copyright Institute for Theoretical Physics, ETH Zurich 2014.
 * Distributed under the Boost Software License, Version 1.0.
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef AMBIENT_CONTAINER_NUMERIC_MATRIX_DETAIL_EXPERIMENTAL
#define AMBIENT_CONTAINER_NUMERIC_MATRIX_DETAIL_EXPERIMENTAL

#include "ambient/container/numeric/bindings/experimental.hpp"

namespace ambient { namespace numeric { namespace kernels {

    namespace detail {

        template<typename T, typename UL, typename OFF>
        void laset2(matrix<T>& a, const T& alfa){
            T* ad = a.data();
            plasma::lapack<T>::laset2(UL::value, a.num_rows()-OFF::value, a.num_cols()-OFF::value, alfa, ad + OFF::value*a.num_rows(), a.num_rows());
        }
        
        template<typename ALFA, typename T>
        void add_vectors(matrix<T>& a, const size_t& aoffset, const matrix<T>& b, const size_t& boffset, const size_t& size){
            const T* bd = b.data() + boffset;
            T* ar = a.data() + aoffset;
        
            for(int k = 0; k < size; k++) 
                ar[k] = ALFA::value*ar[k] + bd[k];
        }
            
        template<typename T>
        void labrd_update_col(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, matrix<T>& tq, matrix<T>& d, const int& i){
            static const double mone = -1.;
            static const double one = 1.;
            static const double zero = 0.;
            static const int lone = 1;
        
            int m  = num_rows(say);
            int n  = num_cols(sax);
            int ri = m-i;
            int rj = n-i-1;
        
            T* sayd = say.data(); int ldsay = say.num_rows();
      const T* saxd = sax.data(); int ldsax = sax.num_rows();
            T* syd  = sy.data();  int ldsy = sy.num_rows();
      const T* sxd  = sx.data();  int ldsx = sx.num_rows();
            T* tqd  = tq.data();
            T* dd   = d.data();
            
            if(i == 0){
                mkl::lapack_exp<T>::larfg(&ri, sayd, &sayd[1], &lone, tqd);
                *dd = *sayd;
                *sayd = 1.0;
                return;
            }
        
            ambient::memptf<T, ambient::memcpy>(sayd, ldsay, dim2(i, i-1), 
                                                saxd, ldsax, dim2(i, i-1), 
                                                dim2( num_cols(say)-i, 1));
        
            mkl::blas<T>::gemv("N", &ri, &i, &mone, &sayd[ i ], &ldsay, &syd[ i ], &ldsy, &one, &sayd[i + i*ldsay], &lone);
            mkl::blas<T>::gemv("N", &ri, &i, &mone, &sxd[ i ], &ldsx, &sayd[ i*ldsay ], &lone, &one, &sayd[i + i*ldsay], &lone);
            
            mkl::lapack_exp<T>::larfg( &ri, &sayd[i+i*ldsay], &sayd[std::min(i+1, m-1)+i*ldsay], &lone, &tqd[i] );
            dd[i] = sayd[i+i*ldsay];
            sayd[i+i*ldsay] = 1.000;
        }
        
        template<typename T, typename IB>
        void labrd_reduce_col(matrix<T>& say, const matrix<T>& sax, matrix<T>& sy, const matrix<T>& sx, const int& i){
            static const double mone = -1.;
            static const double one = 1.;
            static const double zero = 0.;
            static const int lone = 1;
        
            int m  = num_rows(say);
            int n  = num_cols(sax);
            int ri = m-i;
            int rj = n-i-1;
            int ari = IB::value-i-1;
        
            T* sayd = say.data(); int ldsay = say.num_rows();
      const T* saxd = sax.data(); int ldsax = sax.num_rows();
            T* syd  = sy.data();  int ldsy = sy.num_rows();
      const T* sxd  = sx.data();  int ldsx = sx.num_rows();
            
            mkl::blas<T>::gemv("T", &ri, &ari, &one, &sayd[i + (i+1)*ldsay], &ldsay, &sayd[i+i*ldsay], &lone, &zero, &syd[i+1 + i*ldsy], &lone); // part of big gemv
        
            mkl::blas<T>::gemv("T", &ri, &i, &one, &sayd[i], &ldsay, &sayd[i+i*ldsay], &lone, &zero, &syd[i*ldsy], &lone);
            mkl::blas<T>::gemv("N", &rj, &i, &mone, &syd[ i+1 ], &ldsy, &syd[i*ldsy], &lone, &one, &syd[i+1 + i*ldsy], &lone);
        
            mkl::blas<T>::gemv("T", &ri, &i, &one, &sxd[i], &ldsx, &sayd[i + i*ldsay], &lone, &zero, &syd[ i*ldsy ], &lone);
            mkl::blas<T>::gemv("T", &i, &rj, &mone, &saxd[ (i+1)*ldsax ], &ldsax, &syd[i*ldsy], &lone, &one, &syd[ i+1 + i*ldsy], &lone);
        }
        
        template<typename T>
        void labrd_update_row(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, matrix<T>& tp, matrix<T>& e, const int& i){
            static const double mone = -1.;
            static const double one = 1.;
            static const double zero = 0.;
            static const int lone = 1;
        
            int m   = num_rows(say);
            int n   = num_cols(sax);
            int ri  = m-i;
            int rj  = n-i-1;
            int rij = m-i-1;
            int r3  = i+1;
        
      const T* sayd = say.data(); int ldsay = say.num_rows();
            T* saxd = sax.data(); int ldsax = sax.num_rows();
      const T* syd  = sy.data();  int ldsy = sy.num_rows();
            T* sxd  = sx.data();  int ldsx = sx.num_rows();
            T* tpd  = tp.data();
            T* ed   = e.data();
            
            ambient::memptf<T, ambient::memcpy>(saxd, ldsax, dim2(i, i), 
                                                sayd, ldsay, dim2(i, i), 
                                                dim2( 1, ldsax-i ));
            
            mkl::blas<T>::gemv("T", &i, &rj, &mone, &saxd[(i+1)*ldsax], &ldsax, &sxd[i], &ldsx, &one, &saxd[ i + (i+1)*ldsax], &ldsax);
            mkl::blas<T>::gemv("N", &rj, &r3, &mone, &syd[ i+1 ], &ldsy, &saxd[i], &ldsax, &one, &saxd[i + (i+1)*ldsax], &ldsax);
        
            mkl::lapack_exp<T>::larfg(&rj, &saxd[i + (i+1)*ldsax], &saxd[i + std::min(i+2, n-1)*ldsax], &ldsax, &tpd[i] );
            ed[i] = saxd[i + (i+1)*ldsax];
            saxd[i + (i+1)*ldsax] = 1.000;
        }
        
        template<typename T, typename IB>
        void labrd_reduce_row(const matrix<T>& say, matrix<T>& sax, const matrix<T>& sy, matrix<T>& sx, const int& i){
            static const double mone = -1.;
            static const double one = 1.;
            static const double zero = 0.;
            static const int lone = 1;
        
            int m   = num_rows(say);
            int n   = num_cols(sax);
            int ri  = m-i;
            int rj  = n-i-1;
            int rij = m-i-1;
            int r3  = i+1;
            int ari = IB::value-i-1;
        
      const T* sayd = say.data(); int ldsay = say.num_rows();
            T* saxd = sax.data(); int ldsax = sax.num_rows();
      const T* syd  = sy.data();  int ldsy = sy.num_rows();
            T* sxd  = sx.data();  int ldsx = sx.num_rows();
            
            mkl::blas<T>::gemv("T", &rj, &r3, &one, &syd[i+1], &ldsy, &saxd[i+(i+1)*ldsax], &ldsax, &zero, &sxd[i*ldsx], &lone);
            mkl::blas<T>::gemv("N", &rij, &r3, &mone, &sayd[i+1], &ldsay, &sxd[i*ldsx], &lone, &zero, &sxd[i+1+i*ldsx], &lone);
        
            mkl::blas<T>::gemv("N", &i, &rj, &one, &saxd[(i+1)*ldsax], &ldsax, &saxd[ i +(i+1)*ldsax], &ldsax, &zero, &sxd[i*ldsx], &lone);
            mkl::blas<T>::gemv("N", &rij, &i, &mone, &sxd[i+1], &ldsx, &sxd[i*ldsx], &lone, &one, &sxd[i+1+i*ldsx], &lone);
        
            mkl::blas<T>::gemv("N", &ari, &rj, &one, &saxd[i+1 + (i+1)*ldsax], &ldsax, &saxd[ i +(i+1)*ldsax], &ldsax, &one, &sxd[i+1 + i*ldsx], &lone); // part of big gemv
        }
        
        template<typename T, typename TR>
        void larfg(matrix<T>& a, matrix<T>& t, matrix<T>& d, const size_t& k){
            int lda;
            int n;
            T* alfa;
            T* x;
            
            T* ad = a.data();
            T* td = t.data();
            T* dd = d.data();
        
            if(TR::value == PlasmaNoTrans){
                alfa = &ad[k + k*a.num_rows()];
                x = &ad[std::min(k+1, a.num_rows()-1)+k*a.num_rows()];
                n = a.num_rows()-k;
                lda = 1;
            }else{
                alfa = &ad[k + (k+1)*a.num_rows()];
                x = &ad[k + std::min(k+2, a.num_cols()-1)*a.num_rows()];
                n = a.num_cols()-k-1;
                lda = a.num_rows();
            }
        
            mkl::lapack_exp<T>::larfg(&n, alfa, x, &lda, &td[k]);
            
            dd[k] = *alfa;
            *alfa = 1.00;
        }
        
        template<typename T>
        void gebd2(matrix<T>& a, unbound< matrix<T> >& d, unbound< matrix<T> >& e, unbound< matrix<T> >& tq, unbound< matrix<T> >& tp){
            int m = a.num_rows();
            int n = a.num_cols();
            int lda = a.num_rows();
            int info;
        
            T* work = (T*)std::malloc(std::max(m,n)*sizeof(T));
            mkl::lapack_exp<T>::gebd2(&m, &n, a.data(), &lda, d.data(), e.data(), tq.data(), tp.data(), work, &info);
            std::free(work);
        }
        
        template<typename T>
        void gebrd(matrix<T>& a, unbound< matrix<T> >& d, unbound< matrix<T> >& e, unbound< matrix<T> >& q, unbound< matrix<T> >& p){
            static int zero = 0;
            static int one  = 1;
            int m = q.num_rows();
            int n = p.num_rows();
            int lda = a.num_rows();
            int k = std::min(m,n);
            int lwork = -1;
            int info;
            
            T* ad = a.data();
            T* dd = d.data();
            T* ed = e.data();
            T* qd = q.data();
            T* pd = p.data();
        
            T* work = (T*)std::malloc(sizeof(T));
            T* tauq = (T*)std::malloc(sizeof(T)*k); // leak
            T* taup = (T*)std::malloc(sizeof(T)*k); // leak
        
            mkl::lapack_exp<T>::gebrd(&m, &n, ad, &lda, dd, ed, tauq, taup, work, &lwork, &info);
            lwork = (int)work[0];
            work = (T*)realloc(work, sizeof(T)*lwork);
            mkl::lapack_exp<T>::gebrd(&m, &n, ad, &lda, dd, ed, tauq, taup, work, &lwork, &info);
        
            T* ac = (T*)std::malloc(m*n*sizeof(T));
            std::memcpy(ac, ad, m*n*sizeof(T));
        
            lwork = -1;
            mkl::lapack_exp<T>::orgbr("Q",&m,&k,&n, ad, &lda, tauq, work, &lwork, &info);
            lwork = (int)work[0];
            work = (T*)realloc(work, sizeof(T)*lwork);
            mkl::lapack_exp<T>::orgbr("Q",&m,&k,&n, ad, &lda, tauq, work, &lwork, &info);
        
        
            lwork = -1;
            mkl::lapack_exp<T>::orgbr("P",&k,&n,&m, ac, &lda, taup, work, &lwork, &info);
            lwork = (int)work[0];
            work = (T*)realloc(work, sizeof(T)*lwork);
            mkl::lapack_exp<T>::orgbr("P",&k,&n,&m, ac, &lda, taup, work, &lwork, &info);
        
            for(int j = 0; j < n; ++j){
                std::memcpy(&pd[n*j], &ac[lda*j], sizeof(T)*k);
            }
        
            std::free(work); std::free(ac);
        }
        
        template<typename T>
        void gbbrd(matrix<T>& a, unbound< matrix<T> >& d, unbound< matrix<T> >& e, unbound< matrix<T> >& q, unbound< matrix<T> >& p){
            static int zero = 0;
            static int one  = 1;
            int m = q.num_rows();
            int n = p.num_rows();
            int k = a.num_rows();
            int kl = (m <  n) ? k-1 : 0;
            int ku = (m >= n) ? k-1 : 0;
            int info;
            
            T* work = (T*)std::malloc(std::max(m,n)*2*sizeof(T));
            T* ad = a.data();
            T* dd = d.data();
            T* ed = e.data();
            T* qd = q.data();
            T* pd = p.data();
        
            mkl::lapack_exp<T>::gbbrd("B", &m, &n, &zero, &kl, &ku, ad, &k, dd, ed, 
                    qd, &m, pd, &n, NULL, &one, work, &info);
        
            std::free(work);
        }
        
        template<typename T>
        void bdsqr(matrix<T>& d, matrix<T>& e, matrix<T>& u, matrix<T>& v){
            static int zero = 0;
            static int one  = 1;
            int n = d.num_rows();
            int nv = v.num_cols();
            int mv = v.num_rows();
            int mu = u.num_rows();
            int info;
            
            T* work = (T*)std::malloc(n*4*sizeof(T));
            mkl::lapack_exp<T>::bdsqr("U", &n, &nv, &mu, &zero, d.data(), e.data(), 
                                    v.data(), &mv, u.data(), &mu, NULL, &one, work, &info); std::free(work);
            // Uncomment for dc numerically loose algorithm:
            // LAPACKE_dbdsdc(102, 'U', 'N', n, d.data(), e.data(),  u.data(), one, v.data(), one, NULL, NULL);
        }
        
        template<typename T, typename UL, typename IB>
        void copy_band(const matrix<T>& src, matrix<T>& dst, const size_t& dj){
      const T* sd = src.data();
            T* dd = dst.data();
            size_t ldd = dst.num_rows();
            size_t m = src.num_rows();
            size_t n = src.num_cols();
            size_t offset = std::min(ldd-1,(size_t)IB::value);
        
            dd += dj*ldd;
            if(UL::value == PlasmaUpper){
                for(int j = 0; j < n; ++j)
                for(int i = 0; i <= j && i < m; ++i)
                dd[j*ldd+i-j+offset] = sd[i+m*j]; 
            }else{
                for(int j = 0; j < n; ++j)
                for(int i = j; i < m; ++i)
                dd[j*ldd+i-j] = sd[i+m*j]; 
            }
        }
        
        template<typename T, typename ADD, class VA, class VB, class VC, class VF>
        void gemv_scale(const matrix<T>& a, const size_t& aoffset, 
                        const matrix<T>& b, const size_t& boffset,
                              matrix<T>& c, const size_t& coffset,
                        const matrix<T>& f, const size_t& foffset,
                        const size_t& rows, const size_t& cols)
        {
            const T* ad = a.data();
            const T* bd = b.data();
            const T* fd = f.data();
            T* cd = c.data();
            int lda = num_rows(a);
            int ldb = VB::inc(b);
            int ldc = VC::inc(c);
            int m = rows;
            int n = cols;
            static const double one = 1.;
            const double* alfa = &fd[foffset];
            const double* beta = (ADD::value == 1) ? &one : alfa;
        
            if(*VA::code() == 'T') std::swap(m,n);
            mkl::blas<T>::gemv(VA::code(), &m, &n, alfa, &ad[aoffset], &lda, &bd[boffset], &ldb, beta, &cd[coffset], &ldc);
        }
        
        template<typename T, typename ALFA, typename BETA, class ViewA, class ViewB, class ViewC>
        void gemv(const matrix<T>& a, const size_t& aoffset, 
                  const matrix<T>& b, const size_t& boffset,
                        matrix<T>& c, const size_t& coffset,
                  const size_t& rows, const size_t& cols)
        {
            const T* ad = a.data();
            const T* bd = b.data();
            T* cd = c.data();
            int lda = num_rows(a);
            int ldb = ViewB::inc(b);
            int ldc = ViewC::inc(c);
            static const double salfa(ALFA::value); 
            static const double sbeta(BETA::value);
            int m = rows;
            int n = cols;
        
            if(*ViewA::code() == 'T') std::swap(m,n);
            mkl::blas<T>::gemv(ViewA::code(), &m, &n, &salfa, &ad[aoffset], &lda, &bd[boffset], &ldb, &sbeta, &cd[coffset], &ldc);
        }
        
        template<typename T>
        void norm_vector(const matrix<T>& a, unbound< matrix<typename real_type<T>::type> >& b){
            int m = num_rows(a);
            int n = num_cols(a);
            const T* ad = a.data();
            typename real_type<T>::type* bd = b.data();
            std::cout << " (m,n) " << m << " " << n << std::endl;
            for(int i(0); i<n; ++i)
                for(int j(0); j<m; ++j) // writing garbage, aren't we? :)
                    bd[i] += mkl::helper_complex<T>::real(ad[i*m+j]*mkl::helper_complex<T>::conj(ad[i*m+j])); 
        }   
            
        template<typename T>
        void max_vector(const matrix<typename real_type<T>::type>& a, future<typename real_type<T>::type>& ret){
            const typename real_type<T>::type* ad = a.data();
            typename real_type<T>::type tmp = ad[0];
            int n = num_cols(a);
                
            for(int i(0); i<n; ++i)
                tmp = std::max(tmp,ad[i]);
            ret.set(tmp);
        }
        
        template<typename T>
        void sqrt_inplace(matrix<typename real_type<T>::type>& a){
            size_t size = a.num_rows()*a.num_cols();
            typename real_type<T>::type* ar = a.data();
            for(int i=0; i < size; ++i)
                ar[i] =sqrt(ar[i]);   
        }
        
        template<typename T>
        void init_gaussian(unbound< matrix<T> >& a){
            size_t size = ambient::get_square_dim(a);
            T* ad = a.data();
            // for(size_t i = 0; i < size; ++i) ad[i] = gaussian_generator(rng);
        }

    }

    AMBIENT_EXPORT_TEMPLATE(detail::laset2, laset2)
    AMBIENT_EXPORT_TEMPLATE(detail::add_vectors, add_vectors)
    AMBIENT_EXPORT_TEMPLATE(detail::labrd_update_col, labrd_update_col)
    AMBIENT_EXPORT_TEMPLATE(detail::labrd_reduce_col, labrd_reduce_col)
    AMBIENT_EXPORT_TEMPLATE(detail::labrd_update_row, labrd_update_row)
    AMBIENT_EXPORT_TEMPLATE(detail::labrd_reduce_row, labrd_reduce_row)
    AMBIENT_EXPORT_TEMPLATE(detail::larfg, larfg)
    AMBIENT_EXPORT_TEMPLATE(detail::gebd2, gebd2)
    AMBIENT_EXPORT_TEMPLATE(detail::gebrd, gebrd)
    AMBIENT_EXPORT_TEMPLATE(detail::gbbrd, gbbrd)
    AMBIENT_EXPORT_TEMPLATE(detail::bdsqr, bdsqr)
    AMBIENT_EXPORT_TEMPLATE(detail::copy_band, copy_band)
    AMBIENT_EXPORT_TEMPLATE(detail::gemv_scale, gemv_scale)
    AMBIENT_EXPORT_TEMPLATE(detail::gemv, gemv)
    AMBIENT_EXPORT_TEMPLATE(detail::norm_vector, norm_vector)
    AMBIENT_EXPORT_TEMPLATE(detail::max_vector, max_vector)
    AMBIENT_EXPORT_TEMPLATE(detail::sqrt_inplace, sqrt_inplace)
    AMBIENT_EXPORT_TEMPLATE(detail::init_gaussian, init_gaussian)

} } }

#endif
