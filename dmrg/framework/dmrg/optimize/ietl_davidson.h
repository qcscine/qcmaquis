/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2001-2011 by Rene Villiger <rvilliger@smile.ch>,
 *                            Prakash Dayal <prakash@comp-phys.org>,
 *                            Matthias Troyer <troyer@comp-phys.org>
 *                            Bela Bauer <bauerb@phys.ethz.ch>
 *
 * This software is part of the ALPS libraries, published under the ALPS
 * Library License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Library License along with
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

/* $Id: jacobi.h,v 1.6 2003/09/05 08:12:38 troyer Exp $ */

#ifndef IETL_DAVIDSON_H
#define IETL_DAVIDSON_H

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h> 
 
#include <ietl/cg.h>
#include <ietl/gmres.h>

#include <complex>
#include <vector>

#include <boost/function.hpp>

namespace ietl
{
    //enum DesiredEigenvalue { Largest, Smallest };
    
    template <class MATRIX, class VS>
    class davidson
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        
        
        davidson(const MATRIX& matrix, 
                 const VS& vec,
                 DesiredEigenvalue desired = Largest);
        ~davidson();
        
        template <class GEN, class SOLVER, class PRECOND, class ITER>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen, 
                                                                    SOLVER& solver,
                                                                    PRECOND& mdiag,
                                                                    ITER& iter);
        
    private:
        void get_extremal_eigenvalue(magnitude_type& theta, std::vector<double>& s, fortran_int_t dim);
        void get_extremal_eigenvalue(magnitude_type& theta, std::vector<std::complex<double> >& s, fortran_int_t dim);
        MATRIX const & matrix_;
        VS vecspace_;
        int n_;
        FortranMatrix<scalar_type> M;
        magnitude_type atol_;
        DesiredEigenvalue desired_;
    };
    
    // C L A S S :   J A C O B I _ D A V I D S O N ////////////////////////////////////
    
    template <class MATRIX, class VS>
    davidson<MATRIX, VS>::davidson(const MATRIX& matrix, const VS& vec, DesiredEigenvalue desired) : 
    matrix_(matrix),
    vecspace_(vec),
    M(1,1),
    desired_(desired)
    {
        n_ = vec_dimension(vecspace_);
    }
    
    template <class MATRIX, class VS>
    davidson<MATRIX, VS>::~davidson()
    {
        
    }
    
//    template<class Vector>
//    Vector orthogonalize(Vector input, std::vector<Vector> const & against)
//    {
//        
//    }
    
    template <class MATRIX, class VS> 
    template <class GEN, class SOLVER, class PRECOND, class ITER>
    std::pair<typename davidson<MATRIX,VS>::magnitude_type, typename davidson<MATRIX,VS>::vector_type> 
    davidson<MATRIX, VS>::calculate_eigenvalue(const GEN& gen,
                                               SOLVER& solver,
                                               PRECOND& mdiag,
                                               ITER& iter)
    {
        typedef alps::numeric::matrix<scalar_type> matrix_t;

        vector_type t  = new_vector(vecspace_);
        vector_type u  = new_vector(vecspace_);
        vector_type uA = new_vector(vecspace_);
        // vector_type vA = new_vector(vecspace_);
        vector_type r  = new_vector(vecspace_);

        std::vector<scalar_type> s(iter.max_iterations());
        std::vector<vector_type> V;
        std::vector<vector_type> VA;

        unsigned int i,j;
        M.resize(iter.max_iterations(), iter.max_iterations());
        magnitude_type theta, tau;
        magnitude_type kappa = 0.25;
        magnitude_type rel_tol;
        atol_ = iter.absolute_tolerance();
        
        // Start with t=v_o, starting guess
        ietl::generate(t,gen);
        ietl::project(t,vecspace_);
        
        // Start iteration
        do
        {
            // Modified Gram-Schmidt Orthogonalization with Refinement
            tau = ietl::two_norm(t);
            for (i = 0; i < V.size(); i++)
                t -= ietl::dot(V[i], t) * V[i];
            if (ietl::two_norm(t) < kappa * tau)
                for (i = 0; i < V.size(); i++)
                    t -= ietl::dot(V[i], t) * V[i];
            
            // Project out orthogonal subspace
            ietl::project(t, vecspace_);
            
            // v_m = t / |t|_2,  v_m^A = A v_m
            V.push_back(t/ietl::two_norm(t));
            VA.resize(V.size());
            ietl::mult(matrix_, V[V.size() - 1], VA[V.size() - 1]);
            
            
            /////////////////////////////////////////////////7
 
            //  for (i=1; i<=iter.iterations()+1; i++)
            //      M(i-1,iter.iterations()) = ietl::dot(V[i-1], VA[iter.iterations()]);
            //  
            //  // compute the largest eigenpair (\theta, s) of M (|s|_2 = 1)
            //  get_extremal_eigenvalue(theta, s, iter.iterations()+1);
            //  
            //  // u = V s
            //  u = V[0] * s[0];
            //  for (j=1;j<=iter.iterations();j++)
            //      u += V[j] * s[j];
            //  
            //  // u^A = V^A s
            //  // ietl::mult(matrix_,u,uA);
            //  uA = VA[0] * s[0];
            //  for (j=1;j<=iter.iterations();++j)
            //      uA += VA[j] * s[j];
            //  
            //  ietl::project(uA, vecspace_);
            
            //  // r = u^A - \theta u
            //  r = uA-theta*u;
            //  //std::cout << "Iteration " << iter.iterations() << ", resid = " << ietl::two_norm(r) << std::endl;
            //  
            //  // if (|r|_2 < \epsilon) stop
            //  ++iter;
            //  if (iter.finished(ietl::two_norm(r),theta))
            //      break;
            //  
            //  // solve (approximately) a t orthogonal to u from
            //  //   (I-uu^\star)(A-\theta I)(I- uu^\star)t = -r
            //  rel_tol = 1. / pow(2.,double(iter.iterations()+1));
            //  vector_type told = t;
            //  solver(u, theta, r, t, rel_tol);

            //  // std::cout << "Orthogonal? " << ietl::dot(t, u) << std::endl;
            //  // std::cout << "Expansion? " << ietl::dot(t, told) << std::endl;
            /////////////////////////////////////////////////7

            std::size_t iter_dim = V.size();
            matrix_t Mp(iter_dim, iter_dim), Mevecs(iter_dim, iter_dim);
            std::vector<magnitude_type> Mevals(iter_dim);

            for (i = 0; i < iter_dim; ++i)
            for (j = i; j < iter_dim; ++j)
            {
                Mp(i,j) = ietl::dot(V[i], VA[j]);
                Mp(j,i) = Mp(i,j);
            }

            boost::numeric::bindings::lapack::heevd('V', Mp, Mevals);
            Mevecs = Mp;

            std::vector<vector_type> Vp = V, VAp = VA;
            for (i = 0; i < iter_dim; ++i)
                V[i] *= Mevecs(i,i);
            for (i = 0; i < iter_dim; ++i)
                VA[i] *= Mevecs(i,i);

            for (i = 0; i < iter_dim; ++i)
            for (j = 0; j < iter_dim; ++j)
                if(i != j)
                {
                    V[j] += Vp[i] * Mevecs(i,j);
                    VA[j] += VAp[i] * Mevecs(i,j);
                }

            r = VA[0] - V[0] * Mevals[0];

            theta = Mevals[0];
            u = V[0];

            //  // if (|r|_2 < \epsilon) stop
            ++iter;
            if (iter.finished(ietl::two_norm(r), Mevals[0]))
                break;

            mdiag.precondition(r, V[0], Mevals[0]);
            std::swap(t,r);  
            t /= ietl::two_norm(t);

            if (V.size() >= 20)
            {
                V.resize(2);
                VA.resize(2);
            }

            /////////////////////////////////////////////////7

        } while (true);
        
        // accept lambda=theta and x=u
        return std::make_pair(theta, u);
    }
    
    template <class MATRIX, class VS>
    void davidson<MATRIX, VS>::get_extremal_eigenvalue(magnitude_type& theta, std::vector<double>& s, fortran_int_t dim)
    {
        FortranMatrix<scalar_type> M_(dim,dim);
        for (int i=0;i<dim;i++) for (int j=0;j<=i;j++)
            M_(j,i) = M(j,i);
        double abstol = atol_;
        char jobz='V';     char range='I';   char uplo='U';
        fortran_int_t n=dim;
        fortran_int_t lda=dim;       
        fortran_int_t il, iu;
        if (desired_ == Largest)
            il = iu = n;
        else
            il = iu = 1;

        fortran_int_t m;
        fortran_int_t ldz=n;
        fortran_int_t lwork=8*n;
        fortran_int_t info;

        double vl, vu;

        double *w = new double[n];
        double *z = new double[n];
        double *work = new double[lwork];
        fortran_int_t *iwork = new fortran_int_t[5*n];
        fortran_int_t *ifail = new fortran_int_t[n];

        LAPACK_DSYEVX(&jobz, &range, &uplo, &n, M_.data(), &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);
        
        theta = w[0];
        for (int i=0;i<n;i++)
            s[i] = z[i];
        
        delete [] w;
        delete [] z;
        delete [] work;
        delete [] iwork;
        delete [] ifail;

        //typedef alps::numeric::matrix<scalar_type> matrix_t;
        //matrix_t M_(dim, dim), evecs_(dim, dim);
        //std::vector<magnitude_type> evals_(dim);
        //for (int i=0; i < dim; i++)
        //for (int j=i; j < dim; j++)
        //{
        //    M_(i,j) = M(i,j);
        //    M_(j,i) = M(i,j);
        //}

        //heev(M_, evecs_, evals_);
        //theta = evals_[dim-1];
        //std::copy(evecs_.col(dim-1).first, evecs_.col(dim-1).second, s.begin());
    }
        
    template <class MATRIX, class VS>
    void davidson<MATRIX, VS>::get_extremal_eigenvalue(magnitude_type& theta, std::vector<std::complex<double> >& s, fortran_int_t dim)
    {
        FortranMatrix<scalar_type> M_(dim,dim);
        for (int i=0;i<dim;i++) for (int j=0;j<=i;j++)
            M_(j,i) = M(j,i);
        double abstol = atol_;
        char jobz='V';     char range='I';   char uplo='U';
        fortran_int_t n=dim;
        fortran_int_t lda=dim;
        fortran_int_t il, iu;
        if (desired_ == Largest)
            il = iu = n;
        else
            il = iu = 1;
        fortran_int_t m;
        fortran_int_t ldz=n;
        fortran_int_t lwork=8*n;
        fortran_int_t info; 
        double vl, vu;

        double * w = new double[n];
        std::complex<double> * z = new std::complex<double>[n];
        std::complex<double> * work = new std::complex<double>[lwork];
        fortran_int_t *iwork = new fortran_int_t[5*n];
        fortran_int_t *ifail = new fortran_int_t[n];
        double * rwork = new double[7*n];

        LAPACK_ZHEEVX(&jobz, &range, &uplo, &n, M_.data(), &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, &lwork, rwork, iwork, ifail, &info);
            
        theta = w[0];
        for (int i=0;i<n;i++)
            s[i] = z[i];
        
        
        delete [] w;
        delete [] z;
        delete [] work;
        delete [] iwork;
        delete [] ifail;
        delete [] rwork;
    }
}
#endif
