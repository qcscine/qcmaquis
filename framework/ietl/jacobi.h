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
 *               2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef IETL_JACOBI_H
#define IETL_JACOBI_H

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
    enum DesiredEigenvalue { Largest, Smallest };
    
    template <class MATRIX, class VS>
    class jcd_left_preconditioner
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        
        jcd_left_preconditioner(const MATRIX& matrix, const VS& vec, const int& max_iter);
        void operator()(const vector_type& u, const magnitude_type& theta, const vector_type& r, vector_type& t, const magnitude_type& rel_tol);
        
    private:
        void sysv(const char& uplo, fortran_int_t n, fortran_int_t nrhs, double a[], fortran_int_t lda, fortran_int_t ipiv[], double b[], fortran_int_t ldb, double work[], fortran_int_t lwork, fortran_int_t& info)
        { LAPACK_DSYSV(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info); };
        void sysv(const char& uplo, fortran_int_t n, fortran_int_t nrhs, std::complex<double> a[], fortran_int_t lda, fortran_int_t ipiv[], std::complex<double> b[], fortran_int_t ldb, std::complex<double> work[], fortran_int_t lwork, fortran_int_t& info)
        { LAPACK_ZHESV(&uplo, &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info); };
        MATRIX K;
        VS vecspace_;
        int n_;
        int max_iter_;
    };
    
    template <class MATRIX, class VS>
    class jcd_simple_solver
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        
        jcd_simple_solver(const MATRIX& matrix, const VS& vec);
        void operator()(const vector_type& u, const magnitude_type& theta, const vector_type& r, vector_type& t, const magnitude_type& rel_tol);
        
    private:
        MATRIX matrix_;
        VS vecspace_;
        int n_;
    };
    //
    // JCD_SOLVER_OPERATOR
    // -------------------
    // Class crehated indirectly by JCD_GMRES_SOLVER.
    // Has two functions, the () and the sum
    template<class Matrix, class VS, class Vector>
    class jcd_solver_operator
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        
        jcd_solver_operator(const vector_type& u,
            const magnitude_type& theta,
            const vector_type& r,
            const Matrix & m)
        : u_(u), r_(r), theta_(theta), m_(m) { }
        void operator()(vector_type const & x, vector_type & y) const;
    private:
        vector_type const & u_, r_;
        magnitude_type const & theta_;
        Matrix const & m_;
    };
    // Multiplication of solvers
    template<class Matrix, class VS, class Vector>
    void mult(jcd_solver_operator<Matrix, VS, Vector> const & m,
        typename jcd_solver_operator<Matrix, VS, Vector>::vector_type const & x,
        typename jcd_solver_operator<Matrix, VS, Vector>::vector_type & y)
    {
        m(x,y);
    }
    // Actual routine called when an object of JCD_SOLVER_OPERATOR is crehated
    template<class Matrix, class VS, class Vector>
    void jcd_solver_operator<Matrix, VS, Vector>::operator()(vector_type const & x, vector_type & y) const
    {
        // Calculate (1-uu*)(A-theta*1)(1-uu*), where u is the trial vector and
        // A is the operator to be diagonalized
        vector_type t, t2, t3;
        // t2 = (1-uu*) x
        scalar_type ust = dot(u_, x);
        t2 = x - ust*u_;
        // y = (A-theta*1) t2
        mult(m_, t2, t3);
        y = t3 - theta_*t2;
        // t = (1-uu*) y
        ust = dot(u_, y);
        t = y - ust*u_;
        y = t;
    }
    //
    // JCD_MODIFIED_SOLVER_OPERATOR
    // ----------------------------
    // Class crehated indirectly by JCD_GMRES_MODIFIED_SOLVER.
    // Has two functions, the () and the sum
    template<class Matrix, class VS, class Vector>
    class jcd_modified_solver_operator
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;

        jcd_modified_solver_operator(const vector_type& v,
                                     const magnitude_type& theta,
                                     const vector_type& r,
                                     const Matrix & m,
                                     const vector_type& z)
                : z_(z), r_(r), v_(v), theta_(theta), m_(m) { }
        void operator()(vector_type const & x, vector_type & y) const;
    private:
        vector_type const z_, r_, v_ ;
        magnitude_type const & theta_;
        Matrix const & m_;
    };
    // Multiplication of solvers
    template<class Matrix, class VS, class Vector>
    void mult(jcd_modified_solver_operator<Matrix, VS, Vector> const & m,
              typename jcd_modified_solver_operator<Matrix, VS, Vector>::vector_type const & x,
              typename jcd_modified_solver_operator<Matrix, VS, Vector>::vector_type & y)
    {
        m(x,y);
    }
    // Actual routine called when an object of JCD_SOLVER_OPERATOR is crehated
    template<class Matrix, class VS, class Vector>
    void jcd_modified_solver_operator<Matrix, VS, Vector>::operator()(vector_type const & x, vector_type & y) const
    {
        // Calculate (1-uu*)(A-theta*1)(1-uu*), where u is the trial vector and
        // A is the operator to be diagonalized
        vector_type t, t2, t3;
        // t2 = (1-uu*) x
        ietl::mult(m_, x, t);
        scalar_type ust = dot(z_, t);
        t2 = x - ust*v_;
        // y = (A-theta*1) t2
        mult(m_, t2, t3);
        y = t3 - theta_*t2;
        // t = (1-uu*) y
        ust = dot(z_, y);
        t = y - ust*z_ ;
        y = t;
    }
    //
    // JCD_GMRES_SOLVER
    // ----------------
    // kept for backward compatibility
    // NOTE : this is an "empty" class, when the object is called as jcd_gmres_solver(...),
    //        the other class
    template<class Matrix, class VS>
    class jcd_gmres_solver
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        // Constructor
        jcd_gmres_solver(Matrix const & matrix, VS const & vec,
            std::size_t max_iter = 5, bool verbose = false)
        : matrix_(matrix)
        , vecspace_(vec)
        , n_(vec_dimension(vec))
        , max_iter_(max_iter)
        , verbose_(verbose) { }
        // Function called during calculate_eigenvalue
        void operator()(const vector_type& u,
                        const magnitude_type& theta,
                        const vector_type& r, vector_type& t,
                        const magnitude_type& rel_tol)
        {
            jcd_solver_operator<Matrix, VS, vector_type> op(u, theta, r, matrix_);
            ietl_gmres gmres(max_iter_, verbose_);
            vector_type inh = -r;
            // initial guess for better convergence
            scalar_type dru = ietl::dot(r,u);
            scalar_type duu = ietl::dot(u,u);
            t = -r + dru/duu*u;
            if (max_iter_ > 0)
                t = gmres(op, inh, t, rel_tol);
        }
    private:
        Matrix const & matrix_;
        VS vecspace_;
        std::size_t n_, max_iter_;
        bool verbose_;
    };
    //
    // JCD_GMRES_MODIFIED_SOLVER
    // -------------------------
    // kept for backward compatibility
    // NOTE : this is an "empty" class, when the object is called as jcd_gmres_solver(...),
    //        the other class
    template<class Matrix, class VS>
    class jcd_gmres_modified_solver
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        // Constructor
        jcd_gmres_modified_solver(Matrix const & matrix, VS const & vec,
                                  std::size_t max_iter = 5, bool verbose = false)
                : matrix_(matrix)
                , vecspace_(vec)
                , n_(vec_dimension(vec))
                , max_iter_(max_iter)
                , verbose_(verbose) { }
        // Function called during calculate_eigenvalue
        void operator()(const vector_type& v,
                        const magnitude_type& theta,
                        const vector_type& r, vector_type& t,
                        const magnitude_type& rel_tol)
        {
            vector_type loc ;
            ietl::mult(matrix_ , v , loc);
            jcd_modified_solver_operator<Matrix, VS, vector_type> op(v, theta, r, matrix_, loc);
            ietl_gmres gmres(max_iter_, verbose_);
            vector_type inh = -r;
            // initial guess for better convergence
            scalar_type dru = ietl::dot(r,v);
            scalar_type duu = ietl::dot(v,v);
            t = -r + dru/duu*v;
            if (max_iter_ > 0)
                t = gmres(op, inh, v, rel_tol);
        }
    private:
        Matrix const & matrix_;
        VS vecspace_;
        std::size_t n_, max_iter_;
        bool verbose_;
    };
    
    template<class Matrix, class VS>
    class jcd_solver
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        
        template<class Solver>
        jcd_solver(Matrix const & matrix,
                   VS const & vec,
                   Solver solv)
        : matrix_(matrix)
        , vecspace_(vec)
        , solv_(solv)
        , n_(vec_dimension(vec)) { }
        
        template<class Solver>
        void replace_solver(Solver solv) { solv_ = solv; }
        
        void operator()(const vector_type& u,
                        const magnitude_type& theta,
                        const vector_type& r, vector_type& t,
                        const magnitude_type& rel_tol)
        {
            jcd_solver_operator<Matrix, VS, vector_type> op(u, theta, r, matrix_);
            
            vector_type inh = -r;
            
            // initial guess for better convergence
            scalar_type dru = ietl::dot(r,u);
            scalar_type duu = ietl::dot(u,u);
            t = -r + dru/duu*u;
            if (max_iter_ > 0)
                t = solv_(op, inh, t, rel_tol);
        }
        
    private:
        Matrix const & matrix_;
        VS vecspace_;
        boost::function<vector_type(jcd_solver_operator<Matrix, VS, vector_type> const &, vector_type const &, vector_type const &, double)> solv_;
        std::size_t n_, max_iter_;
        bool verbose_;
    };
    
    template <class MATRIX, class VS>
    jcd_simple_solver<MATRIX, VS>::jcd_simple_solver(const MATRIX& matrix, const VS& vec) :
    matrix_(matrix),
    vecspace_(vec)
    {
        n_ = vec_dimension(vecspace_);
    }
    
    template <class MATRIX, class VS>
    void jcd_simple_solver<MATRIX, VS>::operator()(const vector_type& u, const magnitude_type& theta, const vector_type& r, vector_type& t, const magnitude_type& rel_tol)
    {
       //for (int i=0;i<n_;i++)
         //  t[i] = -r[i] / ( matrix_(i,i) - theta );
        // std::cout << "Preconditioner, theta = " << theta << std::endl;
        // t = -1*r / (1-theta);
        t = -r + ietl::dot(r,u)/ietl::dot(u,u)*u;
    }
    
    template <class MATRIX, class VS>
    jcd_left_preconditioner<MATRIX, VS>::jcd_left_preconditioner(const MATRIX& matrix, const VS& vec, const int& max_iter) :
    K(matrix),
    vecspace_(vec),
    max_iter_(max_iter)
    {
        n_ = vecspace_.vec_dimension();
    }
    
    template <class MATRIX, class VS>
    void jcd_left_preconditioner<MATRIX, VS>::operator()(const vector_type& u, const magnitude_type& theta, const vector_type& r, vector_type& t, const magnitude_type& rel_tol)
    {
        // define variables
        FortranMatrix<scalar_type> Ktmp(n_,n_);
        vector_type u_hat = new_vector(vecspace_);
        vector_type vec1  = new_vector(vecspace_);
        vector_type vec2  = new_vector(vecspace_);
        magnitude_type mu, norm;
        
        // initialize variables
        for (int i=0;i<n_;i++) for (int j=0;j<n_;j++)
            Ktmp(i,j) = K(i,j);
        for (int i=0;i<n_;i++)
            Ktmp(i,i) -= theta;
        
        // define variables for LAPACK
        char uplo='U';
        fortran_int_t n=n_;
        fortran_int_t nrhs=1;
        fortran_int_t lda=n;
        fortran_int_t ldb=n;
        fortran_int_t lwork=n*n;
        fortran_int_t info;
        
        fortran_int_t *ipiv;
        ipiv = new fortran_int_t[n];
        scalar_type *work; 
        work = new scalar_type[lwork]; 
        
        // Solve u_hat from K*u_hat = u,  mu = u^star * u_hat
        ietl::copy(u,u_hat);
        sysv(uplo, n, nrhs, Ktmp.data(), lda, ipiv, u_hat.data(), ldb, work, lwork, info);
        mu = std::real(ietl::dot(u,u_hat));
        
        // compute r_tilde = K_tilde^{-1} * r as
        //   solve r_hat from K*r_hat = r
        //   r_tilde = r_hat - \frac{u^\star r_hat}{mu} u_hat
        ietl::copy(r,vec1);
        for (int i=0;i<n_;i++)
            Ktmp(i,i) = K(i,i);
        uplo='L';
        sysv(uplo, n, nrhs, Ktmp.data(), lda, ipiv, vec1.data(), ldb, work, lwork, info);
        vec1 -= ietl::dot(u,vec1)/mu*u_hat;
        vec1*=-1.;
        
        // aplly a Krylov subspace method with t_0=0; operator K_tilde ^{-1} A_tilde
        // and right hand side -r_tilde; given v, z=K_tilde ^{-1} A_tilde v is computed as
        //   y = (A-\theta I)v
        //   solve y_hat from K y_hat = y
        for (int i=0; i<max_iter_; i++) {
            t = vec1/ietl::two_norm(vec1);
            ietl::mult(K,t,vec1);
            for (int j=0;j<n_;j++) for (int k=0;k<n_;k++)
                Ktmp(j,k) = K(j,k);
            for (int j=0;j<n_;j++)
                Ktmp(j,j) -= theta;
            sysv(uplo, n, nrhs, Ktmp.data(), lda, ipiv, vec1.data(), ldb, work, lwork, info);
            vec1 -= -ietl::dot(u,vec1)/mu * u_hat;
            norm = ietl::dot(t,vec1);
            vec2 = vec1-norm*t;
            if ( ietl::two_norm(vec2) < rel_tol * std::abs(norm) )
                break;
        }
        
        delete [] ipiv;
        delete [] work;
    }
    //
    // Main class with Jacobi-Davidson
    // -------------------------------
    // Private attributes:
    // 1) matrix
    // 2) vector space
    // 3) n_, dimension of the vector space
    // 4) ???
    // 5) tolerance criteria
    // 6) which eigenvalue to compute
    //
    template <class MATRIX, class VS>
    class jacobi_davidson
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;


        jacobi_davidson(const MATRIX& matrix,
                        const VS& vec,
                        DesiredEigenvalue desired = Largest);
        jacobi_davidson(const MATRIX& matrix,
                        const VS& vec,
                        const double& omega,
                        DesiredEigenvalue desired = Largest);
        ~jacobi_davidson();

        template <class GEN, class SOLVER, class ITER>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen,
                                                                    SOLVER& solver,
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
        double omega_;
        bool has_omega_ ;
    };
    //
    // Methods of the Jacobi-Davidson class
    // ------------------------------------
    //
    // 1) constructor : just loads the vector space
    // 2) calculation of the eigenvectors (highest or lowest)
    // 3) get_extremal_eigenvalues : interfaces to FORTRAN routines to compute eigenvalues and eigenvectors
    //
    template <class MATRIX, class VS>
    jacobi_davidson<MATRIX, VS>::jacobi_davidson(const MATRIX& matrix, const VS& vec, DesiredEigenvalue desired) : 
    matrix_(matrix),
    vecspace_(vec),
    M(1,1),
    desired_(desired),
    has_omega_(false)
    {
        n_ = vec_dimension(vecspace_);
    }
    template <class MATRIX, class VS>
    jacobi_davidson<MATRIX, VS>::jacobi_davidson(const MATRIX& matrix, const VS& vec, const double& omega, DesiredEigenvalue desired) :
            matrix_(matrix),
            vecspace_(vec),
            M(1,1),
            desired_(desired),
            omega_(omega),
            has_omega_(false)
    {
        n_ = vec_dimension(vecspace_);
    }
    template <class MATRIX, class VS>
    jacobi_davidson<MATRIX, VS>::~jacobi_davidson() { }
    //
    // Method of the jacobi_davidson class to compute the eigenvalues
    // --------------------------------------------------------------
    template <class MATRIX, class VS> 
    template <class GEN, class SOLVER, class ITER>
    std::pair<typename jacobi_davidson<MATRIX,VS>::magnitude_type, typename jacobi_davidson<MATRIX,VS>::vector_type> 
    jacobi_davidson<MATRIX, VS>::calculate_eigenvalue(const GEN& gen,
                                                      SOLVER& solver,
                                                      ITER& iter)
    {
        // Variable declaration
        std::vector<scalar_type> s(iter.max_iterations());
        std::vector<vector_type> V(iter.max_iterations());
        std::vector<vector_type> VA(iter.max_iterations());
        M.resize(iter.max_iterations(), iter.max_iterations());
        //
        magnitude_type theta, tau, rel_tol;
        magnitude_type kappa = 0.25;
        magnitude_type shift = omega_ ;
        atol_ = iter.absolute_tolerance();
        //
        vector_type tA ;
        //
        ietl::generate(V[0],gen); const_cast<GEN&>(gen).clear();
        ietl::project(V[0],vecspace_);
        //
        // Main loop of the algorithm
        // --------------------------
        do {
            vector_type& t = V[iter.iterations()];
            // This part is basically the same as in the Jacobi method, and is based on the
            // orthogonalization of the new guess vector wrt the old ones.
            tau = ietl::two_norm(t);
            // The normalization is different for JCD and modified JCD algorithm
            if (has_omega_){
                for (int i = 1; i <= iter.iterations(); i++) {
                    t -= ietl::dot(VA[i-1], t) * V[i-1];
                    tA -= ietl::dot(VA[i-1], t) * VA[i-1];
                }
                t /= ietl::two_norm(tA);
                VA.push_back(tA/ietl::two_norm(tA));
            }
            else {
                for (int i = 1; i <= iter.iterations(); i++)
                    t -= ietl::dot(V[i-1], t) * V[i-1];
                if (ietl::two_norm(t) < kappa * tau)
                    for (int i = 1; i <= iter.iterations(); i++)
                        t -= ietl::dot(V[i-1], t) * V[i-1];
                t /= ietl::two_norm(t) ;
                ietl::mult(matrix_, t, VA[iter.iterations()]);
            }
            // TODO ALB commented for the moment ietl::project(t,vecspace_)
            // Update of the M matrix
            for(int i = 1; i <= iter.iterations()+1; i++)
                M(i-1,iter.iterations()) = ietl::dot(V[i-1], VA[iter.iterations()]);
            // compute the largest eigenpair (\theta, s) of M (|s|_2 = 1)
            get_extremal_eigenvalue(theta,s,iter.iterations()+1);
            // New guess vector is given in output as a function of the V basis, converted here
            // to the original basis. The same for H*v
            vector_type u = V[0] * s[0];
            for(int j = 1; j <= iter.iterations(); ++j)
                u += V[j] * s[j];
            vector_type uA = VA[0] * s[0];
            for(int j = 1; j <= iter.iterations(); ++j)
                uA += VA[j] * s[j];
            // TODO ALB commented for the moment ietl::project(uA,vecspace_);
            // Compute the error vector
            vector_type& r = uA; r -= theta*u;
            ++iter;
            if(iter.finished(ietl::two_norm(r),theta)) return std::make_pair(theta, u);
            // solve (approximately) a t orthogonal to u from
            //   (I-uu^\star)(A-\theta I)(I- uu^\star)t = -r
            rel_tol = 1. / pow(2.,double(iter.iterations()+1));
            if(has_omega_)
                solver(V[0], theta, r, V[iter.iterations()] , rel_tol);
            else
                solver(u, theta, r, V[iter.iterations()], rel_tol);
            //
            V[iter.iterations()].data().iter_index = VA[iter.iterations()-1].data().iter_index;
            storage::migrate(V[iter.iterations()], parallel::scheduler_balanced_iterative(V[iter.iterations()].data()));
        } while(true);
        
    }
    //
    // The following two matrices are interfaces to the Lapack routines to compute eigenvalues and eigenvectors
    // Overloaded to support both real and complex mumber.
    template <class MATRIX, class VS>
    void jacobi_davidson<MATRIX, VS>::get_extremal_eigenvalue(magnitude_type& theta, std::vector<double>& s, fortran_int_t dim)
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
    }
        
    template <class MATRIX, class VS>
    void jacobi_davidson<MATRIX, VS>::get_extremal_eigenvalue(magnitude_type& theta, std::vector<std::complex<double> >& s, fortran_int_t dim)
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
