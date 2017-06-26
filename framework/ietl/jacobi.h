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

#include "dmrg/optimize/partial_overlap.h"

namespace ietl
{
    //
    // JCD_SOLVER_OPERATOR
    // -------------------
    // Class crehated indirectly by JCD_GMRES_SOLVER.
    //
    // Attributes:
    // 1) vectors u_ and r_ (projector and error vectors)
    // 2) matrix m_ to be diagonalized
    // 3) theta_ , approximation to the eigenvalue
    // Methods:
    // jcd_sovler_operator(x,y)     -->   takes in input a vector x and solve the approximate
    //                                    JD equation to give a new estimate for the vector y
    // m = mult(a,b), with a and b  -->   as previously, but m is specified explicitly
    // vectors and m jcd_solver_operator
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
    // Actual routine called when an object of JCD_SOLVER_OPERATOR is created
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
                                     const vector_type& z,
                                     const double omega )
                : z_(z), r_(r), v_(v), theta_(theta), m_(m), omega_(omega) { }
        void operator()(vector_type const & x, vector_type & y) const;
    private:
        vector_type const z_, r_, v_ ;
        magnitude_type const & theta_;
        Matrix const & m_;
        double omega_;
    };
    // Multiplication of solvers
    template<class Matrix, class VS, class Vector>
    void mult(jcd_modified_solver_operator<Matrix, VS, Vector> const & m,
              typename jcd_modified_solver_operator<Matrix, VS, Vector>::vector_type const & x,
              typename jcd_modified_solver_operator<Matrix, VS, Vector>::vector_type & y)
    {
        m(x,y);
    }
    // Actual routine called when an object of JCD_SOLVER_OPERATOR is created
    template<class Matrix, class VS, class Vector>
    void jcd_modified_solver_operator<Matrix, VS, Vector>::operator()(vector_type const & x, vector_type & y) const
    {
        //
        //  Solve the non-linear equation
        //  -----------------------------
        //
        //  (I - z * z^*) ( A - \theta * I ) (I - v * z^* * A) x = -r
        //
        //  where:
        //  1) z = A*v, where v is the lowest eigenvector
        //  2) A is the matrix of the site problem
        //  3) \theta_k is the harmonic Ritz value
        //  4) v is the lowest-energy eigenvector
        //  5) r is the error vector
        //
        vector_type t, t2, t3;
        ietl::mult(m_, x, t);
        t *= -1. ;
        t += omega_*x ;
        scalar_type ust = dot(z_, t);
        t2 = x - ust*v_;
        // y = (A-theta*1) t2
        mult(m_, t2, t3);
        t3 *= -1. ;
        t3 += omega_*t2 ;
        y = t3 - t2/theta_ ;
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
            std::size_t max_iter = 5, bool verbose = true)
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
            vector_type inh = r;
            // initial guess for better convergence
            scalar_type dru = ietl::dot(r,u);
            scalar_type duu = ietl::dot(u,u);
            t = r + dru/duu*u;
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
        jcd_gmres_modified_solver(Matrix const & matrix, VS const & vec, double omega,
                                  std::size_t max_iter = 5, bool verbose = true)
                : matrix_(matrix)
                , vecspace_(vec)
                , n_(vec_dimension(vec))
                , max_iter_(max_iter)
                , omega_(omega)
                , verbose_(verbose) { }
        // Function called during calculate_eigenvalue
        void operator()(const vector_type& v,
                        const magnitude_type& theta,
                        const vector_type& r, vector_type& t,
                        const magnitude_type& rel_tol)
        {
            vector_type loc ;
            ietl::mult(matrix_ , v , loc);
            loc *= -1 ;
            loc += v*omega_;
            jcd_modified_solver_operator<Matrix, VS, vector_type> op(v, theta, r, matrix_, loc, omega_);
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
        double omega_;
        bool verbose_;
    };

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
    template <class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    class jacobi_davidson
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename std::vector<vector_type> vector_space;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        typedef typename partial_overlap<OtherMatrix,SymmGroup>::partial_overlap partial_overlap;
        // The constructor is overloaded depending if the omega parameter is set in input
        // or not.
        jacobi_davidson(const MATRIX& matrix, const VS& vec, const partial_overlap& poverlap, const int& site,
                        const int& n_mo);
        jacobi_davidson(const MATRIX& matrix, const VS& vec, const double& omega, const partial_overlap& poverlap,
                        const int& site, const int& n_mo);
        ~jacobi_davidson();
        template <class GEN, class SOLVER, class ITER>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen,
                                                                    SOLVER& solver,
                                                                    ITER& iter);

    private:
        // Multiply the vector by a matrix
        vector_type apply_operator (const vector_type& x);
        void update_vecspace(vector_space &V, vector_space &VA, const int i);
        template <class ITER >
        bool check_convergence(const vector_type& u, const vector_type& uA , const magnitude_type theta ,
                               ITER& iter, vector_type& eigvec, magnitude_type& eigval);
        vector_type compute_error (const vector_type& u , const vector_type& uA,
                                   magnitude_type theta) ;
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,  // Input
                                    vector_type& output,       vector_type& outputA,        magnitude_type& theta) ;   // Output
        void get_eigenvalue(std::vector<double>& eigval, std::vector<class std::vector<double> >& eigvecs,
                            fortran_int_t dim, fortran_int_t i1, fortran_int_t i2);
        MATRIX const & matrix_;
        VS vecspace_;
        FortranMatrix<scalar_type> M;
        magnitude_type atol_;
        double omega_;
        bool shift_and_invert_ ;
        partial_overlap poverlap_ ;
        int site_ , n_mo_ , n_ ;
    };
    //
    // Methods of the Jacobi-Davidson class
    // ------------------------------------
    //
    // 1) constructor : just loads the vector space
    // 2) calculation of the eigenvectors (highest or lowest)
    // 3) get_extremal_eigenvalues : interfaces to FORTRAN routines to compute eigenvalues and eigenvectors
    //
    template <class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    jacobi_davidson<MATRIX, VS, OtherMatrix, SymmGroup>::jacobi_davidson(const MATRIX& matrix,
                                                                         const VS& vec,
                                                                         const partial_overlap& poverlap,
                                                                         const int& site,
                                                                         const int& n_mo):
        matrix_(matrix),
        vecspace_(vec),
        M(1,1),
        omega_(0.),
        poverlap_(poverlap),
        site_(site),
        n_mo_(n_mo),
        shift_and_invert_(false)
    {
        n_ = vec_dimension(vecspace_);
    }
    template <class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    jacobi_davidson<MATRIX, VS, OtherMatrix, SymmGroup>::jacobi_davidson(const MATRIX& matrix,
                                                                         const VS& vec,
                                                                         const double& omega,
                                                                         const partial_overlap& poverlap,
                                                                         const int& site,
                                                                         const int& n_mo) :
            matrix_(matrix),
            vecspace_(vec),
            M(1,1),
            omega_(omega),
            poverlap_(poverlap),
            site_(site),
            n_mo_(n_mo),
            shift_and_invert_(true)
    {
        n_ = vec_dimension(vecspace_);
    }
    template <class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    jacobi_davidson<MATRIX, VS, OtherMatrix, SymmGroup>::~jacobi_davidson() { }
    //
    // Method of the jacobi_davidson class to compute the eigenvalues
    // --------------------------------------------------------------
    template <class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    template <class GEN, class SOLVER, class ITER>
    std::pair<typename jacobi_davidson<MATRIX,VS,OtherMatrix,SymmGroup>::magnitude_type,
              typename jacobi_davidson<MATRIX,VS,OtherMatrix,SymmGroup>::vector_type>
    jacobi_davidson<MATRIX, VS, OtherMatrix, SymmGroup>::calculate_eigenvalue(const GEN& gen,
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
        bool converged ;
        //
        vector_type u, uA, eigvec;
        magnitude_type eigval ;
        //
        ietl::generate(V[0],gen);
        const_cast<GEN&>(gen).clear();
        ietl::project(V[0],vecspace_);
        //
        // Main loop of the algorithm
        // --------------------------
        do {
            update_vecspace(V , VA ,iter.iterations() ) ;
            // Update of the M matrix and compute the eigenvalues and the eigenvectors
            for(int i = 1; i <= iter.iterations()+1; i++)
                M(i-1, iter.iterations()) = ietl::dot(V[i-1], VA[iter.iterations()]);
            diagonalize_and_select(V, VA, iter.iterations()+1, u, uA, theta ) ;
            // Check convergence
            ++iter;
            vector_type r = compute_error(u, uA, theta);
            converged     = check_convergence(u, r, theta, iter, eigvec, eigval);
            if (converged)
                return std::make_pair(eigval, eigvec/ietl::two_norm(eigvec));
            rel_tol = 1. / pow(2.,double(iter.iterations()+1));
            solver(u, theta, r, V[iter.iterations()], rel_tol) ;
            V[iter.iterations()].data().iter_index = VA[iter.iterations()-1].data().iter_index ;
            storage::migrate(V[iter.iterations()], parallel::scheduler_balanced_iterative(V[iter.iterations()].data()));
        } while(true);
    }
    //
    // Compute the action of an operator
    template <class Matrix, class VS, class OtherMatrix, class SymmGroup>
    typename jacobi_davidson<Matrix, VS, OtherMatrix, SymmGroup>::vector_type
    jacobi_davidson<Matrix, VS, OtherMatrix, SymmGroup>::apply_operator(jacobi_davidson<Matrix, VS, OtherMatrix, SymmGroup>::vector_type const & x)
    {
        vector_type y, buf ;
        bool check = this->shift_and_invert_ ;
        if (check){
            vector_type buf ;
            ietl::mult(this->matrix_ , x , buf);
            y = this->omega_*x - buf;
        } else {
            ietl::mult(this->matrix_ , x , y);
        }
        return y;
    };
    // Update the vector space in JCD iteration
    template <class Matrix, class VS, class OtherMatrix, class SymmGroup>
    void jacobi_davidson<Matrix, VS, OtherMatrix, SymmGroup>::update_vecspace(jacobi_davidson<Matrix, VS, OtherMatrix, SymmGroup>::vector_space& V,
                                                                              jacobi_davidson<Matrix, VS, OtherMatrix, SymmGroup>::vector_space& VA ,
                                                                              const int idx )
    {
        bool check = this->shift_and_invert_ ;
        vector_type& t = V[idx];
        if (check){
            vector_type tA = apply_operator(t);
            for (int i = 1; i <= idx ; i++) {
                t -= ietl::dot(VA[i-1], tA) * V[i-1];
                tA -= ietl::dot(VA[i-1], tA) * VA[i-1];
            }
            //ietl::project(t,this->vecspace_) ;
            //ietl::project(tA,this->vecspace_) ;
            t /= ietl::two_norm(tA);
            VA[idx] = tA/ietl::two_norm(tA);
        }
        else {
            for (int i = 1; i <= idx; i++)
                t -= ietl::dot(V[i-1], t) * V[i-1];
            ietl::project(t,this->vecspace_) ;
            t /= ietl::two_norm(t) ;
            VA[idx] = apply_operator(t) ;
        }

    };
    // Compute the error vector
    template <class Matrix, class VS, class OtherMatrix, class SymmGroup>
    typename jacobi_davidson<Matrix, VS, OtherMatrix, SymmGroup>::vector_type
    jacobi_davidson<Matrix,VS,OtherMatrix,SymmGroup>::compute_error(const jacobi_davidson<Matrix,VS,OtherMatrix,SymmGroup>::vector_type &u,
                                                                    const jacobi_davidson<Matrix,VS,OtherMatrix,SymmGroup>::vector_type &uA,
                                                                    jacobi_davidson<Matrix,VS,OtherMatrix,SymmGroup>::magnitude_type theta)
    {
        vector_type r = uA ;
        if (this->shift_and_invert_)
            r -= u / theta;
        else
            r -= theta*u;
        return r ;
    }
    // Check if the JD iteration is arrived at convergence
    template <class Matrix, class VS, class OtherMatrix, class SymmGroup>
    template <class ITER>
    bool jacobi_davidson<Matrix, VS, OtherMatrix, SymmGroup>::check_convergence(const jacobi_davidson<Matrix,VS,OtherMatrix,SymmGroup>::vector_type &u,
                                                                                const jacobi_davidson<Matrix,VS,OtherMatrix,SymmGroup>::vector_type &r,
                                                                                const jacobi_davidson<Matrix,VS,OtherMatrix,SymmGroup>::magnitude_type theta,
                                                                                ITER& iter,
                                                                                jacobi_davidson<Matrix,VS,OtherMatrix,SymmGroup>::vector_type &eigvec,
                                                                                jacobi_davidson<Matrix,VS,OtherMatrix,SymmGroup>::magnitude_type &eigval)
    {
        // Compute the error vector
        bool converged ;
        if(iter.finished(ietl::two_norm(r),1./theta)) {
            if (this->shift_and_invert_) {
                eigvec = u;
                eigval = this->omega_ - 1. / theta;
            } else {
                eigvec = u ;
                eigval = theta ;
            }
            converged = true;
            return converged;
        } else {
            converged = false ;
            return converged ;
        }
    };
    //
    // Driver for diagonalization routine, include (if requested) also the selection of the "optimal" overlap
    // based on the MO criterion
    template<class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void jacobi_davidson<MATRIX, VS, OtherMatrix, SymmGroup>::diagonalize_and_select(const vector_space& MPSTns_input,
                                                                                     const vector_space& MPSTns_input_A,
                                                                                     const fortran_int_t& dim,
                                                                                     vector_type& MPSTns_output,
                                                                                     vector_type& MPSTns_output_A,
                                                                                     magnitude_type &theta)
    {
        // Initialization
        typedef typename std::vector<double> vector_scalars ;
        typedef typename jacobi_davidson::vector_type MPSTns_type ;
        double thresh = 0.50 ;
        vector_scalars eigvals , overlaps ;
        std::vector< vector_scalars > eigvecs ;
        MPSTns_type u_local , uA_local ;
        int imin , imax ;
        int  n_eigen  = ((n_mo_ > dim) ? dim : n_mo_) ;
        // Definition of the dimensions and dynamic memory allocation
        assert (n_eigen > 0) ;
        overlaps.resize(n_eigen) ;
        eigvals.resize(n_eigen) ;
        eigvecs.resize(n_eigen) ;
        for (int i = 0 ; i < n_eigen ; ++i)
            eigvecs[i].resize(dim) ;
        if (n_eigen == 1) {
            imin = imax = 1;
        } else {
            imin = 1 ;
            imax = n_eigen ;
        }
        // Diagonalization
        get_eigenvalue(eigvals, eigvecs, dim , imin, imax) ;
        // Finalization
        for (int i = 0 ; i < n_eigen ; ++i){
            // Conversion to the original basis
            u_local = eigvecs[i][0]*MPSTns_input_A[0] ;
            for(int j = 1; j < dim; ++j)
                u_local += eigvecs[i][j] * MPSTns_input_A[j];
            double scr = poverlap_.overlap(u_local, site_) ;
            overlaps[i] = fabs(scr) ;
        }
        int idx = 0 ;
        for (int i = 1 ; i < n_eigen ; ++i)
            if (overlaps[i] > overlaps[idx])
                idx = i ;
        //if ( overlaps[idx] < thresh )
        //    throw std::runtime_error("Satisfactory overlap not found");
        // Finalization
        MPSTns_output   = eigvecs[idx][0]*MPSTns_input[0] ;
        MPSTns_output_A = eigvecs[idx][0]*MPSTns_input_A[0] ;
        for (int j = 1; j < dim; ++j) {
            MPSTns_output   += eigvecs[idx][j]*MPSTns_input[j] ;
            MPSTns_output_A += eigvecs[idx][j]*MPSTns_input_A[j] ;
        }
        theta = eigvals[idx] ;
        // Print summary
        std::cout << " +---------------------------+" << std::endl ;
        std::cout << "  Maximum Overlap Calculation " << std::endl ;
        std::cout << " +---------------------------+" << std::endl ;
        std::cout << " Selected index - " << idx << std::endl ;
        std::cout << " Overlap value  - " << overlaps[idx] << std::endl ;
    };
    //
    // Interface to LAPACK diagonalization routine
    template <class MATRIX, class VS, class OtherMatrix, class SymmGroup>
    void jacobi_davidson<MATRIX, VS, OtherMatrix, SymmGroup>::get_eigenvalue(std::vector<double>& eigval,
                                                                             std::vector<class std::vector<double> >& eigvec,
                                                                             fortran_int_t dim,
                                                                             fortran_int_t id_min,
                                                                             fortran_int_t id_max)
    {
        // Definition of all the quantities needed by the LAPACK routine
        double abstol = atol_;
        char jobz='V';
        char range='I';
        char uplo='U';
        fortran_int_t n     = dim ;
        fortran_int_t lda   = dim ;
        fortran_int_t ldz   = n   ;
        fortran_int_t lwork = 8*n ;
        fortran_int_t info;
        fortran_int_t neig  = id_max-id_min+1 ;
        double vl, vu;
        double *w = new double[neig];
        double *z = new double[neig*n];
        double *work = new double[lwork];
        fortran_int_t *iwork = new fortran_int_t[5*n];
        fortran_int_t *ifail = new fortran_int_t[n];
        // Convert the matrix from general MATRIX class to a FortranMatrix object
        FortranMatrix<scalar_type> M_(dim,dim);
        for (int i=0 ; i<dim ; i++)
            for (int j=0 ; j<=i ; j++)
                M_(j,i) = M(j,i);
        LAPACK_DSYEVX(&jobz, &range, &uplo, &n, M_.data(), &lda, &vl, &vu, &id_min, &id_max, &abstol, &neig, w, z, &ldz, work, &lwork, iwork, ifail, &info);
        for (int j = 0 ; j < neig ; j++) {
            eigval[j] = w[j];
            for (int i = 0; i < n; i++)
                eigvec[j][i] = z[i + n*j];
        }
        // Free space
        delete [] w     ;
        delete [] z     ;
        delete [] work  ;
        delete [] iwork ;
        delete [] ifail ;
    };
}
#endif
