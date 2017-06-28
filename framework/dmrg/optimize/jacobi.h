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
    // +---------------------+
    //  JACOBI-DAVIDSON CLASS
    // +---------------------+
    // This is a general class for Davidson-type eigensolver.
    // The templates arguments are MATRIX and VS, that are usually:
    // MATRIX    : a SiteProblem object (see optimize.h for additional details)
    // VS        : a VectorSpace object, including a MPSTensor and several other vectors
    //             (for excited states orthogonalization)
    //
    // Includes the following attributes, that are common to all the Davidson eigensolvers;
    // matrix_   : the matrix representation of the operator for the site where the optimization is
    //             carried out
    // vecspace_ : the vector space where the optimization is carried out
    // atol_     : tolerance for assessing convergence of the optimization
    // Other attributes, that are specific of other type of eigensolvers (such as a shift omega or a
    // state to target) are defined in the inherited classes.
    //
    // The methods are the following:
    // contructor           : standard constructor
    // calculate_eigenvalue : method to compute eigenvalue, uses virtual functions
    // update_vspace        : virtual protected function, defined when virtual class is inherited
    // apply_operator       : virtual protected function, defined when virtual class is inherited
    // precondition         : virtual protected function, for the guess iteration of each Davidson step
    //
    template <class MATRIX, class VS, class ITER>
    class jacobi_davidson
    {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename std::vector<vector_type> vector_space;
        typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
        jacobi_davidson(const MATRIX& matrix, const VS& vec, const int& site);
        virtual ~jacobi_davidson() {};
        template <class GEN>
        std::pair<magnitude_type, vector_type> calculate_eigenvalue(const GEN& gen, ITER& iter);
    protected:
        // Private method, interface to the LAPACK diagonalization routine
        void get_eigenvalue(std::vector<double>& eigval, std::vector<class std::vector<double> >& eigvecs, fortran_int_t dim,
                            fortran_int_t i1, fortran_int_t i2) ;
        // Virtual protected methods, to be inherited by derived classes
        virtual bool check_convergence(const vector_type& u, const vector_type& uA , const magnitude_type theta ,
                                       ITER& iter, vector_type& eigvec, magnitude_type& eigval) {};
        virtual vector_type apply_operator (const vector_type& x) {} ;
        virtual void update_vecspace(vector_space &V, vector_space &VA, const int i) {};
        virtual vector_type compute_error (const vector_type& u , const vector_type& uA, magnitude_type theta) {} ;
        virtual void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                            vector_type& output, vector_type& outputA, magnitude_type& theta) {} ;
        virtual void solver(const vector_type& u, const magnitude_type& theta, const vector_type& r, vector_type& t,
                            const magnitude_type& rel_tol ) {} ;
        // Protected attributes
        MATRIX const & matrix_ ;
        VS vecspace_ ;
        int site_ ;
        FortranMatrix<scalar_type> M ;
        magnitude_type atol_ ;
        std::size_t max_iter_ ;
    };

    // -- Constructor --
    template <class MATRIX, class VS, class ITER>
    jacobi_davidson<MATRIX, VS, ITER>::jacobi_davidson(const MATRIX& matrix, const VS& vec, const int& site) :
        matrix_(matrix),
        vecspace_(vec),
        site_(site),
        M(1,1),
        max_iter_(0)
    { } ;
    // -- Calculation of eigenvalue --
    template <class MATRIX, class VS, class ITER>
    template <class GEN>
    std::pair<typename jacobi_davidson<MATRIX,VS,ITER>::magnitude_type,
              typename jacobi_davidson<MATRIX,VS,ITER>::vector_type>
    jacobi_davidson<MATRIX, VS, ITER>::calculate_eigenvalue(const GEN& gen, ITER& iter)
    {
        // Variable declaration
        // Scalars
        magnitude_type eigval, rel_tol, theta;
        atol_ = iter.absolute_tolerance();
        bool converged ;
        // Vectors
        std::vector<vector_type> V(iter.max_iterations());
        std::vector<vector_type> VA(iter.max_iterations());
        vector_type u, uA, eigvec;
        // Initialization
        ietl::generate(V[0],gen);
        const_cast<GEN&>(gen).clear();
        ietl::project(V[0],vecspace_);
        M.resize(iter.max_iterations(), iter.max_iterations());
        //
        // Main loop of the algorithm
        // --------------------------
        do {
            update_vecspace(V , VA ,iter.iterations() ) ;
            // Update of the M matrix and compute the eigenvalues and the eigenvectors
            for(int i = 1; i <= iter.iterations()+1; i++)
                M(i-1, iter.iterations()) = ietl::dot(V[i - 1], VA[iter.iterations()]);
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
    // -- Interface to the LAPACK diagonalization routine
    template <class MATRIX, class VS, class ITER>
    void jacobi_davidson<MATRIX, VS, ITER>::get_eigenvalue(std::vector<double>& eigval, std::vector<class std::vector<double> >& eigvec,
                                                           fortran_int_t dim, fortran_int_t id_min, fortran_int_t id_max)
    {
        // Definition of all the quantities needed by the LAPACK routine
        double abstol = atol_;
        char jobz  = 'V';
        char range = 'I';
        char uplo  = 'U';
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
                M_(j, i) = M(j, i);
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
