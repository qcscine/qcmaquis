/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef IETL_JACOBI_MODIFIED_H
#define IETL_JACOBI_MODIFIED_H

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>
#include <ietl/cg.h>
#include <ietl/gmres.h>
#include <complex>
#include <vector>

#include <boost/function.hpp>
#include "dmrg/optimize/jacobi.h"
#include "dmrg/optimize/jcd_solver.h"
#include "dmrg/optimize/partial_overlap.h"

namespace ietl
{
    //
    // MODIFIED JACOBI-DAVIDSON EIGENSOLVER
    // ------------------------------------
    template <class MATRIX, class VS, class ITER>
    class jacobi_davidson_modified : public jacobi_davidson<MATRIX, VS, ITER>
    {
    public:
        typedef jacobi_davidson<MATRIX, VS, ITER> base;
        typedef typename base::vector_type  vector_type;
        typedef typename base::scalar_type  scalar_type;
        typedef typename base::vector_space vector_space;
        typedef typename base::magnitude_type magnitude_type;
        using base::matrix_ ;
        using base::vecspace_ ;
        using base::M ;
        using base::max_iter_ ;
        using base::get_eigenvalue ;
        //
        jacobi_davidson_modified(const MATRIX& matrix, const VS& vec, const int& site, const magnitude_type& omega)
                : base::jacobi_davidson(matrix, vec, site) , omega_(omega) {} ;
        ~jacobi_davidson_modified() {} ;
    private:
        // Multiply the vector by a matrix
        vector_type apply_operator (const vector_type& x);
        void update_vecspace(vector_space &V, vector_space &VA, const int i);
        bool check_convergence(const vector_type& u, const vector_type& uA , const magnitude_type theta ,
                               ITER& iter, vector_type& eigvec, magnitude_type& eigval);
        vector_type compute_error (const vector_type& u , const vector_type& uA, magnitude_type theta) ;
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                    vector_type& output, vector_type& outputA, magnitude_type& theta) ;
        void solver(const vector_type& u, const magnitude_type& theta, const vector_type& r, vector_type& t,
                    const magnitude_type& rel_tol) ;
    protected:
        // Private attributes
        magnitude_type omega_ ;
    };
    // Compute the action of an operator
    template <class Matrix, class VS, class ITER>
    typename jacobi_davidson_modified<Matrix, VS, ITER>::vector_type jacobi_davidson_modified<Matrix, VS, ITER>::apply_operator(vector_type const & x)
    {
        vector_type y, buf ;
        ietl::mult(this->matrix_ , x , buf);
        y = this->omega_*x - buf;
        return y;
    };
    // Update the vector space in JCD iteration
    template <class Matrix, class VS, class ITER>
    void jacobi_davidson_modified<Matrix, VS, ITER>::update_vecspace(vector_space& V, vector_space& VA, const int idx )
    {
        vector_type& t  = V[idx] ;
        vector_type  tA = apply_operator(t);
        for (int i = 1; i <= idx ; i++) {
            t -= ietl::dot(VA[i-1], tA) * V[i-1];
            tA -= ietl::dot(VA[i-1], tA) * VA[i-1];
        }
        //ietl::project(t,this->vecspace_) ;
        //ietl::project(tA,this->vecspace_) ;
        t /= ietl::two_norm(tA);
        VA[idx] = tA/ietl::two_norm(tA);
    };
    // Compute the error vector
    template <class Matrix, class VS, class ITER>
    typename jacobi_davidson_modified<Matrix, VS, ITER>::vector_type jacobi_davidson_modified<Matrix,VS,ITER>::compute_error(const vector_type &u,
                                                                                                                            const vector_type &uA,
                                                                                                                            magnitude_type theta)
    {
        vector_type r = uA ;
        r -= u / theta;
        return r ;
    }
    // Check if the JD iteration is arrived at convergence
    template <class Matrix, class VS, class ITER>
    bool jacobi_davidson_modified<Matrix, VS, ITER>::check_convergence(const vector_type &u, const vector_type &r, const magnitude_type theta,
                                                                       ITER& iter, vector_type &eigvec, magnitude_type &eigval)
    {
        // Compute the error vector
        bool converged ;
        if(iter.finished(ietl::two_norm(r),1./theta)) {
            eigvec = u;
            eigval = this->omega_ - 1. / theta;
            converged = true;
            return converged;
        } else {
            converged = false ;
            return converged ;
        }
    };
    //
    template<class MATRIX, class VS, class ITER>
    void jacobi_davidson_modified<MATRIX, VS, ITER>::diagonalize_and_select(const vector_space& MPSTns_input,
                                                                            const vector_space& MPSTns_input_A,
                                                                            const fortran_int_t& dim,
                                                                            vector_type& MPSTns_output,
                                                                            vector_type& MPSTns_output_A,
                                                                            magnitude_type &theta)
    {
        // Initialization
        typedef typename std::vector<double> vector_scalars ;
        vector_scalars eigvals ;
        std::vector< vector_scalars > eigvecs ;
        int imin , imax ;
        // Definition of the dimensions and dynamic memory allocation
        eigvals.resize(dim) ;
        eigvecs.resize(1) ;
        eigvecs[0].resize(dim) ;
        imin = imax = 1;
        // Diagonalization
        get_eigenvalue(eigvals, eigvecs, dim , imin, imax) ;
        int idx = 0;
        // Finalization
        MPSTns_output   = eigvecs[idx][0]*MPSTns_input[0] ;
        MPSTns_output_A = eigvecs[idx][0]*MPSTns_input_A[0] ;
        for (int j = 1; j < dim; ++j) {
            MPSTns_output   += eigvecs[idx][j]*MPSTns_input[j] ;
            MPSTns_output_A += eigvecs[idx][j]*MPSTns_input_A[j] ;
        }
        theta = eigvals[0] ;
    };
    //
    template<class MATRIX, class VS, class ITER>
    void jacobi_davidson_modified<MATRIX, VS, ITER>::solver(const vector_type& u, const magnitude_type& theta, const vector_type& r,
                                                            vector_type& t, const magnitude_type& rel_tol)
    {
        vector_type z, inh = r;
        z = apply_operator(u) ;
        jcd_solver_operator_modified<MATRIX, VS, vector_type> op(u, r, matrix_, theta, omega_, z);
        ietl_gmres gmres(max_iter_, false);
        // initial guess for better convergence
        scalar_type dru = ietl::dot(r,u);
        scalar_type duu = ietl::dot(u,u);
        t = r + dru/duu*u;
        if (max_iter_ > 0)
            t = gmres(op, inh, t, rel_tol);
    }
    //
    // MODIFIED JACOBI-DAVIDSON EIGENSOLVER WITH OVERLAP TRACKING
    // ---------------------------------------------------------
    template <class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    class jacobi_davidson_modified_mo : public jacobi_davidson_modified<MATRIX, VS, ITER>
    {
    public:
        typedef jacobi_davidson_modified<MATRIX, VS, ITER> base;
        typedef typename base::vector_type  vector_type;
        typedef typename base::scalar_type  scalar_type;
        typedef typename base::vector_space vector_space;
        typedef typename base::magnitude_type magnitude_type;
        typedef typename partial_overlap<OtherMatrix,SymmGroup>::partial_overlap partial_overlap;
        using base::matrix_ ;
        using base::vecspace_ ;
        using base::get_eigenvalue ;
        using base::site_ ;
        using base::omega_ ;
        //
        jacobi_davidson_modified_mo(const MATRIX& matrix, const VS& vec, const int& site, const magnitude_type& omega,
                                    const partial_overlap& pov, const std::size_t n)
                : base::jacobi_davidson_modified(matrix, vec, site, omega) , pov_(pov) , n_maxov_(n) {} ;
        ~jacobi_davidson_modified_mo() {} ;
    private:
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                    vector_type& output, vector_type& outputA, magnitude_type& theta) ;
        // Private attributes
        partial_overlap pov_  ;
        std::size_t n_maxov_  ;
    };
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::diagonalize_and_select
            (const vector_space& MPSTns_input,
             const vector_space& MPSTns_input_A,
             const fortran_int_t& dim,
             vector_type& MPSTns_output,
             vector_type& MPSTns_output_A,
             magnitude_type &theta)
    {
        // Initialization
        typedef typename std::vector<double> vector_scalars ;
        double thresh = 0.50 ;
        vector_scalars eigvals , overlaps ;
        std::vector< vector_scalars > eigvecs ;
        vector_type u_local , uA_local ;
        int imin , imax ;
        int  n_eigen  = ((n_maxov_ > dim) ? dim : n_maxov_) ;
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
        int idx = 0;
        for (int i = 0; i < n_eigen; ++i) {
            // Conversion to the original basis
            u_local = eigvecs[i][0] * MPSTns_input_A[0];
            for (int j = 1; j < dim; ++j)
                u_local += eigvecs[i][j] * MPSTns_input_A[j];
            double scr = pov_.overlap(u_local, site_);
            overlaps[i] = fabs(scr);
        }
        for (int i = 1; i < n_eigen; ++i)
            if (overlaps[i] > overlaps[idx])
                idx = i;
        // Finalization
        MPSTns_output   = eigvecs[idx][0]*MPSTns_input[0] ;
        MPSTns_output_A = eigvecs[idx][0]*MPSTns_input_A[0] ;
        for (int j = 1; j < dim; ++j) {
            MPSTns_output   += eigvecs[idx][j]*MPSTns_input[j] ;
            MPSTns_output_A += eigvecs[idx][j]*MPSTns_input_A[j] ;
        }
        theta = eigvals[idx] ;
        // Print summary
        std::cout << " +---------------------------+ " << std::endl;
        std::cout << "  Maximum Overlap Calculation  " << std::endl;
        std::cout << " +---------------------------+ " << std::endl;
        std::cout << " Selected index - " << idx << std::endl;
        std::cout << " Overlap value  - " << overlaps[idx] << std::endl;

    };
}
#endif
