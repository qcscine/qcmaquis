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

#ifndef IETL_JACOBI_STANDARD_H
#define IETL_JACOBI_STANDARD_H

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
    template <class MATRIX, class VS, class ITER>
    class jacobi_davidson_standard : public jacobi_davidson<MATRIX, VS, ITER>
    {
    public:
        typedef jacobi_davidson<MATRIX, VS, ITER> base;
        typedef typename base::magnitude_type  magnitude_type;
        typedef typename base::matrix_double   matrix_double;
        typedef typename base::property_vector property_vector;
        typedef typename base::scalar_type     scalar_type;
        typedef typename base::size_t          size_t ;
        typedef typename base::vector_double   vector_double;
        typedef typename base::vector_space    vector_space;
        typedef typename base::vector_type     vector_type;
        using base::atol_ ;
        using base::get_eigenvalue ;
        using base::M ;
        using base::matrix_ ;
        using base::max_iter_ ;
        using base::n_restart_max_ ;
        using base::n_restart_min_ ;
        using base::site_ ;
        using base::vecspace_ ;
        //
        jacobi_davidson_standard(const MATRIX& matrix, const VS& vec, const int& site, const int& nmin, const int& nmax)
                : base::jacobi_davidson(matrix, vec, site, nmin, nmax) {} ;
        ~jacobi_davidson_standard() {} ;
    protected:
        vector_type apply_operator (const vector_type& x);
        void update_vecspace(vector_space &V, vector_space &VA, const int i, property_vector& props);
    private:
        bool check_convergence(const vector_type& u, const vector_type& uA , const magnitude_type theta ,
                               ITER& iter, vector_type& eigvec, magnitude_type& eigval);
        vector_type compute_error (const vector_type& u , const vector_type& uA, magnitude_type theta) ;
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                    vector_type& output, vector_type& outputA, magnitude_type& theta,
                                    matrix_double& eigvecs, vector_double& eigvals) ;
        void solver(const vector_type& u, const magnitude_type& theta, const vector_type& r, vector_type& t,
                    const magnitude_type& rel_tol) ;
        void restart_jd(vector_space& V, vector_space& VA, const property_vector& props, const matrix_double& eigvecs,
                        const vector_double& eigvals ) ;
    };
    //
    // Compute the action of an operator
    template <class Matrix, class VS, class ITER>
    typename jacobi_davidson_standard<Matrix, VS, ITER>::vector_type jacobi_davidson_standard<Matrix, VS, ITER>::apply_operator(vector_type const & x)
    {
        vector_type y ;
        ietl::mult(this->matrix_ , x , y);
        return y;
    };
    // Update the vector space in JCD iteration
    template <class Matrix, class VS, class ITER>
    void jacobi_davidson_standard<Matrix, VS, ITER>::update_vecspace(vector_space& V, vector_space& VA, const int idx , property_vector& props)
    {
        vector_type& t = V[idx] ;
        for (int i = 1; i <= idx; i++)
            t -= ietl::dot(V[i-1], t) * V[i-1];
        ietl::project(t,this->vecspace_) ;
        t /= ietl::two_norm(t) ;
        VA[idx] = apply_operator(t) ;
        props[idx] = ietl::dot(V[idx],VA[idx]) ;
    };
    // Compute the error vector
    template <class Matrix, class VS, class ITER>
    typename jacobi_davidson_standard<Matrix, VS, ITER>::vector_type jacobi_davidson_standard<Matrix,VS,ITER>::compute_error(const vector_type &u,
                                                                                                                             const vector_type &uA,
                                                                                                                             magnitude_type theta)
    {
        vector_type r = uA ;
        r -= theta*u;
        return r ;
    }
    // Check if the JD iteration is arrived at convergence
    template <class Matrix, class VS, class ITER>
    bool jacobi_davidson_standard<Matrix, VS, ITER>::check_convergence(const vector_type &u, const vector_type &r, const magnitude_type theta,
                                                                       ITER& iter, vector_type &eigvec, magnitude_type &eigval)
    {
        // Compute the error vector
        bool converged ;
        if(iter.finished(ietl::two_norm(r),1./theta)) {
            eigvec = u ;
            eigval = theta ;
            converged = true;
            return converged;
        } else {
            converged = false ;
            return converged ;
        }
    };
    //
    template<class MATRIX, class VS, class ITER>
    void jacobi_davidson_standard<MATRIX, VS, ITER>::diagonalize_and_select(const vector_space& MPSTns_input,
                                                                            const vector_space& MPSTns_input_A,
                                                                            const fortran_int_t& dim,
                                                                            vector_type& MPSTns_output,
                                                                            vector_type& MPSTns_output_A,
                                                                            magnitude_type &theta,
                                                                            matrix_double& eigvecs,
                                                                            vector_double& eigvals)
    {
        // Initialization
        int imin , imax , nevec ;
        // Definition of the dimensions and dynamic memory allocation
        if (dim != n_restart_max_) {
            imin = imax = 1;
            nevec = 1 ;
        } else {
            imin  = 1;
            imax  = dim ;
            nevec = imax - imin + 1 ;
        }
        eigvals.resize(dim) ;
        eigvecs.resize(nevec) ;
        for (int i = 0 ; i < nevec ; i++)
            eigvecs[i].resize(dim) ;
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
    void jacobi_davidson_standard<MATRIX, VS, ITER>::solver(const vector_type& u, const magnitude_type& theta, const vector_type& r,
                                                            vector_type& t, const magnitude_type& rel_tol)
    {
        jcd_solver_operator_standard<MATRIX, VS, vector_type> op(u, r, matrix_, theta);
        ietl_gmres gmres(max_iter_, false);
        vector_type inh = r;
        // initial guess for better convergence
        scalar_type dru = ietl::dot(r,u);
        scalar_type duu = ietl::dot(u,u);
        t = r + dru/duu*u;
        if (max_iter_ > 0)
            t = gmres(op, inh, t, rel_tol);
    } ;
    template<class MATRIX, class VS, class ITER>
    void jacobi_davidson_standard<MATRIX, VS, ITER>::restart_jd(vector_space &V, vector_space &VA, const property_vector& props,
                                                                const matrix_double &eigvecs, const vector_double &eigvals)
    {
        // Variable declaration
        typedef std::pair<int, float> couple_val ;
        std::vector<couple_val> vector_values ;
        vector_type V_tmp, VA_tmp ;
        size_t idx ;
        // Build the vector
        for (int i = 0; i < n_restart_max_ ; i++)
            vector_values.push_back(std::make_pair(i,eigvals[i]));
        std::sort(vector_values.begin(),vector_values.end(),[](couple_val a, couple_val b)
        {
            return b.second < b.first ;
        });
        // Finalization
        std::cout << "Pippo" << std::endl ;
        std::cout << n_restart_min_ << std::endl ;
        std::cout << n_restart_max_ << std::endl ;
        for (int i = 0; i < n_restart_min_; i++) {
            idx = vector_values[i].second ;
            std::cout << idx << std::endl ;
            V_tmp  = eigvecs[idx][0] * V[0];
            VA_tmp = eigvecs[idx][0] * VA[0];
            for (int j = 1; j < n_restart_max_; ++j) {
                V_tmp  += eigvecs[idx][j] * V[j];
                VA_tmp += eigvecs[idx][j] * VA[j];
            }
            V[i]  = V_tmp ;
            VA[i] = VA_tmp ;
        }
        V[n_restart_min_] = V[n_restart_max_+1] ;
    }
    //
    // JACOBI-DAVIDSON EIGENSOLVER WITH OVERLAP TRACKING
    // -------------------------------------------------
    template <class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    class jacobi_davidson_standard_mo : public jacobi_davidson_standard<MATRIX, VS, ITER>
    {
    public:
        typedef jacobi_davidson_standard<MATRIX, VS, ITER> base;
        typedef typename base::magnitude_type magnitude_type;
        typedef typename base::matrix_double  matrix_double;
        typedef typename partial_overlap<OtherMatrix,SymmGroup>::partial_overlap partial_overlap;
        typedef typename base::property_vector property_vector ;
        typedef typename base::scalar_type     scalar_type ;
        typedef typename base::size_t          size_t ;
        typedef typename base::vector_double   vector_double ;
        typedef typename base::vector_space    vector_space ;
        typedef typename base::vector_type     vector_type ;
        using base::apply_operator ;
        using base::get_eigenvalue ;
        using base::matrix_ ;
        using base::n_restart_max_ ;
        using base::n_restart_mint _ ;
        using base::site_ ;
        using base::vecspace_ ;
        //
         jacobi_davidson_standard_mo(const MATRIX& matrix, const VS& vec, const int& site, const partial_overlap& pov,
                                     const size_t n, const size_t& nmin=1, const size_t& nmax=20)
                : base::jacobi_davidson_standard(matrix, vec, site, nmin, nmax) , pov_(pov) , n_maxov_(n) {} ;
        ~jacobi_davidson_standard_mo() {} ;
    private:
        void update_vecspace(vector_space& V, vector_space& VA, const int idx , property_vector& props) ;
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                    vector_type& output, vector_type& outputA, magnitude_type& theta,
                                    matrix_double& eigvecs, vector_double& eigvals) ;
        void restart_jd(vector_space& V, vector_space& VA, const property_vector& props, const matrix_double& eigvecs,
                        const vector_double& eigvals ) ;
        // Private attributes
        partial_overlap pov_  ;
        std::size_t n_maxov_  ;
    };
    // Update the vector space in JCD iteration
    template <class Matrix, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<Matrix, VS, ITER, OtherMatrix, SymmGroup>::update_vecspace(vector_space& V, vector_space& VA, const int idx , property_vector& props)
    {
        vector_type& t = V[idx] ;
        for (int i = 1; i <= idx; i++)
            t -= ietl::dot(V[i-1], t) * V[i-1];
        ietl::project(t,this->vecspace_) ;
        t /= ietl::two_norm(t) ;
        VA[idx] = apply_operator(t) ;
        double scr = pov_.overlap(t/ietl::two_norm(t), site_);
        props[idx] = scr ;
    };
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::diagonalize_and_select
                    (const vector_space& MPSTns_input,
                     const vector_space& MPSTns_input_A,
                     const fortran_int_t& dim,
                     vector_type& MPSTns_output,
                     vector_type& MPSTns_output_A,
                     magnitude_type &theta,
                     matrix_double& eigvecs,
                     vector_double& eigvals)
    {
        // Initialization
        double thresh = 0.50 ;
        vector_double overlaps ;
        vector_type u_local , uA_local ;
        int imin , imax ;
        int n_eigen  = ((n_maxov_ > dim) ? dim : n_maxov_) ;
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
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::restart_jd
            (vector_space &V, vector_space &VA, const property_vector& props,
             const matrix_double &eigvecs, const vector_double &eigvals)
    {
        // Variable declaration
        typedef std::pair<int, float> couple_val ;
        std::vector<couple_val> vector_values ;
        vector_type V_tmp, VA_tmp ;
        size_t idx ;
        // Build the vector
        for (int i = 0; i < n_restart_max_ ; i++)
            vector_values.push_back(std::make_pair(i,eigvals[i]));
        std::sort(vector_values.begin(),vector_values.end(),[](couple_val a, couple_val b)
        {
            return b.second > b.first ;
        });
        // Finalization
        for (int i = 0; i < n_restart_min_-1; i++) {
            idx = vector_values[i].second() ;
            V_tmp  = eigvecs[idx][0] * MPSTns_input[0];
            VA_tmp = eigvecs[idx][0] * MPSTns_input_A[0];
            for (int j = 1; j < n_restart_max_; ++j) {
                V_tmp  += eigvecs[idx][j] * MPSTns_input[j];
                VA_tmp += eigvecs[idx][j] * MPSTns_input_A[j];
            }
            V[i]  = V_tmp ;
            VA[i] = VA_tmp ;
        }
        V[n_restart_min_] = V[n_restart_max_] ;
    }
}

#endif
