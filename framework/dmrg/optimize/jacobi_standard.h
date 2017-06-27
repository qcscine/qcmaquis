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
        typedef typename base::vector_type  vector_type;
        typedef typename base::scalar_type  scalar_type;
        typedef typename base::vector_space vector_space;
        typedef typename base::magnitude_type magnitude_type;
        using base::matrix_ ;
        using base::vecspace_ ;
        using base::M ;
        using base::atol_ ;
        using base::max_iter_ ;
        using base::get_eigenvalue ;
        //
        jacobi_davidson_standard(const MATRIX& matrix, const VS& vec) : base::jacobi_davidson(matrix, vec) {} ;
        ~jacobi_davidson_standard() {} ;
    private:
        // Multiply the vector by a matrix
        vector_type apply_operator (const vector_type& x);
        void update_vecspace(vector_space &V, vector_space &VA, const int i);
        bool check_convergence(const vector_type& u, const vector_type& uA , const magnitude_type theta ,
                               ITER& iter, vector_type& eigvec, magnitude_type& eigval);
        vector_type compute_error (const vector_type& u , const vector_type& uA, magnitude_type theta) ;
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,  // Input
                                    vector_type& output, vector_type& outputA, magnitude_type& theta) ;   // Output
        void solver(const vector_type& u, const magnitude_type& theta, const vector_type& r, vector_type& t,
                    const magnitude_type& rel_tol) ;
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
    void jacobi_davidson_standard<Matrix, VS, ITER>::update_vecspace(vector_space& V, vector_space& VA, const int idx )
    {
        vector_type t = V[idx] ;
        for (int i = 1; i <= idx; i++)
            t -= ietl::dot(V[i-1], t) * V[i-1];
        ietl::project(t,this->vecspace_) ;
        t /= ietl::two_norm(t) ;
        VA[idx] = apply_operator(t) ;
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
                                                                            magnitude_type &theta)
    {
        // Initialization
        typedef typename std::vector<double> vector_scalars ;
        vector_scalars eigvals , overlaps ;
        std::vector< vector_scalars > eigvecs ;
        vector_type u_local , uA_local ;
        int imin , imax ;
        // Definition of the dimensions and dynamic memory allocation
        eigvals.resize(1) ;
        eigvecs.resize(dim) ;
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
    }
}

#endif
