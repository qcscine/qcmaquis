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

#ifndef IETL_JACOBI_MODIFIED_MO_H
#define IETL_JACOBI_MODIFIED_MO_H

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h> 
#include <ietl/cg.h>
#include <vector>

#include "dmrg/optimize/partial_overlap.h"

// +----------------------------------------------------------+
//  HARMONIC JACOBI-DAVIDSON EIGENSOLVER WITH OVERLAP TRACKING
// +----------------------------------------------------------+
namespace ietl
{
    template <class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    class jacobi_davidson_modified_mo : public jacobi_davidson_modified<MATRIX, VS, ITER>
    {
    public:
        typedef jacobi_davidson_modified<MATRIX, VS, ITER> base;
        typedef typename partial_overlap<OtherMatrix,SymmGroup>::partial_overlap partial_overlap;
        typedef typename std::vector<partial_overlap> pov_vec_type ;
        typedef typename base::couple_vec      couple_vec ;
        typedef typename base::gt_couple       gt_couple ;
        typedef typename base::magnitude_type  magnitude_type ;
        typedef typename base::matrix_double   matrix_double ;
        typedef typename base::size_t          size_t ;
        typedef typename base::vector_double   vector_double ;
        typedef typename base::vector_space    vector_space ;
        typedef typename base::vector_type     vector_type ;
        //
        using base::apply_operator ;
        using base::get_eigenvalue ;
        using base::i_gmres_guess_ ;
        using base::i_state_ ;
        using base::n_restart_max_ ;
        using base::n_root_found_ ;
        using base::nsites_ ;
        using base::omega_ ;
        using base::overlap_ ;
        using base::site1_ ;
        using base::site2_ ;
        using base::vecspace_ ;
        using base::v_guess_ ;
        //
        jacobi_davidson_modified_mo(const MATRIX& matrix, VS& vec, const magnitude_type& omega, const pov_vec_type& pov, const size_t n,
                                    const int& nmin, const int& nmax, const int& max_iter, const int& nsites, 
                                    const int& site1, const int& site2, const double& ietl_atol, const double& ietl_rtol,
                                    const int& i_gmres_guess, const std::vector<int>& order, const double& atol_init, const double& rtol_init,
                                    const size_t& max_iter_init, const int& root_homing_type)
                : base::jacobi_davidson_modified(matrix, vec, omega, nmin, nmax, max_iter, nsites, site1, site2, ietl_atol, ietl_rtol,
                                                 i_gmres_guess, order, atol_init, rtol_init, max_iter_init)
                , pov_(pov) , n_maxov_(n), root_homing_type_(root_homing_type) {} ;
        ~jacobi_davidson_modified_mo() {} ;
    private:
        bool check_convergence(const vector_type& u, const vector_type& uA, const vector_type& r, const magnitude_type theta ,
                               ITER& iter, vector_type& eigvec, magnitude_type& eigval);
        double compute_overlap(const vector_type& vec_test) ;
        vector_double generate_property(const vector_space& V, const vector_space& VA, const size_t& dim,
                                        const matrix_double& eigvecs, const vector_double& eigvals);
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                    vector_type& output, vector_type& outputA, magnitude_type& theta,
                                    matrix_double& eigvecs, vector_double& eigvals) ;
        void print_endline(void) ;
        void print_header_table(void) ;
        void print_newline_table(const size_t& i, const double& error, const magnitude_type& en, const double& overlap) ;
        void sort_prop(couple_vec& vector_values) ;
        // Private attributes
        double atol_init_ ;
        int root_homing_type_ ;
        pov_vec_type pov_  ;
        size_t n_maxov_  ;
    };
    // Diagonalization routine
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::diagonalize_and_select
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
        vector_double overlaps ;
        vector_type u_local , uA_local ;
        int imin , imax , nevec;
        // Definition of the dimensions and dynamic memory allocation
        if (dim != n_restart_max_ && n_maxov_ > 0) {
            nevec  = ((n_maxov_ > dim) ? dim : n_maxov_) ;
            imin   = 1      ;
            imax   = nevec ;
        } else {
            imin  = 1;
            imax  = dim ;
            nevec = imax - imin + 1 ;
        }
        // Definition of the dimensions and dynamic memory allocation
        assert (nevec > 0) ;
        overlaps.resize(nevec) ;
        eigvals.resize(dim) ;
        eigvecs.resize(nevec) ;
        for (int i = 0 ; i < nevec ; i++)
            eigvecs[i].resize(dim) ;
        // Diagonalization
        get_eigenvalue(eigvals, eigvecs, dim , imin, imax) ;
        int idx = 0;
        double scr ;
        for (int i = 0; i < nevec; ++i) {
            // Conversion to the original basis
            u_local = eigvecs[i][0] * MPSTns_input[0];
            for (int j = 1; j < dim; ++j)
                u_local += eigvecs[i][j] * MPSTns_input[j];
            overlaps[i] = compute_overlap(u_local) ;
        }
        for (int i = 1; i < nevec; ++i)
            if (overlaps[i] > overlaps[idx])
                idx = i;
        overlap_ = overlaps[idx] ;
        // Finalization
        MPSTns_output   = eigvecs[idx][0]*MPSTns_input[0] ;
        MPSTns_output_A = eigvecs[idx][0]*MPSTns_input_A[0] ;
        for (int j = 1; j < dim; ++j) {
            MPSTns_output   += eigvecs[idx][j]*MPSTns_input[j] ;
            MPSTns_output_A += eigvecs[idx][j]*MPSTns_input_A[j] ;
        }
        theta = eigvals[idx] ;
    };
    // Check if the JD iteration is arrived at convergence
    template <class Matrix, class VS, class ITER, class OtherMatrix, class SymmGroup>
    bool jacobi_davidson_modified_mo<Matrix, VS, ITER, OtherMatrix, SymmGroup>::check_convergence(const vector_type &u, const vector_type &uA, const vector_type &r,
                                                                                                  const magnitude_type theta, ITER& iter, vector_type &eigvec, 
                                                                                                  magnitude_type &eigval)
    {
        // Compute the error vector
        bool converged ;
        eigvec = u/ietl::two_norm(u);
        eigval = this->omega_ - theta/ietl::dot(u,u) ;
        if(iter.finished(ietl::two_norm(r),1.0) && overlap_ > 0.2) {
            converged = true;
            return converged;
        } else {
            converged = false ;
            return converged ;
        }
    };
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    typename jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::vector_double
    	     jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::generate_property
             (const vector_space &V, const vector_space& VA, const size_t& dim,
              const matrix_double &eigvecs, const vector_double &eigvals)
    {
        // Variable declaration
        vector_type  tmp_V ;
        vector_double p_tmp(dim) ;
        // Rotates the properties
        for (int i = 0; i < dim ; i++) {
            tmp_V = eigvecs[i][0] * V[0];
            for (int j = 1 ; j < dim ; j++)
                tmp_V  += eigvecs[i][j] * V[j];
            p_tmp[i] = compute_overlap(tmp_V) ;
        }
        return p_tmp ;
    }
    // Routine to compute the overlaps
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    double jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::compute_overlap(const vector_type &vec_test)
    {
        double ret, scr ;
        if (root_homing_type_ == 1) {
            if (nsites_ == 1)
                ret = pov_[i_state_].overlap(vec_test/ietl::two_norm(vec_test), site1_);
            else
                ret = pov_[i_state_].overlap(vec_test/ietl::two_norm(vec_test), site1_, site2_);
        } else {
            ret = ietl::dot(vec_test, v_guess_[i_state_]) / (ietl::two_norm(vec_test)*ietl::two_norm(v_guess_[i_state_]));
        }
        return fabs(ret) ;
    }
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::sort_prop(couple_vec &vector_values)
    {
        std::sort(vector_values.begin(),vector_values.end(),gt_couple()) ;
    }
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::print_header_table(void) {
        print_endline() ;
        std::cout << " Iteration |    Error    |    Energy    |  Overlap  " << std::endl ;
        print_endline() ;
    } ;
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::print_endline(void) {
        std::cout << "-----------+-------------+--------------+-----------" << std::endl ;
    } ;
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::print_newline_table(const size_t& i, const double& error,
                                                                                                    const magnitude_type& en, const double& overlap )
    {
        char buf[100] ;
	    int a = i, n;
        n = sprintf(buf, "%5d      | %1.4E  | %6.5f  |  %1.4f", a, error, en, overlap);
        std::cout << buf << std::endl;
    }
}
#endif
