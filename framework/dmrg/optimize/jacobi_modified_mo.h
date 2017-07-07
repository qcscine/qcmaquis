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
#include <ietl/gmres.h>

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
        typedef typename base::couple_vec      couple_vec ;
        typedef typename base::gt_couple       gt_couple ;
        typedef typename base::magnitude_type  magnitude_type ;
        typedef typename base::matrix_double   matrix_double ;
        typedef typename partial_overlap<OtherMatrix,SymmGroup>::partial_overlap partial_overlap ;
        typedef typename base::size_t          size_t ;
        typedef typename base::vector_double   vector_double ;
        typedef typename base::vector_space    vector_space ;
        typedef typename base::vector_type     vector_type ;
        using base::apply_operator ;
        using base::get_eigenvalue ;
        using base::n_restart_max_ ;
        using base::omega_ ;
        using base::overlap_ ;
        using base::site_ ;
        using base::vecspace_ ;
        //
        jacobi_davidson_modified_mo(const MATRIX& matrix, const VS& vec, const int& site, const magnitude_type& omega,
                                    const partial_overlap& pov, const size_t n, const size_t& nmin=1, const size_t& nmax=20)
                : base::jacobi_davidson_modified(matrix, vec, site, omega, nmin, nmax) , pov_(pov) , n_maxov_(n) {} ;
        ~jacobi_davidson_modified_mo() {} ;
    private:
        vector_double generate_property(const vector_space& V, const vector_space& VA, const size_t& dim,
                                        const matrix_double& eigvecs, const vector_double& eigvals) ;
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                    vector_type& output, vector_type& outputA, magnitude_type& theta,
                                    matrix_double& eigvecs, vector_double& eigvals) ;
        void print_header_table(void) ;
        void print_endline(void) ;
        void print_newline_table(const size_t& i, const double& error, const magnitude_type& en, const double& overlap) ;
        void sort_prop(couple_vec& vector_values) ;
        void update_vecspace(vector_space& V, vector_space& VA, const int idx) ;
        // Private attributes
        partial_overlap pov_  ;
        size_t n_maxov_  ;
    };
    // Update the vector space in JCD iteration
    template <class Matrix, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_modified_mo<Matrix, VS, ITER, OtherMatrix, SymmGroup>::update_vecspace
            (vector_space& V, vector_space& VA, const int idx)
    {
        vector_type& t  = V[idx] ;
        vector_type  tA = apply_operator(t);
        for (int i = 1; i <= idx ; i++) {
            t -= ietl::dot(VA[i-1], tA) * V[i-1];
            tA -= ietl::dot(VA[i-1], tA) * VA[i-1];
        }
        ietl::project(t,this->vecspace_) ;
        ietl::project(tA,this->vecspace_) ;
        t /= ietl::two_norm(tA) ;
        VA[idx]    = tA/ietl::two_norm(tA) ;
    };
    //
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
        double thresh = 0.50 ;
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
        eigvals.resize(nevec) ;
        eigvecs.resize(nevec) ;
        for (int i = 0 ; i < nevec ; ++i)
            eigvecs[i].resize(dim) ;
        // Diagonalization
        get_eigenvalue(eigvals, eigvecs, dim , imin, imax) ;
        int idx = 0;
        for (int i = 0; i < nevec; ++i) {
            // Conversion to the original basis
            u_local = eigvecs[i][0] * MPSTns_input[0];
            for (int j = 1; j < dim; ++j)
                u_local += eigvecs[i][j] * MPSTns_input[j];
            double scr = pov_.overlap(u_local, site_);
            overlaps[i] = fabs(scr);
        }
        for (int i = 1; i < nevec; ++i) {
            if (overlaps[i] > overlaps[idx])
                idx = i;
        }
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
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    typename jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::vector_double
    jacobi_davidson_modified_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::generate_property
            (const vector_space &V, const vector_space& VA, const size_t& dim,
             const matrix_double &eigvecs, const vector_double &eigvals)
    {
        // Variable declaration
        vector_type   tmp_V ;
        vector_double p_tmp(dim) ;
        // Rotates the properties
        for (int i = 0; i < dim ; i++) {
            tmp_V = eigvecs[i][0] * V[0];
            for (int j = 1 ; j < dim ; j++)
                tmp_V  += eigvecs[i][j] * V[j];
            p_tmp[i] = pov_.overlap(tmp_V/ietl::two_norm(tmp_V), site_) ;
        }
        return p_tmp ;
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
        char buf[100];
	int a = i , n ;
        n = sprintf(buf, "%5d      | %1.4E  | %6.5f  |  %1.4f", a, error, en, overlap);
        std::cout << buf << std::endl;
    }
}
#endif
