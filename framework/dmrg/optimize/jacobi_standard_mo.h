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

#ifndef IETL_JACOBI_STANDARD_MO_H
#define IETL_JACOBI_STANDARD_MO_H

#include <ietl/traits.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>
#include <ietl/cg.h>
#include <vector>

#include "dmrg/optimize/jacobi.h"
#include "dmrg/optimize/jacobi_standard.h"
#include "dmrg/optimize/partial_overlap.h"

// +-------------------------------------------------+
//  JACOBI-DAVIDSON EIGENSOLVER WITH OVERLAP TRACKING
// +-------------------------------------------------+
namespace ietl
{
    template <class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    class jacobi_davidson_standard_mo : public jacobi_davidson_standard<MATRIX, VS, ITER>
    {
    public:
        typedef jacobi_davidson_standard<MATRIX, VS, ITER> base;
        typedef typename partial_overlap<OtherMatrix,SymmGroup>::partial_overlap partial_overlap;
        typedef typename std::vector<partial_overlap> pov_vec_type ;
        typedef typename base::couple_vec      couple_vec ;
        typedef typename base::gt_couple       gt_couple ;
        typedef typename base::magnitude_type  magnitude_type;
        typedef typename base::matrix_double   matrix_double;
        typedef typename base::scalar_type     scalar_type ;
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
        using base::nsites_ ;
        using base::overlap_ ;
        using base::site1_ ;
        using base::site2_ ;
        using base::vecspace_ ;
        using base::v_guess_ ;
        //
        jacobi_davidson_standard_mo(const MATRIX& matrix, VS& vec, const pov_vec_type& pov, const size_t n,
                                    const int& side_tofollow, const int& nmin, const int& nmax, const int& max_iter,
                                    const int& nsites, const int& site1, const int& site2, const double& ietl_atol,
                                    const double& ietl_rtol, const int& i_gmres_guess, const std::vector<int>& order,
                                    const int& sa_alg, const int& root_homing_type)
                : base::jacobi_davidson_standard(matrix, vec, nmin, nmax, max_iter, nsites, site1, site2, ietl_atol,
                                                 ietl_rtol, i_gmres_guess, order, sa_alg)
                , pov_(pov) , n_maxov_(n), root_homing_type_(root_homing_type), side_tofollow_(side_tofollow) {} ;
        ~jacobi_davidson_standard_mo() {} ;
    private:
        double compute_overlap(const vector_type& vec_test) ;
        vector_double generate_property(const vector_space& V, const vector_space& VA, const size_t& dim,
                                        const matrix_double& eigvecs, const vector_double& eigvals);
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                    const int& mod, vector_type& output, vector_type& outputA, magnitude_type& theta,
                                    matrix_double& eigvecs, vector_double& eigvals) ;
        void diagonalize_first(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                               vector_type& output, vector_type& outputA, magnitude_type& theta,
                               matrix_double& eigvecs, vector_double& eigvals) ;
        void diagonalize_second(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                vector_type& output, vector_type& outputA, magnitude_type& theta,
                                matrix_double& eigvecs, vector_double& eigvals) ;
        void print_endline(void) ;
        void print_header_table(void) ;
        void print_newline_table(const size_t& i, const double& error, const magnitude_type& en, const double& overlap) ;
        void sort_prop(couple_vec& vector_values) ;
        // Private attributes
        int root_homing_type_ , side_tofollow_ ;
        pov_vec_type pov_  ;
        size_t n_maxov_  ;
    };
    // Diagonalization routine
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::diagonalize_and_select
            (const vector_space& MPSTns_input,
             const vector_space& MPSTns_input_A,
             const fortran_int_t& dim,
             const int& mod,
             vector_type& MPSTns_output,
             vector_type& MPSTns_output_A,
             magnitude_type &theta,
             matrix_double& eigvecs,
             vector_double& eigvals)
    {
        // Initialization
        int imin , imax , nevec;
        vector_double overlaps ;
        vector_type u_local , uA_local ;
        // Definition of the dimensions and dynamic memory allocation
        if (dim != n_restart_max_ && n_maxov_ > 0)
            nevec  = ((n_maxov_ > dim) ? dim : n_maxov_) ;
        else
            nevec = dim ;
        // Definition of the dimensions
        if (mod == 0) {
            imin = 1;
            imax = nevec;
        } else if (mod == 1) {
            imin = dim-nevec+1;
            imax = dim;
        } else {
            throw std::runtime_error("Unrecognized modality in diagonalize_and_select") ;
        }
        overlaps.resize(nevec) ;
        eigvals.resize(nevec) ;
        eigvecs.resize(nevec) ;
        for (int i = 0 ; i < nevec ; ++i)
            eigvecs[i].resize(dim) ;
        // Diagonalization
        get_eigenvalue(eigvals, eigvecs, dim , imin, imax) ;
        // Eigenvalue selection
        for (int i = 0; i < nevec; ++i) {
            // Conversion to the original basis
            u_local = eigvecs[i][0] * MPSTns_input[0];
            for (int j = 1; j < dim; ++j)
                u_local += eigvecs[i][j] * MPSTns_input[j];
            overlaps[i] = compute_overlap(u_local) ;
        }
        // Finalization
        int idx = -1  ;
        if (mod == 0)
            overlap_ = 0. ;
        for (int i = 0; i < nevec; ++i) {
            if (overlaps[i] > overlap_) {
                idx = i;
                overlap_ = overlaps[idx];
            }
        }
        // Finalization
        if (idx != -1) {
            MPSTns_output = eigvecs[idx][0] * MPSTns_input[0];
            MPSTns_output_A = eigvecs[idx][0] * MPSTns_input_A[0];
            for (int j = 1; j < dim; ++j) {
                MPSTns_output += eigvecs[idx][j] * MPSTns_input[j];
                MPSTns_output_A += eigvecs[idx][j] * MPSTns_input_A[j];
            }
            theta = eigvals[idx];
        }
    };
    // Diagonalization routine
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::diagonalize_first
            (const vector_space& MPSTns_input,
             const vector_space& MPSTns_input_A,
             const fortran_int_t& dim,
             vector_type& MPSTns_output,
             vector_type& MPSTns_output_A,
             magnitude_type &theta,
             matrix_double& eigvecs,
             vector_double& eigvals)
    {
        if (side_tofollow_ == 0 || side_tofollow_ == 1)
            diagonalize_and_select(MPSTns_input, MPSTns_input_A, dim, 0, MPSTns_output, MPSTns_output_A, theta, eigvecs, eigvals) ;
    };
    // Diagonalization routine - rediagonalization
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::diagonalize_second
            (const vector_space& MPSTns_input,
             const vector_space& MPSTns_input_A,
             const fortran_int_t& dim,
             vector_type& MPSTns_output,
             vector_type& MPSTns_output_A,
             magnitude_type &theta,
             matrix_double& eigvecs,
             vector_double& eigvals)
    {
        if (side_tofollow_ == 0)
            diagonalize_and_select(MPSTns_input, MPSTns_input_A, dim, 1, MPSTns_output, MPSTns_output_A, theta, eigvecs, eigvals) ;
        else if (side_tofollow_ == -1)
            diagonalize_and_select(MPSTns_input, MPSTns_input_A, dim, 0, MPSTns_output, MPSTns_output_A, theta, eigvecs, eigvals) ;
    };
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    typename jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::vector_double
             jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::generate_property
             (const vector_space &V, const vector_space& VA, const size_t& dim,
               const matrix_double &eigvecs, const vector_double &eigvals)
    {
        // Variable declaration
        vector_type tmp_V  ;
        vector_double p_tmp(dim) ;
        // Rotates the properties
        for (int i = 0; i < dim ; i++) {
            tmp_V = eigvecs[i][0] * V[0];
            for (int j = 1 ; j < dim ; j++)
                tmp_V  += eigvecs[i][j] * V[j];
            if (nsites_ == 1)
                p_tmp[i] = pov_[i_state_].overlap(tmp_V/ietl::two_norm(tmp_V), site1_) ;
            else if (nsites_ == 2)
                p_tmp[i] = pov_[i_state_].overlap(tmp_V/ietl::two_norm(tmp_V), site1_, site2_) ;
        }
        return p_tmp ;
    }
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    double jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::compute_overlap(const vector_type &vec_test)
    {
        double ret1, ret2, ret ;
        // Calculates the two overlaps
        if (root_homing_type_ == 1 || root_homing_type_ == 3) {
            if (nsites_ == 1)
                ret1 = pov_[i_state_].overlap(vec_test/ietl::two_norm(vec_test), site1_);
            else
                ret1 = pov_[i_state_].overlap(vec_test/ietl::two_norm(vec_test), site1_, site2_);
        }
        if (root_homing_type_ == 2 || root_homing_type_ == 3)
            ret2 = ietl::dot(vec_test, v_guess_[i_state_]) / (ietl::two_norm(vec_test)*ietl::two_norm(v_guess_[i_state_]));
        // Finalizes the calculation
        if (root_homing_type_ == 1)
            ret = fabs(ret1) ;
        else if (root_homing_type_ == 2)
            ret = fabs(ret2) ;
        else
            ret = (fabs(ret1) + fabs(ret2)) * 0.5 ;
        return ret ;
    }
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::sort_prop(couple_vec& vector_values)
    {
        std::sort(vector_values.begin(),vector_values.end(),gt_couple()) ;
    }
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::print_header_table(void) {
        print_endline() ;
        std::cout << " Iteration |    Error    |    Energy    |  Overlap  " << std::endl ;
        print_endline() ;
    } ;
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::print_endline(void) {
        std::cout << "-----------+-------------+--------------+-----------" << std::endl ;
    } ;
    //
    template<class MATRIX, class VS, class ITER, class OtherMatrix, class SymmGroup>
    void jacobi_davidson_standard_mo<MATRIX, VS, ITER, OtherMatrix, SymmGroup>::print_newline_table(const size_t& i, const double& error,
                                                                                                    const magnitude_type& en, const double& overlap )
    {
        char buf[100] ;
	    int a = i, n;
        n = sprintf(buf, "%5d      | %1.4E  | %6.5f  |  %1.4f", a, error, en, overlap);
        std::cout << buf << std::endl;
    }
}

#endif
