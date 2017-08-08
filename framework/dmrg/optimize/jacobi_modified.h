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

#include <ietl/cg.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>
#include <ietl/traits.h>
#include <vector>
#include <cmath>

#include "dmrg/optimize/gmres_alb.h"
#include "dmrg/optimize/jacobi.h"
#include "dmrg/optimize/partial_overlap.h"

// +------------------------------------+
//  HARMONIC JACOBI-DAVIDSON EIGENSOLVER
// +------------------------------------+

namespace ietl
{
    template <class MATRIX, class VS, class ITER>
    class jacobi_davidson_modified : public jacobi_davidson<MATRIX, VS, ITER>
    {
    public:
        typedef jacobi_davidson<MATRIX, VS, ITER> base;
        typedef typename base::couple_vec         couple_vec;
        typedef typename base::lt_couple          lt_couple;
        typedef typename base::magnitude_type     magnitude_type;
        typedef typename base::matrix_double      matrix_double;
        typedef typename base::result_collector   result_collector ;
        typedef typename base::scalar_type        scalar_type;
        typedef typename base::size_t             size_t ;
        typedef typename base::vector_double      vector_double;
        typedef typename base::vector_ortho_vec   vector_ortho_vec;
        typedef typename base::vector_pairs       vector_pairs;
        typedef typename base::vector_space       vector_space;
        typedef typename base::vector_type        vector_type;
        //
        using base::get_eigenvalue ;
        using base::i_gmres_guess_ ;
        using base::i_state_ ;
        using base::ietl_atol_ ;
        using base::ietl_rtol_ ;
        using base::M ;
        using base::matrix_ ;
        using base::max_iter_ ;
        using base::n_restart_max_ ;
        using base::n_restart_min_ ;
        using base::n_root_found_ ;
        using base::n_sa_ ;
        using base::order_ ;
        using base::ortho_space_left_ ;
        using base::ortho_space_right_ ;
        using base::site1_ ;
        using base::u_and_uA_ ;
        using base::vecspace_ ;
        using base::v_guess_ ;
        //
        jacobi_davidson_modified(const MATRIX& matrix, VS& vec, const magnitude_type& omega, const size_t& nmin,
                                 const size_t& nmax, const size_t& max_iter, const int& nsites, const int& site1,
                                 const int& site2, const double& ietl_atol, const double& ietl_rtol, const size_t& i_gmres_guess,
                                 const std::vector<int>& order, const double& atol_init)
                : base::jacobi_davidson(matrix, vec, nmin, nmax, max_iter, nsites, site1, site2, ietl_atol, ietl_rtol, i_gmres_guess, order)
                , omega_(omega), atol_init_(atol_init) 
        {
            vector_type a, b ;
            for (size_t idx = 0; idx < n_sa_; idx++) {
                vector_type a = v_guess_[idx] ;
                vector_type b = 0.*v_guess_[idx] ;
                gmres_initializer_modified<MATRIX, vector_type, VS> gmres(this->matrix_, v_guess_[idx], vecspace_, omega_,
                                                                          ortho_space_left_, ortho_space_right_, idx, 100,
                                                                          false);
                v_guess_[idx] = gmres(a, b, 1.0E-5, 1.0E-3);
                v_guess_[idx] /= ietl::two_norm(v_guess_[idx]) ;
            }
        } ;
        ~jacobi_davidson_modified() {} ;
    private:
        // Methods
        bool check_convergence(const vector_type& u, const vector_type& uA, const vector_type& r, const magnitude_type theta ,
                               ITER& iter, vector_type& eigvec, magnitude_type& eigval);
        magnitude_type get_matrix_element(const vector_type& V, const vector_type& VA);
        vector_double generate_property(const vector_space& V, const vector_space& VA, const size_t& dim,
                                        const matrix_double& eigvecs, const vector_double& eigvals) ;
        vector_type compute_error (const vector_type& u , const vector_type& uA, magnitude_type theta) ;
        void diagonalize_and_select(const vector_space& input, const vector_space& inputA,  const fortran_int_t& dim,
                                    vector_type& output, vector_type& outputA, magnitude_type& theta,
                                    matrix_double& eigvecs, vector_double& eigvals) ;
        void print_endline(void) ;
        void print_header_table(void) ;
        void print_newline_table(const size_t& i, const double& error, const magnitude_type& en, const double& overlap) ;
        void solver(const vector_type& u, const vector_type& uA, const magnitude_type& theta, const vector_type& r, vector_type& t) ;
    protected:
        // Methods
        void sort_prop(couple_vec& vector_values) ;
        void update_vecspace(vector_space &V, vector_space &VA, const int i, vector_pairs& res);
        void update_orthospace(void) ;
        void update_u_and_uA(const vector_type& u, const vector_type& uA) ;
        // Attributes
        double atol_init_ ;
        magnitude_type omega_ ;
        vector_type apply_operator (const vector_type& x);
    };
    // Compute the action of an operator
    template <class Matrix, class VS, class ITER>
    typename jacobi_davidson_modified<Matrix, VS, ITER>::vector_type jacobi_davidson_modified<Matrix, VS, ITER>::apply_operator(vector_type const & x)
    {
        vector_type y, buf ;
        ietl::mult(this->matrix_ , x , buf, i_state_);
        y = this->omega_*x - buf;
        return y;
    };
    // Routine doing deflation
    template <class Matrix, class VS, class ITER>
    void jacobi_davidson_modified<Matrix,VS,ITER>::update_orthospace(void)
    {
        for (size_t jcont = 0; jcont < n_root_found_; jcont++) {
            vector_type tmp  = vecspace_.return_orthovec(u_and_uA_[jcont].first.first,  order_[n_root_found_], order_[jcont], site1_) ;
            vector_type tmpA = vecspace_.return_orthovec(u_and_uA_[jcont].first.second, order_[n_root_found_], order_[jcont], site1_) ;
            for (size_t j = 0 ; j < jcont ; j++) {
                tmp  -= ietl::dot(ortho_space_right_[j].first, tmpA) * ortho_space_right_[j].first  ;
                tmpA -= ietl::dot(ortho_space_right_[j].first, tmpA) * ortho_space_right_[j].second ;
            }
            if (ietl::two_norm(tmpA) > 0.) {
                tmp /= ietl::two_norm(tmpA);
                tmpA /= ietl::two_norm(tmpA);
            }
            ortho_space_left_.push_back(std::make_pair(tmpA, tmpA));
            ortho_space_right_.push_back(std::make_pair(tmpA, tmp));
        }
    }
    // Update the vector with the quantity to orthogonalize
    template <class Matrix, class VS, class ITER>
    void jacobi_davidson_modified<Matrix, VS, ITER>::update_u_and_uA(const vector_type &u, const vector_type &uA)
    {
        vector_type tmp1, tmp2, tmp3 ;
        tmp1 = u  / ietl::two_norm(uA) ;
        tmp2 = uA / ietl::two_norm(uA) ;
        tmp3 = apply_operator(uA) ;
        u_and_uA_.push_back(std::make_pair(std::make_pair(tmp1,tmp2),tmp3)) ;
    }
    // Update the vector space in JCD iteration
    template <class Matrix, class VS, class ITER>
    void jacobi_davidson_modified<Matrix, VS, ITER>::update_vecspace(vector_space& V, vector_space& VA, const int idx, vector_pairs& res)
    {
        vector_type t  = V[idx] ;
        vector_type tA = apply_operator(t);
        //
        for (typename vector_ortho_vec::iterator it = ortho_space_right_.begin(); it != ortho_space_right_.end(); it++) {
            t  -= ietl::dot((*it).first, tA) * (*it).second;
            tA -= ietl::dot((*it).first, tA) * (*it).first ;
        }
        //
        for (int i = 1; i <= idx ; i++) {
            t -= ietl::dot(VA[i-1], tA) * V[i-1];
            tA -= ietl::dot(VA[i-1], tA) * VA[i-1];
        }
        //
        V[idx]  = t /ietl::two_norm(tA) ;
        VA[idx] = tA/ietl::two_norm(tA) ;
    };
    // Get the matrix element of the Hamiltionian
    template <class Matrix, class VS, class ITER>
    typename jacobi_davidson_modified<Matrix, VS, ITER>::magnitude_type
             jacobi_davidson_modified<Matrix, VS, ITER>::get_matrix_element(const vector_type& V, const vector_type& VA)
    {
        return ietl::dot(V,VA) ;
    }
    // Compute the error vector
    template <class Matrix, class VS, class ITER>
    typename jacobi_davidson_modified<Matrix, VS, ITER>::vector_type jacobi_davidson_modified<Matrix,VS,ITER>::compute_error(const vector_type &u,
                                                                                                                             const vector_type &uA,
                                                                                                                             magnitude_type theta)
    {
        vector_type r = uA, tmp = u ;
        r -= tmp / theta;
        for (typename vector_ortho_vec::iterator it = ortho_space_right_.begin(); it != ortho_space_right_.end(); it++) 
            r    -= ietl::dot((*it).first, r) * (*it).first ;
        return r ;
    }
    // Check if the JD iteration is arrived at convergence
    template <class Matrix, class VS, class ITER>
    bool jacobi_davidson_modified<Matrix, VS, ITER>::check_convergence(const vector_type &u, const vector_type &uA, const vector_type &r,
                                                                       const magnitude_type theta, ITER& iter, vector_type &eigvec, magnitude_type &eigval)
    {
        // Compute the error vector
        bool converged ;
        eigvec = uA/ietl::two_norm(uA);
        eigval = this->omega_ - 1./theta ;
        if(iter.finished(ietl::two_norm(r),1.0)) {
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
                                                                            magnitude_type &theta,
                                                                            matrix_double& eigvecs,
                                                                            vector_double& eigvals)
    {
        // Initialization
        int imin , imax , nevec ;
        // Definition of the dimensions and dynamic memory allocation
        if (dim != n_restart_max_) {
            imin  = imax = 1;
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
    void jacobi_davidson_modified<MATRIX, VS, ITER>::solver(const vector_type& u,
                                                            const vector_type& uA,
                                                            const magnitude_type& theta,
                                                            const vector_type& r,
                                                            vector_type& t)
    {
        vector_type inh = -r, t2 ;
        scalar_type dru, duu ;
        gmres_modified<MATRIX, vector_type, VS> gmres(this->matrix_, u, vecspace_, uA, theta, ortho_space_left_, ortho_space_right_,
                                                      i_state_, omega_, max_iter_, false);
        // initial guess for better convergence
        if (i_gmres_guess_ == 0 || max_iter_ <= 1 ) {
            dru = ietl::dot(r, u);
            duu = ietl::dot(u, u);
            t = (-r + dru / duu * u);
        } else if (i_gmres_guess_ == 1) {
            t = 0.*r ;
        }
        if (max_iter_ > 0) {
            inh -= ietl::dot(inh, uA) * uA / ietl::dot(uA, uA) ;
            t2 = gmres(inh, t, ietl_atol_, ietl_rtol_) ;
            t = t2 / ietl::two_norm(t2);
        }
    }
    template<class MATRIX, class VS, class ITER>
    typename jacobi_davidson_modified<MATRIX, VS, ITER>::vector_double
    jacobi_davidson_modified<MATRIX, VS, ITER>::generate_property(const vector_space& V, const vector_space& VA, const size_t& ndim,
                                                                  const matrix_double &props, const vector_double &eigvals)
    {
        vector_double prop_values(ndim) ;
        for (int i = 0; i < ndim ; i++)
            prop_values[i] = omega_-1./eigvals[i];
        return prop_values ;
    }
    //
    template<class MATRIX, class VS, class ITER>
    void jacobi_davidson_modified<MATRIX, VS, ITER>::sort_prop(couple_vec& vector_values)
    {
        std::sort(vector_values.begin(),vector_values.end(),lt_couple()) ;
    } ;
    //
    template<class MATRIX, class VS, class ITER>
    void jacobi_davidson_modified<MATRIX, VS, ITER>::print_header_table() {
        print_endline() ;
        std::cout << " Iteration |    Error    |    Energy    " << std::endl ;
        print_endline() ;
    } ;
    //
    template<class MATRIX, class VS, class ITER>
    void jacobi_davidson_modified<MATRIX, VS, ITER>::print_endline() {
        std::cout << "-----------+-------------+--------------" << std::endl ;
    } ;
    //
    template<class MATRIX, class VS, class ITER>
    void jacobi_davidson_modified<MATRIX, VS, ITER>::print_newline_table(const size_t& i,          const double& error,
                                                                         const magnitude_type& en, const double& overlap )
    {
        char buf[39];
	    int a = i , n;
        n = sprintf(buf, "%5d      | %1.4E  | %6.5f ", a, error, en);
        std::cout << buf << std::endl;
    }
}
#endif
