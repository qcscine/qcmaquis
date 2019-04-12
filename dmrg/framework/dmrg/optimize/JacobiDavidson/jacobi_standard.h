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
 * THE SOFTWARE IS PROVIDED "AS IS", WITHoUT WARRANTY OF ANY KIND, EXPRESS OR
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

#include <ietl/cg.h>
#include <ietl/fmatrix.h>
#include <ietl/ietl2lapack.h>
#include <ietl/traits.h>
#include <vector>
#include "jacobi.h"
#include <boost/range/adaptor/reversed.hpp>

// +---------------------------+
//  JACOBI-DAVIDSON EIGENSOLVER
// +---------------------------+
namespace ietl
{
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    class jacobi_davidson_standard : public jacobi_davidson<MATRIX, VS, SymmGroup, ITER>
    {
    public:
        typedef jacobi_davidson<MATRIX, VS, SymmGroup, ITER> base;
        //
        typedef typename base::CorrectionEquation CorrectionEquation;
        typedef typename base::couple_vec         couple_vec;
        typedef typename base::Finalizer          Finalizer;
        typedef typename base::lt_couple          lt_couple;
        typedef typename base::magnitude_type     magnitude_type;
        typedef typename base::MicroOptimizer     MicroOptimizer;
        typedef typename base::Orthogonalizer     Orthogonalizer;
        typedef typename base::real_type          real_type ;
        typedef typename base::scalar_type        scalar_type;
        typedef typename base::vector_double      vector_double;
        typedef typename base::vector_space       vector_space;
        typedef typename base::vector_type        vector_type;
        //
        using base::atol_ ;
        using base::bounds ;
        using base::candidates_collector_ ;
        using base::corrector_ ;
        using base::diagonal_elements_ ;
        using base::eigen_collector_ ;
        using base::energy_ref_ ;
        using base::estimate_extremes ;
        using base::finalizer_ ;
        using base::highest_eigen_ ;
        using base::i_state_ ;
		using base::is_folded_ ;
        using base::lowest_eigen_ ;
        using base::M ;
        using base::matrix_ ;
        using base::micro_iterator_ ;
        using base::n_lanczos_ ;
        using base::n_restart_max_ ;
        using base::n_restart_min_ ;
        using base::n_root_found_ ;
        using base::n_sa_ ;
        using base::order_ ;
        using base::orthogonalizer_ ;
        using base::sa_alg_ ;
        using base::site1_ ;
        using base::site2_ ;
        using base::u_and_uA_ ;
        using base::vecspace_ ;
        using base::V_ ;
        using base::VA_ ;
        using base::diagonalize;
        //
        jacobi_davidson_standard(MATRIX& matrix, VS& vec, CorrectionEquation& corrector,std::unique_ptr<MicroOptimizer>& micro_iterator,
                                 Finalizer& finalizer, std::unique_ptr<Orthogonalizer> & ortho, const size_t& nmin, const size_t& nmax,
                                 const size_t& n_block, const double& block_thresh, const int& site1, const int& site2,
                                 const std::vector<std::size_t>& order, const int& sa_alg, const int& n_lanczos,
                                 const bool& do_chebychev, const magnitude_type& chebyshev_shift, const bool& do_H_squared,
                                 const bool& reshuffle_variance, const bool& track_variance, const bool& is_folded,
								 const double& energy_thresh)
                : base::jacobi_davidson(matrix, vec, corrector, micro_iterator, finalizer, ortho, nmin, nmax, n_block,
                                        block_thresh, site1, site2, order, sa_alg, n_lanczos, do_chebychev, chebyshev_shift,
                                        do_H_squared, reshuffle_variance, track_variance, is_folded, energy_thresh)
        {
            bounds.do_lower     = true ;
            bounds.lower_bounds = std::make_pair(1,1) ;
        } ;
        ~jacobi_davidson_standard() {} ;
    protected:
        vector_type apply_operator (const vector_type& x);
        void update_finalizer() ;
        void update_parameters() ;
        void update_vecspace(vector_space& to_add);
        void set_interval(const std::size_t& dim) {} ;
        virtual void get_eigenvalue(std::vector<double>& eigval, std::vector<class std::vector<double> >& eigvecs, size_t dim);
        virtual void do_diagonalize()
        {
            // only diagonalise for the lowest eigenvalues here
            diagonalize(true);
        }

    private:
        bool check_convergence(const size_t& idx, ITER& iter);
        vector_double generate_property() ;
        void print_endline() ;
        void print_header_table() ;
        void print_newline_table(const size_t& i, const real_type& error, const scalar_type& en, const size_t& idx,
                                 const bool& converged) ;
        void solver(const vector_type& r, vector_space& t) ;
        void sort_prop(couple_vec& vector_values) ;
    } ;
    // New version for generation of the guess
    template <class Matrix, class VS, class SymmGroup, class ITER>
    void  jacobi_davidson_standard<Matrix, VS, SymmGroup, ITER>::update_parameters()
    {
        energy_ref_ = std::real(ietl::dot(V_[0], VA_[0]) / ietl::dot(V_[0], V_[0])) ;
        corrector_.update_u(V_[0]) ;
    };
    // Computes the action of an operator
    template <class Matrix, class VS, class SymmGroup, class ITER>
    typename jacobi_davidson_standard<Matrix, VS, SymmGroup, ITER>::vector_type
             jacobi_davidson_standard<Matrix, VS, SymmGroup, ITER>::apply_operator(vector_type const & x)
    {
        vector_type y = x , y2 ;
        ietl::mult(this->matrix_, y, y2, i_state_, is_folded_) ;
        return y2 ;
    };
    // Update the finalizer object
    template <class Matrix, class VS, class SymmGroup, class ITER>
    void jacobi_davidson_standard<Matrix, VS, SymmGroup, ITER>::update_finalizer()
    {
        finalizer_.set_candidate(eigen_collector_) ;
    }
    // Update the vector space in JCD iteration
    template <class Matrix, class VS, class SymmGroup, class ITER>
    void jacobi_davidson_standard<Matrix, VS, SymmGroup, ITER>::update_vecspace(vector_space& to_add)
    {
        for (auto&& i : boost::adaptors::reverse(to_add))
        {
            vector_type&& t = std::move(i);
            vecspace_.project(t) ;
            vector_type&& tA = apply_operator(t) ;
            vecspace_.project(tA) ;
            orthogonalizer_->orthogonalize(t, tA) ;
            if (std::abs(ietl::two_norm(t)) > 1.0E-10) {
                orthogonalizer_->normalize(t, tA) ;
                V_.emplace_back(t) ;
                VA_.emplace_back(tA) ;
            }
        }
        to_add.clear();
    };
    // Check if the JD iteration is arrived at convergence
    template <class Matrix, class VS, class SymmGroup, class ITER>
    bool jacobi_davidson_standard<Matrix, VS, SymmGroup, ITER>::check_convergence(size_t const& idx,
                                                                                  ITER& iter)
    {
        // Compute the error vector
        bool converged ;
        vector_type jnk ;
        ietl::mult(this->matrix_, eigen_collector_[idx].u_, jnk, i_state_, false) ;
        magnitude_type energy = std::real(ietl::dot(eigen_collector_[idx].u_, jnk) / ietl::dot(eigen_collector_[idx].u_,
                                                                                     eigen_collector_[idx].u_) ) ;
        vector_type r = jnk - energy*eigen_collector_[idx].u_ ;
        vecspace_.project(r) ;

        return iter.finished(ietl::two_norm(r),std::fabs(energy));

    };
    // -- SOLVER --
    template<class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER>::solver(const vector_type& r, vector_space& t)
    {
        vector_type t2 = 0.*r ;
        micro_iterator_->set_error(r) ;
        t2 = micro_iterator_->PerformOptimization(t2) ;
        t.emplace_back(t2 / ietl::two_norm(t2));
    } ;
    // -- GENERATE_PROPERTY --
    template<class MATRIX, class VS, class SymmGroup, class ITER>
    typename jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER>::vector_double
             jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER>::generate_property()
    {
        // Variables declaration
        std::size_t dim = candidates_collector_.size() ;
        vector_double vector_values(dim) ;
        // Property calculation
        for (size_t i = 0; i < dim ; i++)
            vector_values[i]  = candidates_collector_[i].theta_ ;
        return vector_values ;
    } ;
    //
    template<class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER>::sort_prop(couple_vec& vector_values)
    {
        std::sort(vector_values.begin(),vector_values.end(),lt_couple()) ;
    }
    //
    template<class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER>::print_header_table()
    {
        print_endline() ;
        std::cout << " Iteration | Sub. Dim. |    Error    |    Energy    " << std::endl ;
        print_endline() ;
    } ;
    //
    template<class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER>::print_endline()
    {
        std::cout << "-----------+-----------+-------------+-------------" << std::endl ;
    } ;
    //
    template<class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER>::print_newline_table
            (const size_t& i, const real_type& error, const scalar_type& en, const size_t& idx, const bool& converged)

    {
        char buf[60];
        int a = i , n ;
        n = sprintf(buf, "%5d      |%5d      | %3.4E  | %6.5F ", a, this->get_dim(V_) , error, en);
        std::cout << buf ;
        if (converged)
            std::cout << " CONVERGED" ;
        std::cout << std::endl;
    }

    // +-----------------------------------------------+
    //  Interface to the LAPACK diagonalization routine
    //  here we use the symmetric diagonalisation
    //  and calculate only one eigenvalue
    // +-----------------------------------------------+
    template <class MATRIX, class VS, class SymmGroup, class ITER>
    void jacobi_davidson_standard<MATRIX, VS, SymmGroup, ITER>::get_eigenvalue(std::vector<double>& eigval,
                                                                      std::vector<class std::vector<double> >& eigvecs,
                                                                      std::size_t dim)
    {
        FortranMatrix<magnitude_type> M_(dim,dim);
        for (int i=0;i<dim;i++)
        for (int j=0;j<=i;j++)
            M_(j,i) = M(j,i);
        double abstol = atol_;
        char jobz='V';     char range='I';   char uplo='U';
        fortran_int_t n=dim;
        fortran_int_t lda=dim;
        fortran_int_t il=1, iu=1;

        fortran_int_t m;
        fortran_int_t ldz=n;
        fortran_int_t lwork=8*n;
        fortran_int_t info;

        double vl, vu;
        double *w = new double[n];
        double *z = new double[n*n];
        double *work = new double[lwork];
        fortran_int_t *iwork = new fortran_int_t[5*n];
        fortran_int_t *ifail = new fortran_int_t[n];

        LAPACK_DSYEVX(&jobz, &range, &uplo, &n, M_.data(), &lda, &vl, &vu, &il, &iu, &abstol, &m, w, z, &ldz, work, &lwork, iwork, ifail, &info);

        // use only one eigenvalue
        eigval[0] = w[0];
        for (int i = 0; i < n; i++)
            eigvecs[0][i] = z[i];


        delete [] w;
        delete [] z;
        delete [] work;
        delete [] iwork;
        delete [] ifail;
    }
}

#endif
