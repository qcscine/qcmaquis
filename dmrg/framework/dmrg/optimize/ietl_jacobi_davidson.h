/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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

#ifndef IETL_JD_SOLVER_H
#define IETL_JD_SOLVER_H

#include "ietl_lanczos_solver.h"
#include "ietl/jd.h"

#include "dmrg/optimize/utils/iter_jacobi.h"
#include "dmrg/optimize/jacobi.h"
#include "dmrg/optimize/jacobi_standard.h"
#include "dmrg/optimize/jacobi_modified.h"
//Leon: Comment partial_overlap because it doesn't support SU2U1 yet
//#include "dmrg/optimize/jacobi_standard_mo.h"
//#include "dmrg/optimize/jacobi_modified_mo.h"
//#include "dmrg/optimize/partial_overlap.h"
#include "dmrg/optimize/vectorset.h"

//
// Jacobi-Davidson diagonalization
// -------------------------------

template<class Matrix, class SymmGroup, class BoundDatabase>
std::vector< std::pair<typename maquis::traits::real_type<typename Matrix::value_type>::type, class MPSTensor<Matrix,SymmGroup> > >
solve_ietl_jcd(SiteProblem<Matrix, SymmGroup> & sp,
               VectorSet<Matrix, SymmGroup> & initial,
               BaseParameters & params ,
//               std::vector< partial_overlap<Matrix, SymmGroup> >  poverlap_vec,
               int nsites,
               int site1,
               int site2,
               int root_homing_type,
               bool do_shiftandinvert,
               std::vector< std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > > vec_sa_left,
               std::vector< std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > > vec_sa_right,
               std::vector< int > const & order,
               BoundDatabase bound_database,
               int sa_alg,
               std::vector< double > omega_vec,
               std::vector< class MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector< class MPSTensor<Matrix, SymmGroup> >())
{
    // -- Initialization --
    typedef MPSTensor<Matrix, SymmGroup> Vector;

    double rtol  = params["ietl_diag_rtol"] ;
    double atol  = params["ietl_diag_atol"] ;
    int n_tofollow = params["maximum_overlap_nstates"] ;
    int n_restart  = params["ietl_diag_restart"] ;
    int side_tofollow = params["maximum_overlap_side"] ;
    std::vector < std::pair<double, Vector> > r0 ;
    ietl::jcd_iterator<double> iter(params["ietl_diag_maxiter"], rtol, atol);
    // -- Orthogonalization --
    if (initial.MPSTns_SA[0].num_elements() <= ortho_vecs.size())
        ortho_vecs.resize(initial.MPSTns_SA[0].num_elements()-1);
    for (int n = 1; n < ortho_vecs.size(); ++n)
        for (int n0 = 0; n0 < n; ++n0)
            ortho_vecs[n] -= ietl::dot(ortho_vecs[n0], ortho_vecs[n])/ietl::dot(ortho_vecs[n0],ortho_vecs[n0])*ortho_vecs[n0];
    // --  JD ALGORITHM --
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial.MPSTns_SA[0], ortho_vecs[n]) << std::endl ;
    }
    // -- GMRES STARTING GUESS --
    size_t i_gmres_guess ;
    if (params["ietl_gmres_guess"] == "error")
        i_gmres_guess = 0 ;
    else if (params["ietl_gmres_guess"] == "zero")
        i_gmres_guess = 1 ;
    else
        throw std::runtime_error("Guess for the GMRES procedure not recognized") ;
    // +----------------------+
    //  EIGENVALUE CALCULATION
    // +----------------------+
    SingleSiteVS<Matrix, SymmGroup> vs(initial, ortho_vecs, vec_sa_left, vec_sa_right, bound_database);
    if ( !do_shiftandinvert ) {
        if ( root_homing_type == 0) {
            ietl::jacobi_davidson_standard<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, ietl::jcd_iterator<double> >
                jd(sp, vs, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"], params["ietl_gmres_maxiter"],
                   nsites, site1, site2, params["ietl_gmres_abstol"], params["ietl_gmres_reltol"], i_gmres_guess, order, sa_alg) ;
            r0 = jd.calculate_eigenvalue(initial, iter) ;
        } else {
/*
			    ietl::jacobi_davidson_standard_mo<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, ietl::jcd_iterator<double> , Matrix, SymmGroup >
                jd(sp, vs, poverlap_vec, n_tofollow, side_tofollow, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"], params["ietl_gmres_maxiter"],
                   nsites, site1, site2, params["ietl_gmres_abstol"], params["ietl_gmres_reltol"], i_gmres_guess, order, sa_alg, root_homing_type) ;
            r0 = jd.calculate_eigenvalue(initial, iter);*/
            throw std::runtime_error("root homing type!=0/partial overlap not implemented for SU2U1");
        }
    } else {
        if ( root_homing_type == 0 ) {
            ietl::jacobi_davidson_modified<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, ietl::jcd_iterator<double> >
                jd(sp, vs, omega_vec, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"], params["ietl_gmres_maxiter"],
                   nsites, site1, site2, params["ietl_gmres_abstol"], params["ietl_gmres_reltol"], i_gmres_guess, order, sa_alg, params["ietl_gmres_init_atol"],
                   params["ietl_gmres_init_rtol"], params["ietl_gmres_init_maxiter"] );
            r0 = jd.calculate_eigenvalue(initial, iter) ;
        } else {
          /*  ietl::jacobi_davidson_modified_mo<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, ietl::jcd_iterator<double> , Matrix, SymmGroup>
                    jd(sp, vs, omega_vec, poverlap_vec, n_tofollow, side_tofollow, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"], params["ietl_gmres_maxiter"],
                       nsites, site1, site2, params["ietl_gmres_abstol"], params["ietl_gmres_reltol"], i_gmres_guess, order, sa_alg, params["ietl_gmres_init_atol"],
                       params["ietl_gmres_init_rtol"], params["ietl_gmres_init_maxiter"], root_homing_type);
            r0 = jd.calculate_eigenvalue(initial, iter);*/
            throw std::runtime_error("root homing type!=0/partial overlap not implemented for SU2U1");
        }
    }
    // Prints final header
    std::cout << "\n Summary of the results " << std::endl ;
    std::cout <<   " ---------------------- " << std::endl ;
    maquis::cout << " Number of iterations - " << iter.iterations() << std::endl;
    return r0;
}

#endif
