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

#include "dmrg/optimize/jacobi_standard.h"
#include "dmrg/optimize/jacobi_modified.h"
#include "dmrg/optimize/jacobi_standard_mo.h"
#include "dmrg/optimize/jacobi_modified_mo.h"
#include "dmrg/optimize/partial_overlap.h"
#include "dmrg/optimize/vectorset.h"

//
// Jacobi-Davidson diagonalization
// -------------------------------

template<class Matrix, class SymmGroup>
std::vector< std::pair< double , class MPSTensor<Matrix,SymmGroup> > >
solve_ietl_jcd(SiteProblem<Matrix, SymmGroup> & sp,
               VectorSet<Matrix, SymmGroup> const & initial,
               BaseParameters & params ,
               std::vector< partial_overlap<Matrix, SymmGroup> >  poverlap_vec,
               int nsites,
               int site1,
               int root_homing_type,
               std::vector< class MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector< class MPSTensor<Matrix, SymmGroup> >(),
               int site2=0)
{
    // -- Initialization --
    typedef MPSTensor<Matrix, SymmGroup> Vector;
    double rtol  = params["ietl_diag_rtol"] ;
    double atol  = params["ietl_diag_atol"] ;
    double omega = params["ietl_si_omega"] ;
    int n_tofollow = params["maximum_overlap_nstates"] ;
    int n_restart  = params["ietl_diag_restart"] ;
    std::vector < std::pair<double, Vector> > r0 ;
    ietl::basic_iteration<double> iter(params["ietl_diag_maxiter"], rtol, atol);
    // -- Orthogonalization --
    if (initial.MPSTns_averaged.num_elements() <= ortho_vecs.size())
        ortho_vecs.resize(initial.MPSTns_averaged.num_elements()-1);
    for (int n = 1; n < ortho_vecs.size(); ++n)
        for (int n0 = 0; n0 < n; ++n0)
            ortho_vecs[n] -= ietl::dot(ortho_vecs[n0], ortho_vecs[n])/ietl::dot(ortho_vecs[n0],ortho_vecs[n0])*ortho_vecs[n0];
    // --  JD ALGORITHM --
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial.MPSTns_averaged, ortho_vecs[n]) << std::endl ;
    }
    // -- GMRES STARTING GUESS --
    size_t i_gmres_guess ;
    if (params["ietl_gmres_guess"] == "error")
        i_gmres_guess = 0 ;
    else if (params["ietl_gmres_guess"] == "zero")
        i_gmres_guess = 1 ;
    else
        throw std::runtime_error("Guess for the GMRES procedure not recognized") ;
    //
    SingleSiteVS<Matrix, SymmGroup> vs(initial, ortho_vecs);
    if (fabs(omega) < 1.0E-15) {
        if ( root_homing_type == 0) {
            ietl::jacobi_davidson_standard<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, ietl::basic_iteration<double> >
                jd(sp, vs, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"], params["ietl_gmres_maxiter"],
                   nsites, site1, site2, params["ietl_gmres_abstol"], i_gmres_guess) ;
            r0 = jd.calculate_eigenvalue(initial, iter) ;
        } else {
            ietl::jacobi_davidson_standard_mo<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, ietl::basic_iteration<double> , Matrix, SymmGroup >
                jd(sp, vs, poverlap_vec, n_tofollow, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"], params["ietl_gmres_maxiter"],
                   nsites, site1, site2, params["ietl_gmres_abstol"], i_gmres_guess, root_homing_type) ;
            r0 = jd.calculate_eigenvalue(initial, iter);
        }
    } else {
        if ( root_homing_type == 0 ) {
            ietl::jacobi_davidson_modified<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, ietl::basic_iteration<double> >
                    jd(sp, vs, omega, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"], params["ietl_gmres_maxiter"],
                       nsites, site1, site2, params["ietl_gmres_abstol"], i_gmres_guess);
            r0 = jd.calculate_eigenvalue(initial, iter) ;
        } else {
            //ietl::jacobi_davidson_modified_mo<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, ietl::basic_iteration<double> , Matrix, SymmGroup>
            //        jd(sp, vs, omega, poverlap, n_tofollow , params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"], params["ietl_gmres_maxiter"],
            //           nsites, site1, site2, params["ietl_gmres_abstol"], i_gmres_guess, root_homing_type);
            //r0 = jd.calculate_eigenvalue(initial, iter);
        }
    }
    std::cout << "\n Summary of the results " << std::endl ;
    std::cout <<   " ---------------------- " << std::endl ;
    maquis::cout << " Number of iterations - " << iter.iterations() << std::endl;
    return r0;
}


#endif
