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
#include "dmrg/optimize/JacobiDavidson/jacobi.h"
#include "dmrg/optimize/POverlap/partial_overlap.h"
#include "dmrg/optimize/vectorset.h"
#include "dmrg/optimize/CorrectionEquation/correctionequation.h"
#include "dmrg/optimize/Finalizer/finalizer.h"
#include "dmrg/optimize/Orthogonalizer/orthogonalizer.h"
#include "dmrg/optimize/utils/orthogonalizer_collector.h"

//
// Jacobi-Davidson diagonalization
// -------------------------------

template<class Matrix, class SymmGroup, class BoundDatabase, class MicroOptimizer>
std::vector< std::pair<typename maquis::traits::real_type<typename Matrix::value_type>::type, class MPSTensor<Matrix,SymmGroup> > >
solve_ietl_jcd(SiteProblem<Matrix, SymmGroup> & sp,
               VectorSet<Matrix, SymmGroup> & initial,
               CorrectionEquation<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >& correction_equation ,
               std::shared_ptr<MicroOptimizer>& micro_optimizer,
               Finalizer< SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >& finalizer ,
               std::shared_ptr<Orthogonalizer<SingleSiteVS<Matrix, SymmGroup> > >& orthogonalizer,
               BaseParameters & params ,
               std::vector< partial_overlap<Matrix, SymmGroup> > poverlap_vec,
               int site1,
               int site2,
               int root_homing_type,
               bool do_root_homing, bool do_shiftandinvert, bool do_chebyshev,
               double chebyshev_shift,
               bool do_H_squared, bool reshuffle_variance, bool track_variance, bool is_folded_,
               std::vector< std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > > vec_sa_left,
               std::vector< std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > > vec_sa_right,
               std::vector< std::size_t > const & order,
               BoundDatabase bound_database,
               int sa_alg,
               std::vector< double > omega_vec,
               std::vector< class MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector< class MPSTensor<Matrix, SymmGroup> >())
{
    // -- Initialization --
    typedef MPSTensor<Matrix, SymmGroup> Vector;
    int n_tofollow = params["maximum_overlap_nstates"] ;
    int side_tofollow = params["maximum_overlap_side"] ;
    std::vector < std::pair<double, Vector> > r0 ;
    ietl::jcd_iterator<double> iter(params["ietl_diag_maxiter"], params["ietl_diag_rtol"], params["ietl_diag_atol"]);
    int n_lanczos = params["lanczos_max_iterations"] ;
    // -- Orthogonalization --
    if (initial.MPSTns_SA[0].num_elements() <= ortho_vecs.size())
        ortho_vecs.resize(initial.MPSTns_SA[0].num_elements()-1);
    for (int n = 1; n < ortho_vecs.size(); ++n)
        for (int n0 = 0; n0 < n; ++n0)
            ortho_vecs[n] -= ietl::dot(ortho_vecs[n0], ortho_vecs[n])/ietl::dot(ortho_vecs[n0],ortho_vecs[n0])*ortho_vecs[n0];
    // --  JD ALGORITHM --
    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial.MPSTns_SA[0], ortho_vecs[n]) << std::endl ;
    // -- CORRECTOR SETUP --
    orthogonalizer_collector< MPSTensor<Matrix, SymmGroup> > ortho_collector(ortho_vecs) ;
    SingleSiteVS<Matrix, SymmGroup> vs(initial, vec_sa_left, vec_sa_right, ortho_collector) ;
    finalizer.set_Hamiltonian(sp) ;
    // +----------------------+
    //  EIGENVALUE CALCULATION
    // +----------------------+
    if ( !do_shiftandinvert ) {
        if ( !do_root_homing ) {
            ietl::jacobi_davidson_standard<SiteProblem<Matrix, SymmGroup>,
                                           SingleSiteVS<Matrix, SymmGroup>,
                                           SymmGroup,
                                           ietl::jcd_iterator<double> >
            jd(sp, vs, correction_equation, micro_optimizer, finalizer, orthogonalizer, params["ietl_diag_restart_nmin"],
               params["ietl_diag_restart_nmax"], params["ietl_diag_block"], params["ietl_block_thresh"], site1,
               site2, order, sa_alg, n_lanczos, do_chebyshev, chebyshev_shift, do_H_squared, reshuffle_variance,
               track_variance, is_folded_, params["homing_energy_threshold"]) ;
            r0 = jd.calculate_eigenvalue(iter) ;
        } else {
            ietl::jacobi_davidson_standard_mo<SiteProblem<Matrix, SymmGroup>,
                                              SingleSiteVS<Matrix, SymmGroup>,
                                              ietl::jcd_iterator<double>,
                                              Matrix,
                                              SymmGroup >
            jd(sp, vs, correction_equation, micro_optimizer, finalizer, orthogonalizer, poverlap_vec, n_tofollow,
               side_tofollow, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"], params["ietl_diag_block"],
               params["ietl_block_thresh"], site1, site2, order, sa_alg, n_lanczos, do_chebyshev, chebyshev_shift,
               do_H_squared, reshuffle_variance, track_variance, is_folded_, params["homing_energy_threshold"], root_homing_type,
               params["uA_homing_ratio"]) ;
            r0 = jd.calculate_eigenvalue(iter);
        }
    } else {
        if ( !do_root_homing ) {
            ietl::jacobi_davidson_modified<SiteProblem<Matrix, SymmGroup>,
                                           SingleSiteVS<Matrix, SymmGroup>,
                                           SymmGroup,
                                           ietl::jcd_iterator<double> >
            jd(sp, vs, correction_equation, micro_optimizer, finalizer, orthogonalizer, omega_vec, params["ietl_diag_restart_nmin"],
               params["ietl_diag_restart_nmax"], params["ietl_diag_block"], params["ietl_block_thresh"], site1, site2,
               order, sa_alg, n_lanczos, do_chebyshev, chebyshev_shift, do_H_squared, reshuffle_variance, track_variance,
               is_folded_, params["homing_energy_threshold"]);
            r0 = jd.calculate_eigenvalue(iter) ;
        } else {
            ietl::jacobi_davidson_modified_mo<SiteProblem<Matrix, SymmGroup>,
                                              SingleSiteVS<Matrix, SymmGroup>,
                                              ietl::jcd_iterator<double> ,
                                              Matrix,
                                              SymmGroup>
            jd(sp, vs, correction_equation, micro_optimizer, finalizer, orthogonalizer, omega_vec,
               poverlap_vec, n_tofollow, side_tofollow, params["ietl_diag_restart_nmin"],
               params["ietl_diag_restart_nmax"], params["ietl_diag_block"], params["ietl_block_thresh"],
               site1, site2, order, sa_alg, n_lanczos, do_chebyshev, chebyshev_shift, do_H_squared, reshuffle_variance,
               track_variance, is_folded_, params["homing_energy_threshold"], root_homing_type, params["uA_homing_ratio"]);
            r0 = jd.calculate_eigenvalue(iter);
        }
    }
    // Prints final header
    std::cout << "\n Summary of the results " << std::endl ;
    std::cout <<   " ---------------------- " << std::endl ;
    maquis::cout << " Number of iterations - " << iter.iterations() << std::endl;
    return r0;
}

#endif
