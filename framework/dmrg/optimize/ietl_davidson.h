/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef IETL_DAVIDSON_SOLVER_H
#define IETL_DAVIDSON_SOLVER_H

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/optimize/vectorset.h"

#include "davidson_standard.h"
#include "davidson_modified.h"
#include "davidson_standard_mo.h"
#include "davidson_modified_mo.h"

template<class Matrix, class SymmGroup, class BoundDatabase>
std::vector< std::pair< double , class MPSTensor<Matrix,SymmGroup> > >
solve_ietl_davidson(SiteProblem<Matrix, SymmGroup> & sp,
                    VectorSet<Matrix, SymmGroup> & initial,
                    BaseParameters & params,
                    std::vector<partial_overlap<Matrix, SymmGroup> > poverlap_vec ,
                    int nsites,
                    int site1,
                    int site2,
                    int root_homing_type,
                    std::vector< std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > > vec_sa_left,
                    std::vector< std::vector< std::vector<block_matrix<typename storage::constrained<Matrix>::type, SymmGroup> > > > vec_sa_right,
                    const std::vector< int > & order, 
                    BoundDatabase bound_database,
                    int sa_alg,
                    std::vector<class MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector< class MPSTensor<Matrix, SymmGroup> >())
{
    // Initialization
    typedef MPSTensor<Matrix, SymmGroup> Vector ;
    SingleSiteVS<Matrix, SymmGroup> vs(initial, ortho_vecs, vec_sa_left, vec_sa_right, bound_database);
    int n_restart = params["ietl_diag_restart_nmax"] ;
    double atol   = params["ietl_diag_atol"];
    double rtol   = params["ietl_diag_rtol"];
    double omega  = params["ietl_si_omega"] ;
    ietl::basic_iteration<double> iter(params["ietl_diag_maxiter"], rtol, atol);
    std::vector < std::pair<double, Vector> > r0 ;
    // Check if the number of MPSTensors is higher than the one of the orthogonal vectors
    // and performs the GS orthogonalization
    if (initial.MPSTns_SA[0].num_elements() <= ortho_vecs.size())
        ortho_vecs.resize(initial.MPSTns_SA[0].num_elements()-1);
    for (int n = 1; n < ortho_vecs.size(); ++n)
        for (int n0 = 0; n0 < n; ++n0)
            ortho_vecs[n] -= ietl::dot(ortho_vecs[n0], ortho_vecs[n])/ietl::dot(ortho_vecs[n0],ortho_vecs[n0])*ortho_vecs[n0] ;
    // Check orthogonality
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial.MPSTns_SA[0], ortho_vecs[n]) << std::endl ;
    }
    // -- Calculation of eigenvalues
    // TODO Alb - here the choice is done based on the numerical value of omega, might be done better
    if (fabs(omega) > 1.0E-15) {
        if ( root_homing_type > 0 ) {
            ietl::davidson_modified_mo<SiteProblem <Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, Matrix, SymmGroup >
                    davidson(sp, vs, omega, poverlap_vec, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"],
                             nsites, site1, site2, root_homing_type);
            r0 = davidson.calculate_eigenvalue(iter);
        } else {
            ietl::davidson_modified< SiteProblem< Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
                    davidson(sp, vs, omega, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"],
                             nsites, site1, site2);
            r0 = davidson.calculate_eigenvalue(iter);
        }
    } else {
        if ( root_homing_type > 0 ) {
            ietl::davidson_standard_mo<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, Matrix, SymmGroup >
                    davidson(sp, vs, poverlap_vec, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"],
                             nsites, site1, site2, root_homing_type);
            r0 = davidson.calculate_eigenvalue(iter);
        } else {
            ietl::davidson_standard< SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
                    davidson(sp, vs, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"],
                             nsites, site1, site2);
            r0 = davidson.calculate_eigenvalue(iter);
        }
    }
    // Check again orthogonality in output
    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0[0].second, ortho_vecs[n]) << std::endl;
    maquis::cout << "Davidson used " << iter.iterations() << " iterations." << std::endl;
    return r0;
}



#endif
