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
#include "dmrg/optimize/partial_overlap.h"
#include "ietl_lanczos_solver.h"

#include "davidson_standard.h"
#include "davidson_modified.h"
#include "davidson_standard_mo.h"
#include "davidson_modified_mo.h"

template<class Matrix, class SymmGroup>
std::pair< double , class MPSTensor<Matrix,SymmGroup> >
solve_ietl_davidson(SiteProblem<Matrix, SymmGroup> & sp,
                    MPSTensor<Matrix, SymmGroup> const & initial,
                    BaseParameters & params,
                    int site,
                    partial_overlap<Matrix, SymmGroup> poverlap,
                    std::vector<class MPSTensor<Matrix, SymmGroup> > ortho_vecs = std::vector< class MPSTensor<Matrix, SymmGroup> >()) {
    // Initialization
    typedef MPSTensor<Matrix, SymmGroup> Vector ;
    SingleSiteVS<Matrix, SymmGroup> vs(initial, ortho_vecs);
    int n_restart = params["ietl_diag_restart_nmax"] ;
    double atol   = params["ietl_diag_atol"];
    double rtol   = params["ietl_diag_rtol"];
    double omega  = params["ietl_si_omega"] ;
    ietl::basic_iteration<double> iter(params["ietl_diag_maxiter"], rtol, atol);
    std::pair<double, Vector> r0 ;
    // Check if the number of MPSTensors is higher than the one of the orthogonal vectors
    // and performs the GS orthogonalization
    if (initial.num_elements() <= ortho_vecs.size())
    ortho_vecs.resize(initial.num_elements()-1);
    for (int n = 1; n < ortho_vecs.size(); ++n)
        for (int n0 = 0; n0 < n; ++n0)
            ortho_vecs[n] -= ietl::dot(ortho_vecs[n0], ortho_vecs[n])/ietl::dot(ortho_vecs[n0],ortho_vecs[n0])*ortho_vecs[n0];
    // Check orthogonality
    for (int n = 0; n < ortho_vecs.size(); ++n) {
        maquis::cout << "Input <MPS|O[" << n << "]> : " << ietl::dot(initial, ortho_vecs[n]) << std::endl;
    }
    // -- Calculation of eigenvalues
    // TODO Alb - here the choice is done based on the numerical value of omega, might be done better
    if (fabs(omega) > 1.0E-15) {
        if ( poverlap.is_defined()) {
            ietl::davidson_modified_mo<SiteProblem <Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, Matrix, SymmGroup >
                    davidson(sp, vs, site, omega, poverlap, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"]);
            r0 = davidson.calculate_eigenvalue(initial, iter);
        } else {
            ietl::davidson_modified<SiteProblem < Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
                    davidson(sp, vs, site, omega, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"]);
            r0 = davidson.calculate_eigenvalue(initial, iter);
        }
    } else {
        if ( poverlap.is_defined()) {
            ietl::davidson_standard_mo<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup>, Matrix, SymmGroup >
                    davidson(sp, vs, site, poverlap, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"]);
            r0 = davidson.calculate_eigenvalue(initial, iter);
        } else {
            ietl::davidson_standard<SiteProblem<Matrix, SymmGroup>, SingleSiteVS<Matrix, SymmGroup> >
                    davidson(sp, vs, site, params["ietl_diag_restart_nmin"], params["ietl_diag_restart_nmax"]);
            r0 = davidson.calculate_eigenvalue(initial, iter);
        }
    }
    // Check again orthogonality in output
    for (int n = 0; n < ortho_vecs.size(); ++n)
        maquis::cout << "Output <MPS|O[" << n << "]> : " << ietl::dot(r0.second, ortho_vecs[n]) << std::endl;
    maquis::cout << "Davidson used " << iter.iterations() << " iterations." << std::endl;
    return r0;
}



#endif
