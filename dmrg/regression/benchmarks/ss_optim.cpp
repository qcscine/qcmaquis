/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Bela Bauer <bauerb@comp-phys.org>
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

#ifdef USE_AMBIENT
#include <mpi.h>
#endif
#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/shared_ptr.hpp>

#include <alps/hdf5.hpp>

#include "matrix_selector.hpp" /// define matrix
#include "symm_selector.hpp"   /// define grp

#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"

#include "dmrg/mp_tensors/mps.h"

#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/optimize/ietl_jacobi_davidson.h"

#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"

#include "utils/timings.h"


template<class Matrix, class SymmGroup>
struct SiteProblem
{
    SiteProblem(Boundary<Matrix, SymmGroup> const & left_,
                Boundary<Matrix, SymmGroup> const & right_,
                MPOTensor<Matrix, SymmGroup> const & mpo_)
    : left(left_)
    , right(right_)
    , mpo(mpo_)
    {
    }
    
    Boundary<Matrix, SymmGroup> const & left;
    Boundary<Matrix, SymmGroup> const & right;
    MPOTensor<Matrix, SymmGroup> const & mpo;
    double ortho_shift;
};


int main(int argc, char ** argv)
{
    try {
        DmrgOptions opt(argc, argv);
        if (!opt.valid) return 0;
        DmrgParameters parms = opt.parms;
        
        maquis::cout.precision(10);
        
        /// Timers
        Timer tim_model      ("Parsing model"),  tim_load        ("Load MPS");
        Timer tim_l_boundary ("Left boundary"),  tim_r_boundary  ("Right boundary");
        Timer tim_optim_jcd  ("Optim. JCD"   ),  tim_optim_alpha ("Optim. alpha");
        
        /// Parsing model
        tim_model.begin();
        Lattice lattice(parms);
        Model<matrix, grp> model(lattice, parms);
        
        MPO<matrix, grp> mpo = make_mpo(lattice, model, parms);
        tim_model.end();
        
        
        /// Initialize & load MPS
        tim_load.begin();
        int L = lattice.size();
        MPS<matrix, grp> mps;
        load(parms["chkpfile"].str(), mps);
        int _site;
        {
            alps::hdf5::archive ar(parms["chkpfile"].str()+"/props.h5");
            ar["/status/site"] >> _site;
        }
        int site, lr;
        if (_site < L) {
            site = _site;
            lr = 1;
        } else {
            site = 2*L-_site-1;
            lr = -1;
        }
        maquis::cout << "Optimization at site " << site << " in " << lr << " direction." << std::endl;
        tim_load.end();
        
        /// Canonize MPS
        mps.canonize(site);
        
        /// Compute left boundary
        tim_l_boundary.begin();
        Boundary<matrix, grp> left = mps.left_boundary();
        for (size_t i=0; i<site; ++i)
            left = contraction::overlap_mpo_left_step(mps[i], mps[i], left, mpo[i]);
        tim_l_boundary.end();
        
        /// Compute right boundary
        tim_r_boundary.begin();
        Boundary<matrix, grp> right = mps.right_boundary();
        for (int i=L-1; i>site; --i)
            right = contraction::overlap_mpo_right_step(mps[i], mps[i], right, mpo[i]);
        tim_r_boundary.end();


        /// Optimization
        std::vector<MPSTensor<matrix, grp> > ortho_vecs;
        std::pair<double, MPSTensor<matrix, grp> > res;
        SiteProblem<matrix, grp> sp(left, right, mpo[site]);
        
        /// Optimization: JCD
        tim_optim_jcd.begin();
        res = solve_ietl_jcd(sp, mps[site], parms, ortho_vecs);
        mps[site] = res.second;
        maquis::cout.precision(10);
        maquis::cout << "Energy " << lr << " " << res.first << std::endl;
        tim_optim_jcd.end();
        
        double alpha = parms["alpha_main"];
        double cutoff = parms["truncation_final"];
        std::size_t Mmax = parms["max_bond_dimension"];
        
        /// Optimization: grow alpha
        std::pair<std::size_t, double> trunc;
        tim_optim_alpha.begin();
        if (lr == +1) {
            if (site < L-1) {
                maquis::cout << "Growing, alpha = " << alpha << std::endl;
                mps.grow_l2r_sweep(mpo[site], left, right, site, alpha, cutoff, Mmax);
            } else {
                block_matrix<matrix, grp> t = mps[site].normalize_left(DefaultSolver());
                if (site < L-1)
                    mps[site+1].multiply_from_left(t);
            }
        } else if (lr == -1) {
            if (site > 0) {
                maquis::cout << "Growing, alpha = " << alpha << std::endl;
                mps.grow_r2l_sweep(mpo[site], left, right, site, alpha, cutoff, Mmax);
            } else {
                block_matrix<matrix, grp> t = mps[site].normalize_right(DefaultSolver());
                if (site > 0)
                    mps[site-1].multiply_from_right(t);
            }
        }
        tim_optim_alpha.end();
        
        /// Compute new boundary
        // TODO: optional here...
        
        
    } catch (std::exception & e) {
        maquis::cerr << "Exception caught:" << std::endl << e.what() << std::endl;
        exit(1);
    }
}

