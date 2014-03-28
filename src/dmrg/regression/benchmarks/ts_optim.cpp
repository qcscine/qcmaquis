/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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
#include "dmrg/mp_tensors/twositetensor.h"

#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/optimize/ietl_jacobi_davidson.h"

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/utils/placement.h"

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
        if (argc != 2 && argc != 3)
        {
            maquis::cout << "Usage: <parms> [<model_parms>]" << std::endl;
            exit(1);
        }
        
        maquis::cout.precision(10);
        
        /// Load parameters
        std::ifstream param_file(argv[1]);
        if (!param_file)
            throw std::runtime_error("Could not open parameter file.");
        DmrgParameters parms(param_file);
        
        /// Load model parameters from second input (if needed)
        std::string model_file;
        if (parms.is_set("model_file")) model_file = parms["model_file"].str();
        if (argc == 3)                  model_file = std::string(argv[2]);
        if (!model_file.empty()) {
            std::ifstream model_ifs(model_file.c_str());
            if (!model_ifs)
                throw std::runtime_error("Could not open model_parms file.");
            parms << ModelParameters(model_ifs);
        }
        
        /// Timers
        #ifdef USE_AMBIENT
        ambient::timer 
        #else
        Timer
        #endif
        tim_model      ("Parsing model"),  tim_load       ("Load MPS"),
        tim_l_boundary ("Left boundary"),  tim_r_boundary ("Right boundary"),
        tim_optim_jcd  ("Optim. JCD"   ),  tim_truncation ("Truncation"),
        tim_ts_obj     ("Two site obj.");
        
        /// Parsing model
        tim_model.begin();
        Lattice lattice(parms);
        Model<matrix, grp> model(lattice, parms);
        
        MPO<matrix, grp> mpo = make_mpo(lattice, model, parms);
        tim_model.end();
        maquis::cout << "Parsing model done!\n";
        

        boost::filesystem::path chkpfile(parms["chkpfile"].str());
        
        /// Initialize & load MPS
        tim_load.begin();
        int L = lattice.size();
        MPS<matrix, grp> mps;
        load(parms["chkpfile"].str(), mps);
        int _site;
        {
            alps::hdf5::archive ar(chkpfile / "props.h5");
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
        tim_load.end();
        maquis::cout << "Load MPS done!\n";
        maquis::cout << "Optimization at site " << site << " in " << lr << " direction." << std::endl;
        
        /// Canonize MPS
        mps.canonize(site);
        

        std::string boundary_name;
        construct_placements(mpo);
        
        /// Create TwoSite objects
        tim_ts_obj.begin();
        TwoSiteTensor<matrix, grp> tst(mps[site], mps[site+1]);
        MPSTensor<matrix, grp> ts_mps = tst.make_mps();
        MPOTensor<matrix, grp> ts_mpo = make_twosite_mpo<matrix,matrix>(mpo[site], mpo[site+1], mps[site].site_dim(), mps[site+1].site_dim(), true);
        if(lr == +1){
            ts_mpo.placement_l = mpo[site].placement_l;
            ts_mpo.placement_r = get_right_placement(ts_mpo, mpo[site].placement_l, mpo[site+1].placement_r);
        }else{
            ts_mpo.placement_l = get_left_placement(ts_mpo, mpo[site].placement_l, mpo[site+1].placement_r);
            ts_mpo.placement_r = mpo[site+1].placement_r;
        }
        tim_ts_obj.end();
        maquis::cout << "Two site obj done!\n";
        
        /// Compute left boundary
        tim_l_boundary.begin();
        Boundary<matrix, grp> left;
        boundary_name = "left" + boost::lexical_cast<std::string>(site) + ".h5";
        if ( exists(chkpfile / boundary_name) ) {
            maquis::cout << "Loading existing left boundary." << std::endl;
            storage::archive ar(chkpfile.string() +"/"+ boundary_name);
            ar["/tensor"] >> left;
        } else {
            left = mps.left_boundary();
            for (size_t i=0; i<site; ++i)
                left = contraction::overlap_mpo_left_step(mps[i], mps[i], left, mpo[i]);
        }
        #ifdef USE_AMBIENT
        if(exists(chkpfile / boundary_name) || lr == -1)
            for(size_t b = 0; b < left.aux_dim(); ++b){
                select_proc(ambient::scope::permute(b,ts_mpo.placement_l));
                storage::migrate(left[b]);
            }
        #endif
        tim_l_boundary.end();
        maquis::cout << "Left boundary done!\n";
        
        /// Compute right boundary
        tim_r_boundary.begin();
        Boundary<matrix, grp> right;
        boundary_name = "right" + boost::lexical_cast<std::string>(site+2) + ".h5";
        if ( exists(chkpfile / boundary_name) ) {
            maquis::cout << "Loading existing right boundary." << std::endl;
            storage::archive ar(chkpfile.string() +"/"+ boundary_name);
            ar["/tensor"] >> right;
        } else {
            right = mps.right_boundary();
            for (int i=L-1; i>site+1; --i)
                right = contraction::overlap_mpo_right_step(mps[i], mps[i], right, mpo[i]);
        }
        #ifdef USE_AMBIENT
        if(exists(chkpfile / boundary_name) || lr == +1) 
            for(size_t b = 0; b < right.aux_dim(); ++b){
                select_proc(ambient::scope::permute(b,ts_mpo.placement_r));
                storage::migrate(right[b]);
            }
        #endif
        tim_r_boundary.end();
        maquis::cout << "Right boundary done!\n";
        
        // Clearing unneeded MPS Tensors
        for (int k = 0; k < mps.length(); k++){
            if(k == site || k == site+1) continue;
            if(lr == -1 && site > 0   && k == site-1) continue; 
            if(lr == +1 && site < L-2 && k == site+2) continue; 
            mps[k].data().clear();
        }
        
        
        std::vector<MPSTensor<matrix, grp> > ortho_vecs;
        std::pair<double, MPSTensor<matrix, grp> > res;
        SiteProblem<matrix, grp> sp(left, right, ts_mpo);

        /// Optimization: JCD
        tim_optim_jcd.begin();
        res = solve_ietl_jcd(sp, ts_mps, parms, ortho_vecs);
        tst << res.second;
        maquis::cout.precision(10);
        maquis::cout << "Energy " << lr << " " << res.first << std::endl;
        tim_optim_jcd.end();
        maquis::cout << "Optim. JCD done!\n";
        
        double alpha = parms["alpha_main"];
        double cutoff = parms["truncation_final"];
        std::size_t Mmax = parms["max_bond_dimension"];
        
        /// Truncation of MPS
        tim_truncation.begin();
        truncation_results trunc;
        if (lr == +1)
        {
            boost::tie(mps[site], mps[site+1], trunc) = tst.split_mps_l2r(Mmax, cutoff);
            
            block_matrix<matrix, grp> t;
            t = mps[site+1].normalize_left(DefaultSolver());
            if (site+1 < L-1) mps[site+2].multiply_from_left(t);
        }
        if (lr == -1){
            boost::tie(mps[site], mps[site+1], trunc) = tst.split_mps_r2l(Mmax, cutoff);
            
            block_matrix<matrix, grp> t;
            t = mps[site].normalize_right(DefaultSolver());
            if (site > 0) mps[site-1].multiply_from_right(t);
        }
        tim_truncation.end();
        maquis::cout << "Truncation done!\n";

        /// Compute new boundary
        // TODO: optional here...
        
    } catch (std::exception & e) {
        maquis::cerr << "Exception caught:" << std::endl << e.what() << std::endl;
        exit(1);
    }
}

