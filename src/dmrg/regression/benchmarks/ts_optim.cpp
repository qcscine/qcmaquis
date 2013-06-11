/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Bela Bauer <bauerb@comp-phys.org>
 *
 *****************************************************************************/

#include <iostream>
#include <sstream>
#include <fstream>

#include <boost/shared_ptr.hpp>

#include <alps/hdf5.hpp>

#include "matrix_selector.hpp" /// define matrix
#include "symm_selector.hpp"   /// define grp

#include "dmrg/models/factory.h"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/twositetensor.h"

#include "dmrg/mp_tensors/ietl_lanczos_solver.h"
#include "dmrg/mp_tensors/ietl_jacobi_davidson.h"

#include "dmrg/utils/DmrgParameters2.h"

#include "utils/timings.h"


template<class Matrix, class SymmGroup>
struct SiteProblem
{
    SiteProblem(Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left_,
                Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right_,
                MPOTensor<Matrix, SymmGroup> const & mpo_)
    : left(left_)
    , right(right_)
    , mpo(mpo_)
    {
#ifdef AMBIENT
        mpo.persist();
#endif
    }
    
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & left;
    Boundary<typename storage::constrained<Matrix>::type, SymmGroup> const & right;
    MPOTensor<Matrix, SymmGroup> const & mpo;
    double ortho_shift;
};


int main(int argc, char ** argv)
{
    try {
        if (argc != 3)
            throw std::runtime_error("Usage: <parms> <model_parms>");
        
        /// Loading parameters
        std::ifstream param_file(argv[1]);
        if (!param_file)
            throw std::runtime_error("Could not open parameter file.");
        DmrgParameters parms(param_file);
        
        /// Loading model
        std::ifstream model_file(argv[2]);
        if (!model_file)
            throw std::runtime_error("Could not open model file.");
        ModelParameters model_parms(model_file);
        
        /// Timers
        Timer tim_model      ("Parsing model"),  tim_load       ("Load MPS");
        Timer tim_l_boundary ("Left boundary"),  tim_r_boundary ("Right boundary");
        Timer tim_optim_jcd  ("Optim. JCD"   ),  tim_truncation ("Truncation");
        Timer tim_ts_obj     ("Two site obj.");
        
        /// Parsing model
        tim_model.begin();
        boost::shared_ptr<Lattice> lattice;
        boost::shared_ptr<Model<matrix, grp> > model;
        model_parser<matrix, grp>(parms.get<std::string>("lattice_library"), parms.get<std::string>("model_library"), model_parms,
                                  lattice, model);
        
        Hamiltonian<matrix, grp> H = model->H();
        MPO<matrix, grp> mpo = make_mpo(lattice->size(), H);
        tim_model.end();
        
        
        /// Initialize & load MPS
        tim_load.begin();
        int L = lattice->size();
        MPS<matrix, grp> mps(L);
        int _site;
        {
            alps::hdf5::archive ar(parms.get<std::string>("chkpfile"));
            ar["/state"]       >> mps;
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
        for (int i=L-1; i>site+1; --i)
            right = contraction::overlap_mpo_right_step(mps[i], mps[i], right, mpo[i]);
        tim_r_boundary.end();
        
        /// Create TwoSite objects
        tim_ts_obj.begin();
        TwoSiteTensor<matrix, grp> tst(mps[site], mps[site+1]);
        MPSTensor<matrix, grp> ts_mps = tst.make_mps();
        MPOTensor<matrix, grp> ts_mpo = make_twosite_mpo<matrix,matrix>(mpo[site], mpo[site+1], mps[site].site_dim());
        tim_ts_obj.end();
        
        
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
        
        double alpha = parms.get<double>("alpha_main");
        double cutoff = parms.get<double>("truncation_final");
        std::size_t Mmax = parms.get<int>("max_bond_dimension");
        
        /// Truncation of MPS
        Logger iteration_log;
        tim_truncation.begin();
        if (lr == +1)
        {
            boost::tie(mps[site], mps[site+1]) = tst.split_mps_l2r(Mmax, cutoff, &iteration_log);
            
            block_matrix<matrix, grp> t;
            t = mps[site+1].normalize_left(DefaultSolver());
            if (site+1 < L-1) mps[site+2].multiply_from_left(t);
        }
        if (lr == -1){
            boost::tie(mps[site], mps[site+1]) = tst.split_mps_r2l(Mmax, cutoff, &iteration_log);
            
            block_matrix<matrix, grp> t;
            t = mps[site].normalize_right(DefaultSolver());
            if (site > 0) mps[site-1].multiply_from_right(t);
        }
        tim_truncation.end();

        /// Compute new boundary
        // TODO: optional here...
        
        
    } catch (std::exception & e) {
        maquis::cerr << "Exception caught:" << std::endl << e.what() << std::endl;
        exit(1);
    }
}

