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

#include <mpi.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/shared_ptr.hpp>
#include <alps/hdf5.hpp>

#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> amatrix;

#ifdef USE_AMBIENT
#include "dmrg/block_matrix/detail/ambient.hpp"
typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > pmatrix;
typedef pmatrix matrix;
#else
typedef amatrix matrix;
#endif

#include "dmrg/block_matrix/symmetry.h"
#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_U1)
typedef U1 grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#else
#error "No valid symmetry defined."
#endif

#include "dmrg/models/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/twositetensor.h"

#include "dmrg/optimize/ietl_lanczos_solver.h"
#include "dmrg/optimize/ietl_jacobi_davidson.h"

#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"

#include "utils/timings.h"


#include <boost/mpi.hpp>
namespace mpi = boost::mpi;


enum optim_normalize_t {with_normalization, no_normalization};

template<class Matrix, class SymmGroup>
struct SiteProblem
{
    SiteProblem(Boundary<Matrix, SymmGroup> const & left_,
                Boundary<Matrix, SymmGroup> const & right_,
                MPOTensor<Matrix, SymmGroup> const & mpo_)
    : left(left_)
    , right(right_)
    , mpo(mpo_)
    { }
    
    Boundary<Matrix, SymmGroup> const & left;
    Boundary<Matrix, SymmGroup> const & right;
    MPOTensor<Matrix, SymmGroup> const & mpo;
    double ortho_shift;
};

template<class Matrix, class DiagMatrix, class SymmGroup>
void split_ts(TwoSiteTensor<Matrix, SymmGroup> const& tst, 
              MPSTensor<Matrix, SymmGroup> & mps1, 
              MPSTensor<Matrix, SymmGroup> & mps2, 
              block_matrix<DiagMatrix, SymmGroup> & inv, 
              double cutoff, std::size_t Mmax)
{
    tst.make_both_paired();
    
    typedef typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type dmt;
    block_matrix<Matrix, SymmGroup> u, v, m1, m2;
    block_matrix<dmt, SymmGroup> s;
    
    truncation_results trunc = svd_truncate(tst.data(), u, v, s, cutoff, Mmax, true);
    
    gemm(u, s, m1);
    gemm(s, v, m2);
    
    mps1 = MPSTensor<Matrix,SymmGroup>(tst.local_site_dim(0), tst.row_dim(), m1.right_basis(), m1, LeftPaired);
    mps2 = MPSTensor<Matrix,SymmGroup>(tst.local_site_dim(1), m2.left_basis(), tst.col_dim(), m2, RightPaired);
    
    inv = s;
    typename dmt::diagonal_iterator it, end;
    for (size_t k=0; k<inv.n_blocks(); ++k)
        for (boost::tie(it, end) = inv[k].diagonal(); it != end; ++it)
            *it = 1. / *it;
}

template<class Matrix, class DiagMatrix, class SymmGroup>
void central_optim(unsigned site, block_matrix<DiagMatrix, SymmGroup> & inv,
                   MPS<Matrix, SymmGroup> & mps, MPO<Matrix, SymmGroup> const& mpo,
                   Boundary<Matrix, SymmGroup> const & left, Boundary<Matrix, SymmGroup> const & right,
                   Boundary<Matrix, SymmGroup> & newleft, Boundary<Matrix, SymmGroup> & newright,
                   double cutoff, std::size_t Mmax, BaseParameters & parms)
{
    
    /// Create TwoSite objects
    mps[site].multiply_from_right(inv);
    TwoSiteTensor<Matrix, SymmGroup> tst(mps[site], mps[site+1]);
    MPSTensor<Matrix, SymmGroup> ts_mps = tst.make_mps();
    MPOTensor<Matrix, SymmGroup> ts_mpo = make_twosite_mpo<Matrix,Matrix>(mpo[site], mpo[site+1], mps[site].site_dim(), mps[site+1].site_dim());
    maquis::cout << "Two site obj done!\n";
    
    
    /// Create site problem
    std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs;
    std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
    SiteProblem<Matrix, SymmGroup> sp(left, right, ts_mpo);
    
    /// Optimization: JCD
    res = solve_ietl_jcd(sp, ts_mps, parms, ortho_vecs);
    tst << res.second;
    maquis::cout.precision(10);
    maquis::cout << "Energy M " << res.first << std::endl;
    maquis::cout << "Optim. JCD done!\n";
    
    /// Truncation of MPS
    truncation_results trunc;
    split_ts(tst, mps[site], mps[site+1], inv, cutoff, Mmax);
    maquis::cout << "Truncation done!\n";
    
    /// Move boundaries
    {
        MPSTensor<Matrix, SymmGroup> t = mps[site];
        t.multiply_from_right(inv);
        newleft = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(t, t, left, mpo[site]);
    }
    {
        MPSTensor<Matrix, SymmGroup> t = mps[site+1];
        t.multiply_from_left(inv);
        newright = contraction::Engine<matrix, matrix, grp>::overlap_mpo_right_step(t, t, right, mpo[site+1]);
    }
}


template<class Matrix, class SymmGroup>
void dmrg_optim(unsigned site, unsigned local_site, int lr, int L,
                MPS<Matrix, SymmGroup> & mps, MPO<Matrix, SymmGroup> const& mpo,
                std::vector<Boundary<Matrix, SymmGroup> > & left, std::vector<Boundary<Matrix, SymmGroup> > & right,
                double cutoff, std::size_t Mmax, BaseParameters & parms,
                optim_normalize_t normalize)
{
    /// Create TwoSite objects
    TwoSiteTensor<Matrix, SymmGroup> tst(mps[site], mps[site+1]);
    MPSTensor<Matrix, SymmGroup> ts_mps = tst.make_mps();
    MPOTensor<Matrix, SymmGroup> ts_mpo = make_twosite_mpo<Matrix,Matrix>(mpo[site], mpo[site+1], mps[site].site_dim(), mps[site+1].site_dim());
    maquis::cout << "Two site obj done!\n";
    
    
    /// Create site problem
    std::vector<MPSTensor<Matrix, SymmGroup> > ortho_vecs;
    std::pair<double, MPSTensor<Matrix, SymmGroup> > res;
    SiteProblem<Matrix, SymmGroup> sp(left[local_site], right[local_site+1], ts_mpo);
    
    /// Optimization: JCD
    res = solve_ietl_jcd(sp, ts_mps, parms, ortho_vecs);
    tst << res.second;
    maquis::cout.precision(10);
    maquis::cout << "Energy " << lr << " " << res.first << std::endl;
    maquis::cout << "Optim. JCD done!\n";
    
    /// Truncation of MPS
    truncation_results trunc;
    if (lr == +1){
        boost::tie(mps[site], mps[site+1], trunc) = tst.split_mps_l2r(Mmax, cutoff);
        
        if (normalize == with_normalization) { // site != L/2-1
            block_matrix<Matrix, SymmGroup> t;
            t = mps[site+1].normalize_left(DefaultSolver());
            if (site+1 < L-1) mps[site+2].multiply_from_left(t);
        }
    }
    if (lr == -1){
        boost::tie(mps[site], mps[site+1], trunc) = tst.split_mps_r2l(Mmax, cutoff);

        if (normalize == with_normalization) { //site != L/2+1
            block_matrix<Matrix, SymmGroup> t;
            t = mps[site].normalize_right(DefaultSolver());
            if (site > 0) mps[site-1].multiply_from_right(t);
        }
    }
    maquis::cout << "Truncation done!\n";

    /// Move boundaries
    if (lr == +1) {
        if (site < L+1) left[local_site+1] = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(mps[site], mps[site], left[local_site], mpo[site]);
    }
    if (lr == -1) {
        if (site > 0) right[local_site] = contraction::Engine<matrix, matrix, grp>::overlap_mpo_right_step(mps[site+1], mps[site+1], right[local_site+1], mpo[site+1]);
    }
    
}

int main(int argc, char ** argv)
{
    try {
        mpi::environment env(argc, argv);
        mpi::communicator comm;
        
        DmrgOptions opt(argc, argv);
        if (!opt.valid) return 0;
        DmrgParameters parms = opt.parms;
        
        maquis::cout.precision(10);
    
        /// Parsing model
        Lattice lattice(parms);
        Model<matrix, grp> model(lattice, parms);
        
        MPO<matrix, grp> mpo = make_mpo(lattice, model);
        maquis::cout << "Parsing model done!\n";
        
        /// Initialize & load MPS
        int L = lattice.size(), lr, site;
        MPS<matrix, grp> mps;
        
        unsigned nprocs = comm.size();
        unsigned rank   = comm.rank();
        unsigned chunk_size = L / nprocs;
        unsigned start = chunk_size * rank;
        unsigned end = (rank == nprocs-1) ? L : chunk_size * (rank+1);
        unsigned local_size = end - start;
        unsigned remainder = L - chunk_size * nprocs;
        lr = (rank % 2 == 0) ? +1 : -1;
        
        maquis::cout << "[" << rank << "] Task is " << start << " --> " << end << std::endl;
        
        if (nprocs % 2 != 0) throw std::runtime_error("number of processes must be a multiple of 2");
        
        if (comm.rank() == 0)
            mps = MPS<matrix, grp>(lattice.size(), *(model.initializer(lattice, parms)));
        mpi::broadcast(comm, mps, 0);
        maquis::cout << "[" << rank << "] " << "BCast MPS done!" << std::endl;;
        
        mps.canonize((lr == +1) ? start : end-1);
        
        std::vector<Boundary<matrix, grp> > left(local_size+1), right(local_size+1);
        {
            /// Compute left boundary
            std::vector<Boundary<matrix, grp> > tmp_left(L);
            if (rank == nprocs-1) {
                tmp_left[0] = mps.left_boundary();
                for (size_t i=0; i<L-1; ++i)
                    tmp_left[i+1] = contraction::Engine<matrix, matrix, grp>::overlap_mpo_left_step(mps[i], mps[i], tmp_left[i], mpo[i]);
                maquis::cout << "Left boundary done!\n";
            }
            
            /// Compute right boundary
            std::vector<Boundary<matrix, grp> > tmp_right(L);
            if (rank == 0) {
                tmp_right[L-1] = mps.right_boundary();
                for (int i=L-2; i>=0; --i)
                    tmp_right[i] = contraction::Engine<matrix, matrix, grp>::overlap_mpo_right_step(mps[i+1], mps[i+1], tmp_right[i+1], mpo[i+1]);
                maquis::cout << "Right boundary done!\n";
            }
            
            mpi::scatter(comm, tmp_left,  &left[0],  chunk_size, nprocs-1);
            mpi::scatter(comm, tmp_right, &right[0], chunk_size, 0);
            
            if (remainder > 0 && rank == 0) {
                std::cout << "!!0 doing something here!!" << std::endl;
                comm.send(nprocs-1, 0, &tmp_right[chunk_size*nprocs], remainder);
            }
            if (remainder > 0 && rank == nprocs-1) {
                std::cout << "!!1 doing something here!!" << std::endl;
                comm.recv(0, 0, &right[chunk_size], remainder);
                std::copy(tmp_left.begin() + chunk_size*nprocs, tmp_left.end(), left.begin() + chunk_size);
            }
        }
        
        mpi::request req[4];
        if (rank > 0) {
            req[0] = comm.isend(rank-1, 50, left[0]);
            req[1] = comm.isend(rank-1, 51, right[0]);
        }
        if (rank < nprocs-1) {
            req[2] = comm.irecv(rank+1, 50, left[local_size]);
            req[3] = comm.irecv(rank+1, 51, right[local_size]);
        }
        wait_all(req, req+4);
        
        double      cutoff = parms["truncation_final"];
        std::size_t Mmax   = parms["max_bond_dimension"];
        int nsweeps = parms["nsweeps"];
        
        typedef alps::numeric::associated_real_diagonal_matrix<matrix>::type dmt;
        block_matrix<dmt, grp> inv = identity_matrix<block_matrix<dmt,grp> >(mps[end-1].col_dim());
        for (int sweep = 0; sweep < nsweeps; ++sweep) {
            
            if (rank > 0 && rank < nprocs-1) {
                if (rank % 2 != 0) {
                    comm.recv(rank+1, 0, mps[end]);
                    comm.recv(rank+1, 1, right[local_size]);
                    
                    maquis::cout << "[" << rank << "] " << "doing : " << end-1 << " dir M" << std::endl;
                    central_optim(end-1, inv, mps, mpo, left[local_size-1], right[local_size], left[local_size], right[local_size-1], cutoff, Mmax, parms);
                    
                    comm.send(rank+1, 2, mps[end]);
                    comm.send(rank+1, 3, left[local_size]);
                } else {
                    comm.send(rank-1, 0, mps[start]);
                    comm.send(rank-1, 1, right[0]);
                    
                    comm.recv(rank-1, 2, mps[start]);
                    comm.recv(rank-1, 3, left[0]);
                }
            }
            
            if (rank % 2 != 0) {
                for (site = end-2; site >= start; --site) {
                    maquis::cout << "[" << rank << "] " << "doing : " << site << " dir " << -1 << std::endl;
                    optim_normalize_t normalize = (site != start) ? with_normalization : no_normalization;
                    dmrg_optim(site, site - start, -1, L, mps, mpo, left, right, cutoff, Mmax, parms, normalize);
                }
            } else {
                for (site = start; site < end-1; ++site) {
                    maquis::cout << "[" << rank << "] " << "doing : " << site << " dir " << +1 << std::endl;
                    optim_normalize_t normalize = (site != end-2) ? with_normalization : no_normalization;
                    dmrg_optim(site, site - start, +1, L, mps, mpo, left, right, cutoff, Mmax, parms, normalize);
                }
            }
            
            
            if (rank % 2 == 0) {
                comm.recv(rank+1, 0, mps[end]);
                comm.recv(rank+1, 1, right[local_size]);
                
                maquis::cout << "[" << rank << "] " << "doing : " << end-1 << " dir M" << std::endl;
                central_optim(end-1, inv, mps, mpo, left[local_size-1], right[local_size], left[local_size], right[local_size-1], cutoff, Mmax, parms);
                
                comm.send(rank+1, 2, mps[end]);
                comm.send(rank+1, 3, left[local_size]);
            } else {
                comm.send(rank-1, 0, mps[start]);
                comm.send(rank-1, 1, right[0]);
                
                comm.recv(rank-1, 2, mps[start]);
                comm.recv(rank-1, 3, left[0]);
            }
            
            
            if (rank % 2 != 0) {
                for (site = start; site < end-1; ++site) {
                    maquis::cout << "[" << rank << "] " << "doing : " << site << " dir " << +1 << std::endl;
                    optim_normalize_t normalize;
                    if (rank == nprocs-1) normalize = with_normalization;
                    else                  normalize = (site != end-2) ? with_normalization : no_normalization;
                    dmrg_optim(site, site - start, +1, L, mps, mpo, left, right, cutoff, Mmax, parms, normalize);
                }
            } else {
                for (site = end-2; site >= int(start); --site) {
                    maquis::cout << "[" << rank << "] " << "doing : " << site << " dir " << -1 << std::endl;
                    optim_normalize_t normalize;
                    if (rank == 0) normalize = with_normalization;
                    else           normalize = (site != start) ? with_normalization : no_normalization;
                    dmrg_optim(site, site - start, -1, L, mps, mpo, left, right, cutoff, Mmax, parms, normalize);
                }
            }
            
        }
        
        /// final merge
        /*
        if (rank > 0 && rank < nprocs-1) {
            if (rank % 2 != 0) {
                comm.recv(rank+1, 0, mps[end]);
                comm.recv(rank+1, 1, right[local_size]);
                
                maquis::cout << "[" << rank << "] " << "doing : " << end-1 << " dir M" << std::endl;
                central_optim(end-1, inv, mps, mpo, left[local_size-1], right[local_size], left[local_size], right[local_size-1], cutoff, Mmax, parms);
                
                comm.send(rank+1, 2, mps[end]);
                // comm.send(rank+1, 3, left[local_size]); // not really needed
            } else {
                comm.send(rank-1, 0, mps[start]);
                comm.send(rank-1, 1, right[0]);
                
                comm.recv(rank-1, 2, mps[start]);
                // comm.recv(rank-1, 3, left[0]); // not really needed
            }
        }
        */
        
        MPS<matrix, grp> full_mps(L);
        mpi::gather(comm, &mps[start], chunk_size, &full_mps[0], 0);
        if (remainder > 0 && rank == 0       ) comm.recv(nprocs-1, 100, &full_mps[chunk_size*nprocs], remainder);
        if (remainder > 0 && rank == nprocs-1) comm.send(0, 100, &mps[chunk_size], remainder);
        
        std::vector<block_matrix<dmt, grp> > all_inv(nprocs);
        mpi::gather(comm, inv, &all_inv[0], 0);
        
        if (rank == 0) {
            for (unsigned n=0; n<nprocs-1; ++n)
                full_mps[chunk_size*(n+1) - 1].multiply_from_right(all_inv[n]);

            double nn = norm(full_mps);
            maquis::cout << "Norm " << nn << std::endl;
            double energy = expval(full_mps, mpo);
            maquis::cout << "Energy final " << energy << std::endl;
            maquis::cout << "Energy/norm final " << energy/nn << std::endl;
            
            full_mps.canonize(0);
            nn = norm(full_mps);
            maquis::cout << "Norm " << nn << std::endl;
            energy = expval(full_mps, mpo);
            maquis::cout << "Energy final " << energy << std::endl;
            maquis::cout << "Energy/norm final " << energy/nn << std::endl;

            
            std::string chkpfile = parms["chkpfile"];
            /// save state to chkp dir
            if (!boost::filesystem::exists(chkpfile))
                boost::filesystem::create_directory(chkpfile);
            save(chkpfile, full_mps);
            
            /// save status
            {
                storage::archive ar(chkpfile+"/props.h5", "w");
                ar["/status/sweep"] << (parms["nsweeps"]-1);
                ar["/status/site"]  << -1;
            }
            
        }
    } catch (std::exception & e) {
        maquis::cerr << "Exception caught:" << std::endl << e.what() << std::endl;
        throw;
    }
}

