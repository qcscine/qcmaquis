/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MP_TENSORS_MPO_CONTRACTOR_SS_H
#define MP_TENSORS_MPO_CONTRACTOR_SS_H

#include "dmrg/optimize/ietl_lanczos_solver.h"

#include "dmrg/optimize/ietl_jacobi_davidson.h"

#ifdef HAVE_ARPACK
#include "dmrg/optimize/arpackpp_solver.h"
#endif

#include "dmrg/mp_tensors/siteproblem.h"
#include "dmrg/utils/BaseParameters.h"

#define BEGIN_TIMING(name) \
now = std::chrono::high_resolution_clock::now();
#define END_TIMING(name) \
then = std::chrono::high_resolution_clock::now(); \
    maquis::cout << "Time elapsed in " << name << ": " << std::chrono::duration<double>(then-now).count() << std::endl;


/// TODO: 1) implement two-site time evolution. (single-site is stuck in initial MPS structure)
///       2) implement zip-up compression. E. M. Stoudenmire and S. R. White, New Journal of Physics 12, 055026 (2010).

template<class Matrix, class SymmGroup, class Storage>
class mpo_contractor_ss
{
    typedef contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup> contr;

public:
    mpo_contractor_ss(MPS<Matrix, SymmGroup> const & mps_,
                      MPO<Matrix, SymmGroup> const & mpo_,
                      BaseParameters & parms_)
    : mps(mps_)
    , mpsp(mps_)
    , mpo(mpo_)
    , parms(parms_)
    {
        mps.canonize(0);
        init_left_right(mpo);

        mpsp = mps;
    }

    std::pair<double,double> sweep(int sweep)
    {
        std::chrono::high_resolution_clock::time_point sweep_now = std::chrono::high_resolution_clock::now();

        std::size_t L = mps.length();

        std::pair<double,double> eps;
        block_matrix<Matrix, SymmGroup> norm_boudary;
        norm_boudary.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);

        for (int _site = 0; _site < 2*L; ++_site)
        {
            int site, lr;
            if (_site < L) {
                site = _site;
                lr = 1;
            } else {
                site = 2*L-_site-1;
                lr = -1;
            }

            SiteProblem<Matrix, SymmGroup> sp(mps[site], left_[site], right_[site+1], mpo[site]);
            ietl::mult(sp, mps[site], mpsp[site]);

            if (lr == +1) {
                if (site < L-1) {
                    block_matrix<Matrix, SymmGroup> t;
                    t = mpsp[site].leftNormalizeAndReturn(DefaultSolver());
                    mpsp[site+1].multiply_from_left(t);
                }

                left_[site+1] = contr::overlap_mpo_left_step(mpsp[site], mps[site], left_[site], mpo[site]);
                norm_boudary = contr::overlap_left_step(mpsp[site], MPSTensor<Matrix,SymmGroup>(mpsp[site]), norm_boudary);
            } else if (lr == -1) {
                if (site > 0) {
                    block_matrix<Matrix, SymmGroup> t;
                    t = mpsp[site].rightNormalizeAndReturn(DefaultSolver());
                    mpsp[site-1].multiply_from_right(t);
                }

                right_[site] = contr::overlap_mpo_right_step(mpsp[site], mps[site], right_[site+1], mpo[site]);
                norm_boudary = contr::overlap_right_step(mpsp[site], MPSTensor<Matrix,SymmGroup>(mpsp[site]), norm_boudary);
            }

            if (_site == L-1) {
                double nn = maquis::real( norm_boudary.trace() );
                eps.first = nn - 2.*maquis::real(left_[L][0].trace());

                /// prepare backward sweep
                norm_boudary = block_matrix<Matrix, SymmGroup>();
                norm_boudary.insert_block(Matrix(1, 1, 1), mps[L-1].col_dim()[0].first, mps[L-1].col_dim()[0].first);
            }

            if (_site == 2*L-1) {
                double nn = maquis::real( norm_boudary.trace() );
                eps.second = nn - 2.*maquis::real(right_[0][0].trace());
            }

        }

        return eps; /// note: the actual eps contain a constant, which is not important here.
    }

    void finalize()
    {
        mpsp[0].normalize_right(DefaultSolver());
    }

    MPS<Matrix, SymmGroup> get_original_mps() const { return mps; }
    MPS<Matrix, SymmGroup> get_current_mps() const { return mpsp; }

private:
    void init_left_right(MPO<Matrix, SymmGroup> const & mpo)
    {
        std::size_t L = mps.length();

        left_.resize(mpo.length()+1);
        right_.resize(mpo.length()+1);

        Storage::drop(left_[0]);
        left_[0] = mps.left_boundary();
        Storage::StoreToFile(left_[0]);

        for (int i = 0; i < L; ++i) {
            Storage::drop(left_[i+1]);
            left_[i+1] = contr::overlap_mpo_left_step(mpsp[i], mps[i], left_[i], mpo[i]);
            Storage::StoreToFile(left_[i+1]);
        }

        Storage::drop(right_[L]);
        right_[L] = mps.right_boundary();
        Storage::StoreToFile(right_[L]);

        for(int i = L-1; i >= 0; --i) {
            Storage::drop(right_[i]);
            right_[i] = contr::overlap_mpo_right_step(mpsp[i], mps[i], right_[i+1], mpo[i]);
            Storage::StoreToFile(right_[i]);
        }
    }

    MPS<Matrix, SymmGroup> mps, mpsp;
    MPO<Matrix, SymmGroup> const& mpo;

    BaseParameters & parms;
    std::vector<Boundary<typename storage::constrained<Matrix>::type, SymmGroup> > left_, right_;
};

#endif

