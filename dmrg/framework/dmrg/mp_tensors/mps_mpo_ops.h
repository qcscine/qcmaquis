/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPS_MPO_OPS_H
#define MPS_MPO_OPS_H

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/mp_tensors/special_mpos.h"
#include "dmrg/mp_tensors/contractions.h"

#include "dmrg/utils/utils.hpp"
#include "utils/traits.hpp"

template<class Matrix, class SymmGroup>
std::vector<Boundary<Matrix, SymmGroup> >
left_mpo_overlaps(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo)
{
    assert(mpo.length() == mps.length());
    std::size_t L = mps.length();
    
    std::vector<Boundary<Matrix, SymmGroup> > left_(L+1);
    left_[0] = mps.left_boundary();
    
    for (int i = 0; i < L; ++i) {
        left_[i+1] = contraction::overlap_mpo_left_step(mps[i], mps[i], left_[i], mpo[i]);
    }
    return left_;
}

template<class Matrix, class SymmGroup>
std::vector<Boundary<Matrix, SymmGroup> >
right_mpo_overlaps(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo)
{
    assert(mpo.length() == mps.length());
    std::size_t L = mps.length();
    
    std::vector<Boundary<Matrix, SymmGroup> > right_(L+1);
    right_[L] = mps.right_boundary();
    
    for (int i = L-1; i >= 0; --i) {
        right_[i] = contraction::overlap_mpo_right_step(mps[i], mps[i], right_[i+1], mpo[i]);
    }
    return right_;
}

template<class Matrix, class SymmGroup>
double expval(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo, int d)
{
    if (d == 0) {
        std::vector<Boundary<Matrix, SymmGroup> > left_ = left_mpo_overlaps(mps, mpo);
        assert( check_real(left_[mps.length()][0].trace()) );
        return maquis::real(left_[mps.length()][0].trace());
    } else {
        std::vector<Boundary<Matrix, SymmGroup> > right_ = right_mpo_overlaps(mps, mpo);
        assert( check_real(right_[0][0].trace()) );
        return maquis::real(right_[0][0].trace());
    }
}

template<class Matrix, class SymmGroup>
double expval(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo,
              bool verbose = false)
{
    assert(mpo.length() == mps.length());
    std::size_t L = mps.length();
    
    Boundary<Matrix, SymmGroup> left = mps.left_boundary();
    
    semi_parallel_for (locale::compact(L), locale i = 0; i < L; ++i) {
        if (verbose)
            maquis::cout << "expval site " << (size_t)i << std::endl;
        left = contraction::overlap_mpo_left_step(mps[i], mps[i], left, mpo[i]);
    }
    
    return maquis::real(left[0].trace());
}

template<class Matrix, class SymmGroup>
std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> multi_expval(MPS<Matrix, SymmGroup> const & mps,
                                                                       MPO<Matrix, SymmGroup> const & mpo)
{
    assert(mpo.length() == mps.length());
    std::size_t L = mps.length();
    
    Boundary<Matrix, SymmGroup> left = mps.left_boundary();
    
    for (int i = 0; i < L; ++i) {
        left = contraction::overlap_mpo_left_step(mps[i], mps[i], left, mpo[i]);
    }
    
    return left.traces();
}

template<class Matrix, class SymmGroup>
typename MPS<Matrix, SymmGroup>::scalar_type norm(MPS<Matrix, SymmGroup> const & mps)
{
    std::size_t L = mps.length();
    
    block_matrix<Matrix, SymmGroup> left;
    left.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
    
    semi_parallel_for (locale::compact(L), locale i = 0; i < L; ++i) {
        MPSTensor<Matrix, SymmGroup> cpy = mps[i];
        left = contraction::overlap_left_step(mps[i], cpy, left); // serial
    }
    
    return trace(left);
}

template<class Matrix, class SymmGroup>
typename MPS<Matrix, SymmGroup>::scalar_type overlap(MPS<Matrix, SymmGroup> const & mps1,
                                                     MPS<Matrix, SymmGroup> const & mps2)
{
    assert(mps1.length() == mps2.length());
    
    std::size_t L = mps1.length();
    
    block_matrix<Matrix, SymmGroup> left;
    left.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
    
    semi_parallel_for (locale::compact(L), locale i = 0; i < L; ++i) {
        left = contraction::overlap_left_step(mps1[i], mps2[i], left);
    }
    
    return trace(left);
}

template<class Matrix, class SymmGroup>
std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> multi_overlap(MPS<Matrix, SymmGroup> const & mps1,
                                                                        MPS<Matrix, SymmGroup> const & mps2)
{
    // assuming mps2 to have `correct` shape, i.e. left size=1, right size=1
    //          mps1 more generic, i.e. left size=1, right size arbitrary
    
    
    assert(mps1.length() == mps2.length());
    
    std::size_t L = mps1.length();
    
    block_matrix<Matrix, SymmGroup> left;
    left.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
    
    for (int i = 0; i < L; ++i) {
        left = contraction::overlap_left_step(mps1[i], mps2[i], left);
    }
    
    assert(left.right_basis().sum_of_sizes() == 1);
    std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals;
    vals.reserve(left.left_basis().sum_of_sizes());
    for (int n=0; n<left.n_blocks(); ++n)
        for (int i=0; i<left.left_basis()[n].second; ++i)
            vals.push_back( left[n](i,0) );
        
    return vals;
}

template<class Matrix, class SymmGroup>
std::vector<double>
calculate_bond_renyi_entropies(MPS<Matrix, SymmGroup> & mps, double n,
                               std::vector< std::vector<double> > * spectra = NULL) // to be optimized later
{
    std::size_t L = mps.length();
    std::vector<double> ret;
    
    MPS<Matrix, SymmGroup> const& constmps = mps;
    
    block_matrix<Matrix, SymmGroup> lb;
    
    if (spectra != NULL)
        spectra->clear();
    
    mps.canonize(0);
    for (std::size_t p = 1; p < L; ++p)
    {
        block_matrix<Matrix, SymmGroup> t, u, v;
        block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> s;
        
        constmps[p-1].make_left_paired();
        constmps[p].make_right_paired();
        
        gemm(constmps[p-1].data(), constmps[p].data(), t);
        
        svd(t, u, v, s);
        
        std::vector<double> sv = maquis::dmrg::detail::bond_renyi_entropies(s);
        
        if (spectra != NULL)
            spectra->push_back(sv);
        
        double S = 0;
        if (n == 1) {
            for (std::vector<double>::const_iterator it = sv.begin();
                 it != sv.end(); ++it)
                S += *it * log(*it);
            ret.push_back(-S);
        } else {
            for (std::vector<double>::const_iterator it = sv.begin();
                 it != sv.end(); ++it)
                S += pow(*it, n);
            ret.push_back(1/(1-n)*log(S));
        }
        
        mps.move_normalization_l2r(p-1, p, DefaultSolver());
    }
    
    return ret;
}

template<class Matrix, class SymmGroup>
std::vector<double>
calculate_bond_entropies(MPS<Matrix, SymmGroup> & mps,
                         std::vector< std::vector<double> > * spectra = NULL)
{
    return calculate_bond_renyi_entropies(mps, 1, spectra);
}

template<class Matrix, class SymmGroup>
typename MPS<Matrix, SymmGroup>::scalar_type dm_trace(MPS<Matrix, SymmGroup> const& mps, Index<SymmGroup> const& phys_psi)
{
    typedef typename SymmGroup::charge charge;
    charge I = SymmGroup::IdentityCharge;
    size_t L = mps.length();
    
    Index<SymmGroup> phys_rho = phys_psi * adjoin(phys_psi);
    ProductBasis<SymmGroup> pb(phys_psi, phys_psi, boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                                                       boost::lambda::_1, -boost::lambda::_2));
    
    Matrix identblock(phys_rho.size_of_block(I), 1, 0.);
    for (int s=0; s<phys_psi.size(); ++s)
        for (int ss=0; ss<phys_psi[s].second; ++ss) {
            identblock(pb(phys_psi[s].first, phys_psi[s].first) + ss*phys_psi[s].second+ss, 0) = 1.;
        }
    block_matrix<Matrix, SymmGroup> ident;
    ident.insert_block(identblock, I, I);
    
    Index<SymmGroup> trivial_i;
    trivial_i.insert(std::make_pair(I, 1));
    MPSTensor<Matrix, SymmGroup> mident(phys_rho, trivial_i, trivial_i);
    mident.data() = ident;
    
    MPS<Matrix,SymmGroup> mps_ident(L);
    for (int p=0; p<L; ++p)
        mps_ident[p] = mident;
    
    return overlap(mps, mps_ident);
}


// Specific to Fermi-Hubbard on a Ladder!!
template<class Matrix, class SymmGroup>
void fix_density(MPS<Matrix, SymmGroup> & mps, std::vector<block_matrix<Matrix, SymmGroup> > const & dens_ops, std::vector<std::vector<double> > const & dens)
{
    assert( mps.size() == dens[0].size() );
    assert( dens_ops.size() == dens.size() );
    size_t L = mps.size();
    
    mps.normalize_left();
    mps.canonize(0);
    for (int p=0; p<L; ++p)
    {
        

        Index<SymmGroup> phys = mps[p].site_dim();
        typename SymmGroup::charge empty, up, down, updown;
        empty[0]  = 0;  empty[1]  = 0;
        up[0]     = 1;  up[1]     = 0;
        down[0]   = 0;  down[1]   = 1;
        updown[0] = 1;  updown[1] = 1;
        
        
        block_matrix<Matrix, SymmGroup> rho = contraction::density_matrix(mps[p], mps[p]);
        
        for (size_t j=0; j<dens.size(); ++j) {
            
            MPSTensor<Matrix, SymmGroup> tmp = contraction::local_op(mps[p], dens_ops[j]);
            double cur_dens = mps[p].scalar_overlap(tmp);
            maquis::cout << "Density[" << j << "] (before) = " << cur_dens << std::endl;
        }
        
        double a = trace(rho(down, down)) * trace(rho(updown, updown));
        double b = trace(rho(up, up)) * trace(rho(down, down)) + dens[0][p] * trace(rho(updown, updown)) - dens[1][p] * trace(rho(updown, updown));
        double c = - dens[1][p] * trace(rho(up, up));
        double k2 = ( -b + sqrt(b*b - 4*a*c) ) / (2*a);
        
        double k1 = dens[0][p] / ( trace(rho(up, up)) + k2*trace(rho(updown,updown)) );
        
        double t0 = 0.;
        t0 += k1*trace( rho(up, up) );
        t0 += k2*trace( rho(down, down) );
        t0 += k1*k2*trace( rho(updown, updown) );
        double k0 = (1.-t0) / trace(rho(empty, empty));
        
        maquis::cout << "k0 = " << k0 << std::endl;
        maquis::cout << "k1 = " << k1 << std::endl;
        maquis::cout << "k2 = " << k2 << std::endl;
        assert( k0 > 0 ); // not always the case!!!
        
        block_matrix<Matrix, SymmGroup> rescale = identity_matrix<Matrix>(phys);
        rescale(empty, empty) *= std::sqrt(k0);
        rescale(up, up) *= std::sqrt(k1);
        rescale(down, down) *= std::sqrt(k2);
        rescale(updown, updown) *= std::sqrt(k1*k2);
        
        mps[p] = contraction::local_op(mps[p], rescale);
        
        {
            for (size_t j=0; j<dens.size(); ++j) {
                MPSTensor<Matrix, SymmGroup> tmp = contraction::local_op(mps[p], dens_ops[j]);
                double meas_dens = mps[p].scalar_overlap(tmp) / mps[p].scalar_norm();
                maquis::cout << "Density[" << j << "] (after) = " << meas_dens << ", should be " << dens[j][p] << std::endl;
            }
        }
        
        block_matrix<Matrix, SymmGroup> t_norm = mps[p].normalize_left(DefaultSolver());
        if (p < L-1)
            mps[p+1].multiply_from_left(t_norm);
        
    }

}


#endif
