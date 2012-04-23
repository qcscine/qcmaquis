/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPS_MPO_OPS_H
#define MPS_MPO_OPS_H

#include <alps/numeric/real.hpp>

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/mp_tensors/special_mpos.h"
#include "dmrg/mp_tensors/contractions.h"

#include "dmrg/utils/utils.hpp"
//#include "dmrg/detail/algorithms_impl.h"

template<class Matrix, class SymmGroup>
std::vector<Boundary<Matrix, SymmGroup> >
left_mpo_overlaps(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo)
{
    assert(mpo.length() == mps.length());
    std::size_t L = mps.length();
    
    std::vector<Boundary<Matrix, SymmGroup> > left_(L+1);
    Boundary<Matrix, SymmGroup> left = mps.left_boundary();
    left_[0] = left;
    
    for (int i = 0; i < L; ++i) {
        MPSTensor<Matrix, SymmGroup> bkp = mps[i];
        left = contraction::overlap_mpo_left_step(mps[i], bkp, left, mpo[i]);
        left_[i+1] = left;
//        cout << "Left at " << i+1 << " " << left.data_[0] << endl;
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
    Boundary<Matrix, SymmGroup> right = mps.right_boundary();
    right_[L] = right;
    
    for (int i = L-1; i >= 0; --i) {
        MPSTensor<Matrix, SymmGroup> bkp = mps[i];
        right = contraction::overlap_mpo_right_step(mps[i], bkp, right, mpo[i]);
        right_[i] = right;
//        cout << "right at " << i << " " << right.data_[0] << endl;
    }
    return right_;
}

template<class Matrix, class SymmGroup>
double expval(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo, int d)
{
    if (d == 0) {
        std::vector<Boundary<Matrix, SymmGroup> > left_ = left_mpo_overlaps(mps, mpo);
        assert( check_real(left_[mps.length()].traces()[0]) );
        return alps::numeric::real(left_[mps.length()].traces()[0]);
    } else {
        std::vector<Boundary<Matrix, SymmGroup> > right_ = right_mpo_overlaps(mps, mpo);
        assert( check_real(right_[0].traces()[0]) );
        return alps::numeric::real(right_[0].traces()[0]);
    }
}

template<class Matrix, class SymmGroup>
double expval(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo,
              bool verbose = false)
{
    assert(mpo.length() == mps.length());
    std::size_t L = mps.length();
    
    Boundary<Matrix, SymmGroup> left = mps.left_boundary();
    
    for (int i = 0; i < L; ++i) {
        if (verbose)
            cout << "expval site " << i << std::endl;
        MPSTensor<Matrix, SymmGroup> bkp = mps[i];
        left = contraction::overlap_mpo_left_step(mps[i], bkp, left, mpo[i]);
    }
    
    std::vector<typename Matrix::value_type> traces = left.traces();
    assert( check_real(traces[0]) );
    return alps::numeric::real(traces[0]);
}

template<class Matrix, class SymmGroup>
std::vector<double> multi_expval(MPS<Matrix, SymmGroup> const & mps,
                                 MPO<Matrix, SymmGroup> const & mpo)
{
    assert(mpo.length() == mps.length());
    std::size_t L = mps.length();
    
    Boundary<Matrix, SymmGroup> left = mps.left_boundary();
    
    for (int i = 0; i < L; ++i) {
        MPSTensor<Matrix, SymmGroup> bkp = mps[i];
        left = contraction::overlap_mpo_left_step(mps[i], bkp, left, mpo[i]);
    }
    
    return alps::numeric::real(left.traces());
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type norm(MPS<Matrix, SymmGroup> const & mps)
{
    std::size_t L = mps.length();
    
    block_matrix<Matrix, SymmGroup> left;
    left.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
    
    for (int i = 0; i < L; ++i) {
        MPSTensor<Matrix, SymmGroup> cpy = mps[i];
        left = contraction::overlap_left_step(mps[i], cpy, left);
    }
    
    return trace(left);
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type overlap(MPS<Matrix, SymmGroup> const & mps1,
                                    MPS<Matrix, SymmGroup> const & mps2)
{
    assert(mps1.length() == mps2.length());
    
    std::size_t L = mps1.length();
    
    block_matrix<Matrix, SymmGroup> left;
    left.insert_block(Matrix(1, 1, 1), SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
    
    for (int i = 0; i < L; ++i) {
        left = contraction::overlap_left_step(mps1[i], mps2[i], left);
    }
    
    return trace(left);
}

template<class Matrix, class SymmGroup>
std::vector<double>
calculate_bond_renyi_entropies(MPS<Matrix, SymmGroup> & mps, double n)
{
    std::size_t L = mps.length();
    std::vector<double> ret;
    
    mps.normalize_right();
//    mps.canonize(1);
    
    block_matrix<Matrix, SymmGroup> lb;
    
    for (std::size_t p = 1; p < L; ++p)
    {
        block_matrix<Matrix, SymmGroup> t, u, v;
        block_matrix<typename maquis::types::associated_diagonal_matrix<Matrix>::type, SymmGroup> s;
        
        mps[p-1].make_left_paired();
        mps[p].make_right_paired();
        
        gemm(mps[p-1].data(), mps[p].data(), t);
        
        svd(t, u, v, s);
        
        std::vector<double> sv;
        
        for (std::size_t k = 0; k < s.n_blocks(); ++k)
            detail::iterable_matrix_impl<Matrix,SymmGroup>::caculate_bond_renyi_entropies_impl(s[k],sv);
//        cout << p << " " << sv[0] << " " << sv[1] << endl;
        
        double r = std::accumulate(sv.begin(), sv.end(), double(0));
//        std::transform(sv.begin(), sv.end(), sv.begin(),
//                       boost::lambda::_1 / r);
        
//        cout << r << " " << sv.size() << endl;
//        if (fabs(1-r) < 0.01)
//            std::copy(sv.begin(), sv.end(), std::ostream_iterator<double>(cout, " ")); cout << endl;
        
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
        
//        cout << ret.back() << endl;
        
        t = mps[p-1].normalize_left(SVD);
        mps[p].multiply_from_left(t);
    }
    
    return ret;
}

template<class Matrix, class SymmGroup>
std::vector<double>
calculate_bond_entropies(MPS<Matrix, SymmGroup> & mps)
{
    return calculate_bond_renyi_entropies(mps, 1);
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
            cout << "Density[" << j << "] (before) = " << cur_dens << std::endl;
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
        
        cout << "k0 = " << k0 << std::endl;
        cout << "k1 = " << k1 << std::endl;
        cout << "k2 = " << k2 << std::endl;
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
                cout << "Density[" << j << "] (after) = " << meas_dens << ", should be " << dens[j][p] << std::endl;
            }
        }
        
        block_matrix<Matrix, SymmGroup> t_norm = mps[p].normalize_left(SVD);
        if (p < L-1)
            mps[p+1].multiply_from_left(t_norm);
        
    }

}


#endif
