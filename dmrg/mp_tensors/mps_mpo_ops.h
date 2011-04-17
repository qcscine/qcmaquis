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

#include "mp_tensors/mps.h"
#include "mp_tensors/mpo.h"

#include "mp_tensors/special_mpos.h"

#include "mp_tensors/contractions.h"

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
//        zout << "Left at " << i+1 << " " << left.data_[0] << endl;
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
//        zout << "right at " << i << " " << right.data_[0] << endl;
    }
    return right_;
}

template<class Matrix, class SymmGroup>
double expval(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo, int d)
{
    if (d == 0) {
        std::vector<Boundary<Matrix, SymmGroup> > left_ = left_mpo_overlaps(mps, mpo);
        return left_[mps.length()].traces()[0];
    } else {
        std::vector<Boundary<Matrix, SymmGroup> > right_ = right_mpo_overlaps(mps, mpo);
        return right_[0].traces()[0];
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
            cout << "expval site " << i << endl;
        MPSTensor<Matrix, SymmGroup> bkp = mps[i];
        left = contraction::overlap_mpo_left_step(mps[i], bkp, left, mpo[i]);
    }
    
    std::vector<typename Matrix::value_type> traces = left.traces();
    
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
    
    return left.traces();
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type norm(MPS<Matrix, SymmGroup> const & mps)
{
    std::size_t L = mps.length();
    
    block_matrix<Matrix, SymmGroup> left;
    left.insert_block(Matrix(1, 1, 1), SymmGroup::SingletCharge, SymmGroup::SingletCharge);
    
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
    left.insert_block(Matrix(1, 1, 1), SymmGroup::SingletCharge, SymmGroup::SingletCharge);
    
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
        block_matrix<blas::diagonal_matrix<double>, SymmGroup> s;
        
        mps[p-1].make_left_paired();
        mps[p].make_right_paired();
        
        gemm(mps[p-1].data(), mps[p].data(), t);
        
        svd(t, u, v, s);
        
        std::vector<double> sv;
        
        double r = 0;
        for (std::size_t k = 0; k < s.n_blocks(); ++k)
            for (typename blas::diagonal_matrix<double>::element_iterator it = elements(s[k]).first;
                 it != elements(s[k]).second; ++it)
            {
                double a = fabs(*it);
                if (a > 1e-10)
                    sv.push_back(a*a);
            }
        
//        cout << p << " " << sv[0] << " " << sv[1] << endl;
        
        r = std::accumulate(sv.begin(), sv.end(), double(0));
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

#endif
