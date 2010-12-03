#include <cmath>
#include <iterator>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

#include "general_matrix.hpp"
#include "matrix_interface.hpp"
#include "resizable_matrix_interface.hpp"
#include "general_matrix_algorithms.h"
#include "matrix_algorithms.hpp"
typedef blas::general_matrix<double> Matrix;

#include "indexing.h"
#include "mpstensor.h"
#include "mpotensor.h"
#include "contractions.h"

#include "special_mpos.h"

typedef NullGroup grp;

struct MPS
{
    MPS(int L_, int Mmax)
    : mps_(L_), L(L_)
    {
        std::vector<int> bond_sizes(L+1, Mmax);
        for (std::size_t k = 0; k < L+1; ++k) {
            bond_sizes[k] = std::min(bond_sizes[k], (int)powf(2, k));
            bond_sizes[k] = std::min(bond_sizes[k], (int)powf(2, L-k));
        }
        std::copy(bond_sizes.begin(), bond_sizes.end(), std::ostream_iterator<int>(cout, " ")); cout << endl;
        
        Index<grp> phys; phys.insert(std::make_pair(NullGroup::Plus, 2));
        
        for (int i = 0; i < L; ++i)
        {
            Index<grp> li; li.insert(std::make_pair(NullGroup::Plus, bond_sizes[i]));
            Index<grp> ri; ri.insert(std::make_pair(NullGroup::Plus, bond_sizes[i+1]));
            
            mps_[i] = MPSTensor<Matrix, grp>(phys, li, ri);
        }
        
        for (int i = 0; i < L; ++i)
            mps_[i].normalize_left(SVD);
//        for (int i = 0; i < L; ++i)
//            mps_[i].normalize_right(SVD);
    }
    
    block_matrix<Matrix, grp> start_mtx()
    {
        Index<grp> i; i.insert(std::make_pair(NullGroup::Plus, 1));
        block_matrix<Matrix, grp> ret(i, i);
        ret.fill(functors::constant<double>(1));
        return ret;
    }
    
    std::vector<block_matrix<Matrix, grp> > left_mpo_overlaps(MPOTensor<Matrix, grp> & mpo)
    {
        std::vector<block_matrix<Matrix, grp> > left_(L+1);
        block_matrix<Matrix, grp> left = start_mtx();
        left_[0] = left;
        for (int i = 0; i < L; ++i) {
            MPSTensor<Matrix, grp> bkp = mps_[i];
            left = contraction::overlap_mpo_left_step(mps_[i], bkp, left, mpo);
            left_[i+1] = left;
        }
        return left_;
    }
    
    std::vector<block_matrix<Matrix, grp> > right_mpo_overlaps(MPOTensor<Matrix, grp> & mpo)
    {
        std::vector<block_matrix<Matrix, grp> > right_(L+1);
        block_matrix<Matrix, grp> right = start_mtx();
        right_[L] = right;
        
        for (int i = L-1; i >= 0; --i) {
            MPSTensor<Matrix, grp> bkp = mps_[i];
            right = contraction::overlap_mpo_right_step(mps_[i], bkp, right, mpo);
            right_[i] = right;
        }
        return right_;
    }  
    
    double norm(MPOTensor<Matrix, grp> * mpo = NULL)
    {
        MPOTensor<Matrix, grp> id_mpo;
        if (mpo != NULL)
            id_mpo = *mpo;
        else
            id_mpo = identity_mpo<Matrix>(mps_[0].site_dim());
        
        std::vector<block_matrix<Matrix, grp> > left_ = left_mpo_overlaps(id_mpo);
        return trace(*left_.rbegin());
    }
    
    double stupid_site_expval(MPOTensor<Matrix, grp> & mpo, int w)
    {
        MPOTensor<Matrix, grp> id_mpo = identity_mpo<Matrix>(mps_[0].site_dim());
        
        std::vector<block_matrix<Matrix, grp> >
        left_ = left_mpo_overlaps(id_mpo),
        right_ = right_mpo_overlaps(id_mpo);
        
        MPSTensor<Matrix, grp> bkp = mps_[w];
        block_matrix<Matrix, grp> t = contraction::overlap_mpo_left_step(mps_[w], bkp, left_[w], mpo), t2;
        
        t2 = transpose(t);
        gemm(t2, right_[w+1], t);
        return trace(t);
    }
    
    void normalize_upto(int w)
    {
        for (int i = 0; i < w; ++i)
        {
            block_matrix<Matrix, grp> t = mps_[i].normalize_left(SVD);
            mps_[i+1].multiply_from_left(t);
        }
        
        for (int i = L-1; i > w; --i)
        {
            block_matrix<Matrix, grp> t = mps_[i].normalize_right(SVD);
            mps_[i-1].multiply_from_right(t);
        }
    }
    
    int L;
    std::vector<MPSTensor<Matrix, grp> > mps_;
};

int main()
{
    int L = 20, M = 10;
    MPS mps(L, M);
    
    MPOTensor<Matrix, grp> id_mpo = identity_mpo<Matrix>(mps.mps_[0].site_dim());
    cout << mps.norm() << endl;
    mps.normalize_upto(3);
    cout << mps.norm() << endl;
    
    MPOTensor<Matrix, grp> sz_mpo = s12_sz_mpo<Matrix>(mps.mps_[0].site_dim());
    for (int k = 0; k < L; ++k) {
        cout << mps.stupid_site_expval(id_mpo, k) << " ";
        cout << mps.stupid_site_expval(sz_mpo, k) << endl;
    }
    
    cout << mps.norm(&sz_mpo) << endl;
}
