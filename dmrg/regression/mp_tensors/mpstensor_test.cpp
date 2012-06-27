#include <iostream>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "dmrg/kernels/alps_matrix.hpp"
typedef alps::numeric::matrix<double> Matrix;


#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/special_mpos.h"

void ng()
{
    typedef TrivialGroup grp;
    
    Index<grp> physical, aux1, aux2;
    physical.insert(std::make_pair(TrivialGroup::Plus, 2));
    aux1.insert(std::make_pair(TrivialGroup::Plus, 10));
    aux2.insert(std::make_pair(TrivialGroup::Plus, 11));
    
    MPSTensor<Matrix, grp> mps(physical, aux1, aux1);
    maquis::cout << mps.scalar_norm() << std::endl;
    mps.multiply_by_scalar(1/mps.scalar_norm());
    maquis::cout << mps.scalar_norm() << std::endl;
    
    mps.normalize_left(SVD);
    maquis::cout << "Norm after normalization: " << mps.scalar_norm() << std::endl;
    
    MPSTensor<Matrix, grp> mps2 = mps;
    block_matrix<Matrix, grp> left = identity_matrix<Matrix, grp>(mps.row_dim());
    left = contraction::overlap_left_step(mps, mps2, left);
    maquis::cout << left << std::endl;
    
    mps.normalize_right(SVD);
    mps2 = mps;
    block_matrix<Matrix, grp> right = identity_matrix<Matrix, grp>(mps.row_dim());
    right = contraction::overlap_right_step(mps, mps2, right);
    maquis::cout << right << std::endl;
    
    mps = MPSTensor<Matrix, grp>(physical, aux1, aux1);
    mps2 = MPSTensor<Matrix, grp>(physical, aux2, aux2);
    
    block_matrix<Matrix, grp> ovlp(aux2, aux1);
    contraction::overlap_left_step(mps2, mps, ovlp);
    contraction::overlap_left_step(mps2, mps, ovlp);
    contraction::overlap_right_step(mps2, mps, ovlp);
    contraction::overlap_right_step(mps2, mps, ovlp);
}

void u1()
{
    typedef U1 grp;
    
    Index<grp> physical, aux1, aux2;
    physical.insert(std::make_pair(-1, 1));
    physical.insert(std::make_pair(1, 1));
    
    aux1.insert(std::make_pair(1, 3));
    aux1.insert(std::make_pair(-1, 3));
    
    aux2.insert(std::make_pair(2, 3));
    aux2.insert(std::make_pair(0, 3));
    aux2.insert(std::make_pair(-2, 3));
    
    MPSTensor<Matrix, grp> mps(physical, aux1, aux2);
    
    maquis::cout << mps.scalar_norm() << std::endl;
    mps.multiply_by_scalar(1/mps.scalar_norm());
    maquis::cout << mps.scalar_norm() << std::endl;
    mps.normalize_left(SVD);
    maquis::cout << mps.scalar_norm() << std::endl;
    
//    mps.make_left_paired();
//    maquis::cout << mps << std::endl;
//    
//    mps.make_right_paired();
//    mps.make_left_paired();
//    
//    maquis::cout << mps << std::endl;
    
    MPSTensor<Matrix, grp> mps2 = mps;
//    block_matrix<Matrix, grp> left = identity_matrix<Matrix, grp>(mps.row_dim());
//    left = contraction::overlap_left_step(mps, mps2, left);
//    maquis::cout << left << std::endl;
//    mps.normalize_right(SVD);
//    mps2 = mps;
//    left = contraction::overlap_right_step(mps, mps2, left);
//    maquis::cout << left << std::endl;
    
    MPOTensor<Matrix, U1> ident(1, 1);
    ident(0,0) = block_matrix<Matrix, U1>();
    ident(0,0).insert_block(Matrix(1, 1, 1), -1, -1);
    ident(0,0).insert_block(Matrix(1, 1, 1), 1, 1);
    
    MPOTensor<Matrix, U1> splus(1, 1);
    splus(0,0) = block_matrix<Matrix, U1>();
    splus(0,0).insert_block(Matrix(1, 1, 1), 1, -1);
    splus(0,0).insert_block(Matrix(1, 1, 1), -1, 1);
    
    mps.normalize_left(SVD);
    mps2 = mps;
    Boundary<Matrix, grp> mleft(mps.row_dim(), mps.row_dim(), 1);
    mleft.data_[0] = identity_matrix<Matrix, grp>(mps.row_dim());
//    maquis::cout << mleft.data_[0] << std::endl;
//    mleft = contraction::overlap_mpo_left_step(mps, mps2, mleft, splus);
//    maquis::cout << mleft.data_[0] << std::endl;
    
    mps.normalize_right(SVD);
    mps2 = mps;
    
    mleft.data_[0] = identity_matrix<Matrix, grp>(mps.col_dim());
    mleft.data_[0] = contraction::overlap_right_step(mps, mps2, mleft.data_[0]);
    maquis::cout << mleft.data_[0] << std::endl;
    
    mleft.data_[0] = identity_matrix<Matrix, grp>(mps.col_dim());
    mleft = contraction::overlap_mpo_right_step(mps, mps2, mleft, ident);
    maquis::cout << mleft.data_[0] << std::endl;
}

int main()
{
    ng();
    u1();
}

