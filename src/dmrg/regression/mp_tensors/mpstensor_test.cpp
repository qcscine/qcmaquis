#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

#include "dense_matrix/dense_matrix.hpp"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"
typedef blas::dense_matrix<double> Matrix;

#include "block_matrix/indexing.h"
#include "mp_tensors/mpstensor.h"
#include "mp_tensors/mpotensor.h"
#include "mp_tensors/contractions.h"

int main()
{
    typedef NullGroup grp;
    
    Index<grp> physical, aux1, aux2;
    physical.insert(std::make_pair(NullGroup::Plus, 2));
    aux1.insert(std::make_pair(NullGroup::Plus, 10));
    aux2.insert(std::make_pair(NullGroup::Plus, 11));
    
    MPSTensor<Matrix, grp> mps(physical, aux1, aux1);
    cout << mps << endl;
    cout << mps.scalar_norm() << endl;
    mps.multiply_by_scalar(1/mps.scalar_norm());
    cout << mps.scalar_norm() << endl;
    
    cout << mps.isleftnormalized() << endl;
    mps.normalize_left(SVD);
    cout << mps.isleftnormalized() << endl;
    cout << mps.scalar_norm() << endl;
    
    MPSTensor<Matrix, grp> mps2 = mps;
    block_matrix<Matrix, grp> left = identity_matrix<Matrix, grp>(mps.row_dim());
    left = contraction::overlap_left_step(mps, mps2, left);
    cout << left << endl;
    
    mps.normalize_right(SVD);
    mps2 = mps;
    block_matrix<Matrix, grp> right = identity_matrix<Matrix, grp>(mps.row_dim());
    right = contraction::overlap_right_step(mps, mps2, right);
    cout << right << endl;
    
    mps = MPSTensor<Matrix, grp>(physical, aux1, aux1);
    mps2 = MPSTensor<Matrix, grp>(physical, aux2, aux2);
    
    block_matrix<Matrix, grp> ovlp(aux2, aux1);
    contraction::overlap_left_step(mps2, mps, ovlp);
    contraction::overlap_left_step(mps2, mps, ovlp);
    contraction::overlap_right_step(mps2, mps, ovlp);
    contraction::overlap_right_step(mps2, mps, ovlp);
    
    return 0;
}
