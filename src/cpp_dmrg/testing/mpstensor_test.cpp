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
    left = overlap_left_step(mps, mps2, left);
    cout << left << endl;
    
    mps.normalize_right(SVD);
    mps2 = mps;
    block_matrix<Matrix, grp> right = identity_matrix<Matrix, grp>(mps.row_dim());
    right = overlap_right_step(mps, mps2, right);
    cout << right << endl;
    
    mps = MPSTensor<Matrix, grp>(physical, aux1, aux1);
    mps2 = MPSTensor<Matrix, grp>(physical, aux2, aux2);
    
    block_matrix<Matrix, grp> ovlp(aux2, aux1);
    overlap_left_step(mps2, mps, ovlp);
    overlap_left_step(mps2, mps, ovlp);
    overlap_right_step(mps2, mps, ovlp);
    overlap_right_step(mps2, mps, ovlp);
    
    return 0;
}
