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


int main()
{
    typedef NullGroup grp;
    
    Index<grp> physical, aux;
    physical.insert(std::make_pair(NullGroup::Plus, 2));
    aux.insert(std::make_pair(NullGroup::Plus, 10));
    
    MPSTensor<Matrix, grp> mps(physical, aux, aux);
    cout << mps << endl;
    cout << mps.scalar_norm() << endl;
    mps.multiply_by_scalar(1/mps.scalar_norm());
    cout << mps.scalar_norm() << endl;
    
    cout << mps.isleftnormalized() << endl;
    mps.normalize_left(SVD);
    mps.normalize_right(SVD);
    cout << mps.isleftnormalized() << endl;
    cout << mps.scalar_norm() << endl;
    
    MPSTensor<Matrix, grp> mps2 = mps;
    block_matrix<Matrix, grp> left = identity_matrix<Matrix, grp>(mps.row_dim());
    left = overlap_left_step(mps, mps2, left);
    cout << left << endl;
    
    return 0;
}
