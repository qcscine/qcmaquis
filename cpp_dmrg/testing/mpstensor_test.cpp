#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

#include "mpstensor.h"

#include "general_matrix.hpp"
#include "matrix_interface.hpp"
#include "resizable_matrix_interface.hpp"
#include "matrix_algorithms.hpp"
typedef blas::general_matrix<double> Matrix;

#include "indexing.h"


int main()
{
    typedef NullGroup grp;
    
    Index<grp> physical, aux;
    physical.insert(std::make_pair(NullGroup::Plus, 2));
    aux.insert(std::make_pair(NullGroup::Plus, 10));
    
    MPSTensor<Matrix, grp> mps(physical, aux, aux);
    mps.normalize_left();
    mps.multiply_by_scalar(2);
    
    MPSTensor<Matrix, grp> mps2 = mps;
    
    cout << mps.scalar_norm() << endl;
//    cout << mps.scalar_overlap(mps2) << endl;
    
    return 0;
}
