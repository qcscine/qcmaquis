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

int main()
{
    typedef NullGroup grp;
    
    Index<grp> physical, aux;
    physical.insert(std::make_pair(NullGroup::Plus, 2));
    aux.insert(std::make_pair(NullGroup::Plus, 10));
    
    MPOTensor<Matrix, grp> mpo(physical, aux, aux);
    MPSTensor<Matrix, grp> mps(physical, aux, aux);
    
    block_matrix<Matrix, grp> left(aux*aux, aux);
    
    left = contraction::overlap_mpo_left_step(mps, mps, left, mpo);
    
}