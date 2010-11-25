#include "mpstensor.h"

#include "general_matrix.hpp"
#include "matrix_interface.hpp"
#include "resizable_matrix_interface.hpp"
#include "matrix_algorithms.hpp"
typedef blas::general_matrix<double> Matrix;

#include "indexing.h"


int main()
{
    MPSTensor<Matrix, NullGroup> mps;
    mps.normalize_left();
    
    return 0;
}
