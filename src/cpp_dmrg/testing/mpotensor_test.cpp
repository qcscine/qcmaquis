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
    
    Index<grp> physical, aux;
    physical.insert(std::make_pair(NullGroup::Plus, 2));
    aux.insert(std::make_pair(NullGroup::Plus, 10));
    
    MPOTensor<Matrix, grp> mpo(physical, aux, aux);
    
    cout << mpo(std::make_pair(NullGroup::Plus, 0),
                std::make_pair(NullGroup::Plus, 0),
                std::make_pair(NullGroup::Plus, 0),
                std::make_pair(NullGroup::Plus, 0)) << endl;
    cout << mpo(std::make_pair(NullGroup::Plus, 0),
                std::make_pair(NullGroup::Plus, 0),
                std::make_pair(NullGroup::Plus, 0),
                std::make_pair(NullGroup::Plus, 1)) << endl;
    cout << mpo(std::make_pair(NullGroup::Plus, 1),
                std::make_pair(NullGroup::Plus, 0),
                std::make_pair(NullGroup::Plus, 0),
                std::make_pair(NullGroup::Plus, 1)) << endl;
}