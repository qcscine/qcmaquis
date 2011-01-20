#include <cmath>
#include <iterator>
#include <iostream>

using std::cout;
using std::cout;
using std::endl;

#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/matrix_interface.hpp"
#include "p_dense_matrix/resizable_matrix_interface.hpp"
#include "p_dense_matrix/dense_matrix_algorithms.h"
#include "p_dense_matrix/matrix_algorithms.hpp"
#include "p_dense_matrix/p_dense_matrix_blas.hpp"
typedef blas::p_dense_matrix<double> Matrix;

#include "block_matrix/indexing.h"
#include "mp_tensors/mps.h"
#include "mp_tensors/mpo.h"
#include "mp_tensors/contractions.h"
#include "mp_tensors/mps_mpo_ops.h"

#include "mp_tensors/special_mpos.h"

#include "mp_tensors/ss_optimize.h"

#include "mpos/adjancency.h"
#include "mpos/generate_mpo.h"

typedef U1 grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

int main(int argc, char ** argv)
{
    cout.precision(10);
    
    Index<grp> phys;
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(-1, 1));
    
//    Index<grp> phys;
//    phys.insert(std::make_pair(1, 1));
//    phys.insert(std::make_pair(0, 1));
//    phys.insert(std::make_pair(-1, 1));

    int L = atoi(argv[1]), W = atoi(argv[2]);    
//    double theta = atof(argv[3]);

    int N = L*W, M = 5;
    MPS<Matrix, grp> mps(N, M, phys);
    
    SquareAdj adj(L, W);
    Heisenberg<Matrix> H(1,1); // replace by factory, eventually
//    Spin1BlBq<Matrix> H(cos(theta*M_PI), sin(theta*M_PI));
    
    MPO<Matrix, grp> mpo = mpos::MPOMaker<Matrix, grp>::create_mpo(adj, H);
    
    cout << expval(mps, mpo, 0) << endl;
    cout << expval(mps, mpo, 1) << endl;
    
//    ss_optimize<Matrix, grp>(mps, mpo, 8, 1e-5, 1000);
}
