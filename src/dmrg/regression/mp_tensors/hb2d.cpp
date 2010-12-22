#include <cmath>
#include <iterator>
#include <iostream>

using std::cout;
using std::cout;
using std::endl;

#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"
#include "dense_matrix/dense_matrix_blas.hpp"
typedef blas::dense_matrix<double> Matrix;

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

    int L = atoi(argv[1]), W = atoi(argv[2]);    

    int N = L*W, M = 2;
    MPS<Matrix, grp> mps(N, M, phys);
    
    SquareAdj adj(L, W);
    MPO<Matrix, grp> H =
    mpos::MPOMaker<Matrix, grp>::create_mpo(adj,
                                            mpos::MPOMaker<Matrix, grp>::hb_ops(1, 1));
    
    cout << expval(mps, H, 0) << endl;
    cout << expval(mps, H, 1) << endl;
    
    ss_optimize<Matrix, grp>(mps, H, 2, 1e-6, 1000);
}
