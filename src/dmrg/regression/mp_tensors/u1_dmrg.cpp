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

typedef U1 grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

int main()
{
    cout.precision(10);
    
    Index<grp> phys;
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(-1, 1));
    
    int L = 64, M = 2;
    MPS<Matrix, grp> mps(L, M, phys);
    
//    MPOTensor<Matrix, grp> id_mpo = identity_mpo<Matrix>(mps[0].site_dim());
//    MPOTensor<Matrix, grp> sz_mpo = s12_sz_mpo<Matrix>(mps[0].site_dim());
    
    MPO<Matrix, grp> szsz = s12_heisenberg<Matrix>(L, 1, 1);
    mps.normalize_left();
    cout << expval(mps, szsz, 0) << endl;
    cout << expval(mps, szsz, 1) << endl;
//    mps.normalize_right();
//    MPSTensor<Matrix, grp>::stupid_grow(mps[1], mps[2], 0, 0);
//    cout << expval(mps, szsz, 0) << endl;
//    cout << expval(mps, szsz, 1) << endl;
//    exit(0);
    
    ss_optimize<Matrix, grp>(mps, szsz, 10, 1e-14, 1000);
}
