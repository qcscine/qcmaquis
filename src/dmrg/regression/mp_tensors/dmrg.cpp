#include <cmath>
#include <iterator>
#include <iostream>

using std::cout;
using std::cerr;
using std::endl;

#include "dense_matrix.hpp"
#include "matrix_interface.hpp"
#include "resizable_matrix_interface.hpp"
#include "dense_matrix_algorithms.h"
#include "matrix_algorithms.hpp"
typedef blas::dense_matrix<double> Matrix;

#include "indexing.h"
#include "mps.h"
#include "mpo.h"
#include "contractions.h"
#include "mps_mpo_ops.h"

#include "special_mpos.h"

#include "ss_optimize.h"

typedef NullGroup grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

int main()
{
    Index<grp> phys; phys.insert(std::make_pair(NullGroup::Plus, 2));
    
    int L = 32, M = 100;
    MPS<Matrix, grp> mps(L, M, phys);
    
    MPOTensor<Matrix, grp> id_mpo = identity_mpo<Matrix>(mps[0].site_dim());
    MPOTensor<Matrix, grp> sz_mpo = s12_sz_mpo<Matrix>(mps[0].site_dim());
    
    MPO<Matrix, grp> szsz = s12_ising<Matrix>(L, -1, 1);
//    ss_optimize<Matrix, grp>(mps, MPO<Matrix, grp>(L, id_mpo));
    ss_optimize<Matrix, grp>(mps, szsz, 2);
}
