#include <cmath>
#include <iterator>
#include <iostream>

#include "alps/numeric/matrix/matrix.hpp"
#include "alps/numeric/matrix/matrix_interface.hpp"
#include "alps/numeric/matrix/resizable_matrix_interface.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/matrix/matrix_algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
typedef alps::numeric::matrix<double> Matrix;

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

#include "dmrg/mp_tensors/special_mpos.h"

#include "dmrg/mp_tensors/optimize.h"

typedef TrivialGroup grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

int main()
{
    Index<grp> phys; phys.insert(std::make_pair(TrivialGroup::Plus, 2));
    
    int L = 32, M = 10;
    MPS<Matrix, grp> mps(L, M, phys);
    
    MPOTensor<Matrix, grp> id_mpo = identity_mpo<Matrix>(mps[0].site_dim());
    MPOTensor<Matrix, grp> sz_mpo = s12_sz_mpo<Matrix>(mps[0].site_dim());
    
    MPO<Matrix, grp> szsz = s12_ising<Matrix>(L, -1, 1);
//    ss_optimize<Matrix, grp>(mps, MPO<Matrix, grp>(L, id_mpo));
    ss_optimize<Matrix, grp>(mps, szsz, 2, 1e-14, 1000);
}
