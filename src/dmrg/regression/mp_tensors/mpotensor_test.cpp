#include <iostream>

using std::cout;
using std::cerr;
using std::endl;
#include "utils/zout.hpp"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/contractions.h"

#include "dmrg/mp_tensors/special_mpos.h"

int main()
{
    typedef TrivialGroup grp;
    
    Index<grp> physical, aux, mpo_aux;
    physical.insert(std::make_pair(TrivialGroup::Plus, 2));
    aux.insert(std::make_pair(TrivialGroup::Plus, 10));
    mpo_aux.insert(std::make_pair(TrivialGroup::Plus, 1));
    
    MPOTensor<Matrix, grp> mpo(physical, mpo_aux, mpo_aux);
    MPSTensor<Matrix, grp> mps(physical, aux, aux);
    
    block_matrix<Matrix, grp> left(aux*aux, aux), right(aux*aux, aux), o;
    
    left = contraction::overlap_mpo_left_step(mps, mps, left, mpo);
    left = contraction::overlap_mpo_left_step(mps, mps, left, mpo);
    
    right = contraction::overlap_mpo_right_step(mps, mps, right, mpo);
    right = contraction::overlap_mpo_right_step(mps, mps, right, mpo);
    
    right = transpose(right);
    gemm(left, right, o);
    zout << trace(o) << endl;
}
