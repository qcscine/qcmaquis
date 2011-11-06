#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"
#include "dense_matrix/dense_matrix_blas.hpp"
#include "dense_matrix/aligned_allocator.h"

typedef blas::dense_matrix<double> Matrix;

#include <alps/hdf5.hpp>

#include "block_matrix/indexing.h"
#include "mp_tensors/mps.h"
#include "mp_tensors/mpo.h"
#include "mp_tensors/contractions.h"
#include "mp_tensors/mps_mpo_ops.h"
#include "mp_tensors/mpo_ops.h"
#include "mp_tensors/mps_initializers.h"


#include "app/hamiltonian.h"

#include "app/te_utils.hpp"
#include "mp_tensors/te.h"


using namespace app;

typedef U1 grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

std::ostream& operator<< (std::ostream& os, std::pair<grp::charge, std::size_t> const& p)
{
    os << "(" << p.first << " : " << p.second << ")";
    return os;
}

std::ostream& operator<< (std::ostream& os, index_product_iterator<grp>::value_type const& v)
{
    //std::copy(v.begin(), v.end(), std::ostream_iterator<std::pair<grp::charge, std::size_t> >(os, " "));
    for (int i=0; i<v.size(); ++i)
        os << v[i] << " ";
    return os;
}

std::ostream& operator<< (std::ostream& os, std::pair<MultiIndex<grp>::coord_t, MultiIndex<grp>::coord_t> const& p)
{
    os << p.first << ", " << p.second;
    return os;
}

std::ostream& operator<< (std::ostream& os, MPO<Matrix, grp> const& mpo)
{
    for (int p=0; p<mpo.size(); ++p)
        for (size_t r=0; r<mpo[p].row_dim(); ++r)
            for (size_t c=0; c<mpo[p].col_dim(); ++c)
                if (mpo[p].has(r, c))
                    os << "** Position " << p << " [" << r << "," << c << "]:" << endl << mpo[p](r,c);
    return os;
}

typedef Hamiltonian<Matrix, U1> ham;
typedef ham::op_t op_t;
typedef MultiIndex<grp>::index_id index_id;
typedef MultiIndex<grp>::set_id set_id;


MPO<Matrix, grp> block_to_mpo(Index<grp> const & phys_i,
                              block_matrix<Matrix, grp> ident,
                              block_matrix<Matrix, grp> bond_op,
                              std::size_t dist,
                              std::size_t length)
{
    
    MPO<Matrix, grp> mpo(length);
    
    std::size_t pos1 = 0, pos2 = dist;
    
    bond_op = reshape_2site_op(phys_i, bond_op);
    block_matrix<Matrix, grp> U, V, left, right;
    block_matrix<blas::associated_diagonal_matrix<Matrix>::type, grp> S, Ssqrt;
    svd(bond_op, U, V, S);
    Ssqrt = sqrt(S);
    gemm(U, Ssqrt, left);
    gemm(Ssqrt, V, right);
    
    // reshape and write back
    std::vector<block_matrix<Matrix, grp> > U_list = reshape_right_to_list(phys_i, left);
    std::vector<block_matrix<Matrix, grp> > V_list = reshape_left_to_list(phys_i, right);
    assert(U_list.size() == V_list.size());
    
    MPOTensor<Matrix, grp> left_tensor(1, U_list.size());
    MPOTensor<Matrix, grp> middle_tensor(U_list.size(), U_list.size());
    MPOTensor<Matrix, grp> right_tensor(U_list.size(), 1);
    
    for (std::size_t use_b=0; use_b<U_list.size(); ++use_b)
    {
        left_tensor(0, use_b) = U_list[use_b];
        middle_tensor(use_b, use_b) = ident;
        right_tensor(use_b, 0) = V_list[use_b];
    }
    mpo[pos1] = left_tensor;
    mpo[pos2] = right_tensor;
    for (std::size_t p=pos1+1; p<pos2; ++p)
    {
        mpo[p] = middle_tensor;
    }
    
    return mpo;
}


int main(int argc, char ** argv)
{
    
    Index<U1> phys;
    phys.insert(std::make_pair(0, 1));
    phys.insert(std::make_pair(1, 1));
    
    op_t ident;
    op_t create, destroy, count, sign;
    
    ident.insert_block(Matrix(1, 1, 1), 0, 0);
    ident.insert_block(Matrix(1, 1, 1), 1, 1);
    
    create.insert_block(Matrix(1, 1, 1), 0, 1);
    destroy.insert_block(Matrix(1, 1, 1), 1, 0);
    
    count.insert_block(Matrix(1, 1, 1), 1, 1);

    sign.insert_block(Matrix(1, 1, 1), 0, 0);
    sign.insert_block(Matrix(1, 1, -1), 1, 1);

    cout << "op 0:" << endl << create;
    cout << "op 1:" << endl << destroy;
    
    size_t dist = 2;
    
    MultiIndex<grp> midx;
    std::vector<std::pair<index_id, bool> > left_vec, right_vec;
    for (size_t p=0; p<=dist; ++p) {
        index_id id1 = midx.insert_index(phys);
        index_id id2 = midx.insert_index(phys);
        
        left_vec.push_back( std::make_pair(id1, true) );
        right_vec.push_back( std::make_pair(id2, true) );
    }
    set_id op_set = midx.create_set(left_vec, right_vec);
    
    op_t bond_op, tmp;
    op_kron_long(midx, op_set, create, destroy, ident, dist, bond_op);
    
    op_kron_long(midx, op_set, destroy, create, ident, dist, tmp);
    bond_op += tmp;
//    op_kron_long(midx, op_set, count, ident, ident, dist, tmp);
//    bond_op += tmp;
//    op_kron_long(midx, op_set, ident, count, ident, dist, tmp);
//    bond_op += tmp;
    
    op_t bond_op2;
    op_kron(phys, create, destroy, bond_op2);
    
    op_kron(phys, destroy, create, tmp);
    bond_op2 += tmp;
//    op_kron(phys, count, ident, tmp);
//    bond_op2 += tmp;
//    op_kron(phys, ident, count, tmp);
//    bond_op2 += tmp;
    
    
    MPO<Matrix, grp> block_mpo_easy = block_to_mpo(phys, ident, bond_op2, dist, dist+1);
    cout << "MPO (easy):" << endl << block_mpo_easy;
    
    MPO<Matrix, grp> block_mpo = block_to_mpo(phys, bond_op, dist+1);
    cout << "MPO:" << endl << block_mpo;
    
}
